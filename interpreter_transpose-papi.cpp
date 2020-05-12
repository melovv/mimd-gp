/*
 * tinygp, a minimalist genetic programming system in C++ ( (c) M.Keijzer 2004. License: GPL-2 or later)
 */
#include <sys/auxv.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <string>
#include <chrono>
#include <time.h>
#include <mm_malloc.h>

#include "papi.h"


using namespace std;
using namespace std::chrono;


// Parameters: floats are captured in an unsigned: FP_BASE is the denominator. For more tight control, increase FP_BASE 
//enum {MAX_LEN=1000000, POPSIZE=10, GENERATIONS=50};
enum {MAX_LEN=1000000};
int DEPTH=5, POPSIZE=100, GENERATIONS=1;

// Function (the only typedef) is the thing used in evaluation (eval_var, eval_plus, eval_mult, eval_div are all Functions)
typedef void (*Function)();

/* Global data */
Function			jumptable[58];		    // global jumptable (ascii(9) == 57, accomodating 10 variables)

string::const_iterator		current_node;		    // for prefix jumptable evaluation

float ** input_data;
float ** stack;

float *linear_input, *linear_stack;    


int top=0;
unsigned    nvars=9, nprimitives=0, cont=0;

long int ncases;
long int GPop=0;
float eps = 1e-8;

/* General utility functions */
template <typename T> T random(T mx)	         { return static_cast<T>(1./(1.+RAND_MAX) * rand() * mx); }
inline int arity_min_1(char x, char dummy = '0') { return x < '0'?1:-1; }   // ascii '+','*','-','/' are all smaller than '0' (see tokens below)

/* Evaluation using the prefix-jumptable mechanism (see Mike Keith and Martin
 * C. Martin, Genetic Programming in C++, AIGP-1.) To get rid of the warning
 * (with -Wall) about using char as a subscript, needs a static_cast.  char as
 * subscript is well defined, so I'm leaving it as such. 
 *
 * Evaluation is vectorized, so we're evaluating all cases in one go. For
 * non-trivial numbers of fitness cases (> 100), this effectively eliminates
 * the overhead from the interpreter*/

float div_mask (float val) {
  return (std::abs(val) > eps);
}

/*
  void eval_var()	    { 
  int var = *(current_node+1)-'0';
  top++;
  
  #pragma  ivdep
  #pragma vector always   
  for (unsigned j = 0; j < ncases; ++j) 
  stack[top][j] = input_data[var][j]; 
  GPop++;
  }
  
  //void eval_plus()	{ stack[top-1] = stack[top] + stack[top-1]; top--; }
  
  void eval_plus()    { 
  #pragma  ivdep
  #pragma vector always
  for (unsigned j = 0; j < ncases; ++j) 
  stack[top-1][j] += stack[top][j];
  
  top--;  
  GPop++;
  }
  
  void eval_mult()	{ 
  #pragma  ivdep
  #pragma vector always
  for (unsigned j = 0; j < ncases; ++j) 
  stack[top-1][j] *= stack[top][j];
  top--;  GPop++;}
  
  void eval_minus()	{ 
  #pragma  ivdep
  #pragma vector always
  for (unsigned j = 0; j < ncases; ++j) 
  stack[top-1][j] -= stack[top][j];    
  top--;  GPop++;}
  
  
  void eval_div() {
  #pragma  ivdep
  #pragma vector always
  for (unsigned i = 0; i < ncases; ++i) stack[top-1][i] = std::abs(stack[top][i]) > eps ? stack[top-1][i] / stack[top][i]: 0.f;
  top--;  GPop++;
  }
*/


// ROW-MAJOR
void eval_var()	    { 
  //    __assume_aligned(stack,32);
  //    __assume_aligned(input_data,32);   
  int var = *(current_node+1)-'0';
  top++;
  
  float * __restrict a = stack[top];
  float * __restrict b = input_data[var];    
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  
#pragma  ivdep
#pragma vector always   
  for (int i = 0; i < ncases; ++i) 
    a[i] = b[i]; 
  GPop++;
}

void eval_plus()    { 
  float * __restrict a = stack[top-1];
  float * __restrict b = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  
#pragma ivdep
#pragma vector always    
  for (int i = 0; i < ncases; ++i) {
    a[i] += b[i];
  }
  
  top--;  
  GPop++;
}

void eval_mult()    { 
  float * __restrict a = stack[top-1];
  float * __restrict b = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  
#pragma ivdep
#pragma vector always    
  for (int i = 0; i < ncases; ++i) {
    a[i] *= b[i];
  }
  
  top--;  
  GPop++;
}

void eval_minus()    { 
  float * __restrict a = stack[top-1];
  float * __restrict b = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  
#pragma ivdep
#pragma vector always    
  for (int i = 0; i < ncases; ++i) {
    a[i] -= b[i];
  }
  
  top--;  
  GPop++;
}


void eval_div()    { 
  float * __restrict a = stack[top-1];
  float * __restrict b = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  
#pragma ivdep
#pragma vector always    
  for (int i = 0; i < ncases; ++i) {
    a[i] = std::abs(b[i]) > eps ? a[i] / b[i]: 0.f;
  }
  
  top--;  
  GPop++;
}



/* Two more globals, language and functions */
const char			tokens[]	        = "+*-/";
//const char			tokens[]	        = "+**+";
const Function		binary_functions[4] = { eval_plus, eval_mult, eval_minus, eval_div };	// for initializing the jumptable
//const Function		binary_functions[4] = { eval_plus, eval_mult, eval_minus, eval_plus };	// for initializing the jumptable
//const Function		binary_functions[4] = { eval_plus , eval_mult, eval_mult, eval_plus };


/* calculates the fitness of an individual as the sum of absolute errors */
float eval(const string& t) {
  //nprimitives += t.size(); // we will evaluate all nodes, making this a correct count
  current_node = --t.end();
  cont++;
  top = -1;
  for( unsigned int a = t.size(); a>0; a--)
    jumptable[*current_node--]();
  
  float val=0;
  #pragma ivdep
  #pragma vector always      
  for (unsigned j = 0; j < ncases; ++j) 
    val += std::abs(stack[top][j] - input_data[nvars][j]);
  // cout << val;
  return val; 
}

/* Initialization of an individual, the string(int,char) ctor intializes the
 * string with int copies of the char argument */
string build_tree(int grow=1, unsigned depth_left=DEPTH) {
  if(0==depth_left || (grow and random(2))) return string(1,'0' + random<char>(nvars));
  return string(1,tokens[random(4)]) + build_tree(grow, depth_left-1) + build_tree(grow, depth_left-1); 
}

int main(int argc, char* argv[]) {
  // get command line params, no error checking, so beware
  ncases = atoi(argv[1]);
  //nvars = atoi(argv[2]);
  POPSIZE = atoi(argv[2]);
  //GENERATIONS = atoi(argv[2]);
  GENERATIONS = 1;
  DEPTH = atoi(argv[3]);
  GPop = 0;
  
  //POPSIZE = 100; //100;
  
  //srand(time(NULL));
  //srand(1);
  unsigned int *seed;
  
  seed = (unsigned int *)getauxval(AT_RANDOM);
  srand(*seed);
  
  //srand(ncases + GENERATIONS + DEPTH);
  //srand(1);
  
  //PAPI temp vars
  int EventSet = PAPI_NULL;
  int retval;
  long long values[4];
  
  // init jumptable for all possible variables and functions
  fill_n(jumptable+'0',10,&eval_var);
  for (unsigned i = 0; i < 4; ++i) jumptable[tokens[i]] = binary_functions[i];
  //jumptable[tokens[0]] = binary_functions[0];
  
  // load data into the globals input_data and targets, does not check validity of the stream
  linear_input = (float*) _mm_malloc (sizeof (float) * ncases * (nvars+1), 32);
  input_data = (float**) _mm_malloc (sizeof (float *) * (nvars+1), 32);
  for (unsigned i = 0; i < (nvars+1); i++)
    {
      input_data[i] = &linear_input[i*(ncases)];
      for (unsigned j = 0; j < ncases; j++)
	input_data[i][j] =  -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-(-1.0))));
    }    

  linear_stack = (float*) _mm_malloc (sizeof (float) * (DEPTH*2) * ncases, 32);
  stack = (float**) _mm_malloc (sizeof (float *) * (DEPTH*2), 32);
  for (unsigned i = 0; i < (DEPTH*2); i++)
    stack[i] = &linear_stack[i*(ncases)];
  
  /*
    input_data = (float**) malloc (sizeof(float*) * (nvars+1));
    stack = (float**) malloc (sizeof(float*) * DEPTH*2);
    for (unsigned i = 0; i < nvars+1; ++i) {
    input_data[i] = (float*) malloc (sizeof(float) * ncases);
    for (unsigned j = 0; j < ncases; ++j)
    input_data[i][j] =  -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0-(-1.0))));; 
    }
    
    for (unsigned i = 0; i < DEPTH*2; ++i)
    stack[i] = (float*) malloc (sizeof(float) * ncases);    
  */
  
  // init population, note that genotype and fitness are decoupled.
  vector<string>  trees(POPSIZE);
  valarray<float> fitness(POPSIZE);
  
  /* Init PAPI library */
  retval = PAPI_library_init( PAPI_VER_CURRENT );
  if ( retval != PAPI_VER_CURRENT ) {
    printf("Erro em PAPI_library_init : retval = %d\n", retval);
    exit (1);
  }
  
  if ( ( retval = PAPI_create_eventset( &EventSet ) ) != PAPI_OK ) {
    printf("Erro em PAPI_create_eventset : retval = %d\n", retval);
    exit (1);
  }
  
  if (PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK) {
    printf("Erro em PAPI_L1_TCM\n");
    exit (1);
  }
  
  if (PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK) {
    printf("Erro em PAPI_L2_TCM\n");
    exit (1);
  }
  
  if (PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK) {
    printf("Erro em PAPI_L3_TCM\n");
    exit (1);
  }

  for (int i = 0; i < POPSIZE; ++i) {
    trees[i] = build_tree(1, DEPTH);
    //cout << trees[i].c_str() << std::endl;
 }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();    

  if ( ( retval = PAPI_start( EventSet ) ) != PAPI_OK ) {
    printf("Erro em PAPI_start"); exit (1); }
  //srand(1);
  for (unsigned iter = 0; iter < GENERATIONS; ++iter) {
    //for (unsigned iter = 0; iter < 1000; ++iter) {

    for (int i = 0; i < POPSIZE; ++i) {

      //GPop = 0;	
      // fitness[i] = eval(trees[i] = build_tree()); 
      //fitness[i] = eval(trees[i] = build_tree(1, DEPTH)); 
      //trees[i] = build_tree(1, DEPTH);
      //cout << trees[i].c_str() << "  Size=" << trees[i].size() << std::endl;
      //cout << trees[i].c_str() << "  " << GPop << std::endl;

      fitness[i] = eval(trees[i]);
      
      //cout << fitness[i];
      //cout << trees[i].c_str() << "  Size=" << trees[i].size() << "  " << GPop << std::endl;
    }
  }
  if ( ( retval = PAPI_stop( EventSet, NULL ) ) != PAPI_OK ) {
    printf("Erro em PAPI_stop"); exit (1); }
  if ( ( retval = PAPI_read( EventSet, &values[0] ) ) != PAPI_OK ) {
    printf("Erro em PAPI_read"); exit (1); }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  
  //std::cout << "Time required for " << cont << " evaluations was " <<
  //    std::chrono::duration_cast<std::chrono::milliseconds>((t2 - t1)).count() << 
  //    " milliseconds." << std::endl;    
  
  //printf("\n\ncont=%d\n", cont);
  
  //cout << trees[0].c_str() << std::endl;
  //cout << "GPop = " << GPop << std::endl;
  cout << fixed << "SIMD, " << ncases << ", " << POPSIZE*GENERATIONS << ", "
       << GPop*ncases << ", "
       << (std::chrono::duration_cast<std::chrono::microseconds>((t2 - t1)).count()/1000000.0)
       << ", " << GPop*ncases/(float)(std::chrono::duration_cast<std::chrono::microseconds>((t2 - t1)).count()/1000000.0)
       << ", " << values[0] << ", " << values[1] << ", " << values[2]
       << ", " << (values[0]+values[1]+values[2])
       << ", " << (double)(values[0]+values[1]+values[2])/(double)(GPop*ncases)
       << std::endl;
  
  _mm_free (linear_input);
  _mm_free (input_data);
  
  _mm_free (linear_stack);
  _mm_free(stack);        
}


