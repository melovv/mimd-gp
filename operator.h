/*
 * operators as proposed in paper XXXX YYY
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <string>
#include <chrono>

using namespace std;

// Function (the only typedef) is the thing used in evaluation (eval_var, eval_plus, eval_mult, eval_div are all Functions)
typedef void (*Function)();

/* Global data */
Function			jumptable[255];		    // global jumptable (ascii(9) == 57, accomodating 10 variables)

string::const_iterator		current_node;		    // for prefix jumptable evaluation

float ** input_data;
float ** stack;

float *linear_input, *linear_stack;    


long int ncases;

int top;
long int GPop;
float eps;

// ROW-MAJOR
void xxxeval_var()	    {top++; GPop++;}


void eval_var()	    { 
//    __assume_aligned(stack,32);
//    __assume_aligned(input_data,32);   
    int var = *(current_node+1)-'0';
    top++;

	//cout << "top=" << top << std::endl;
    float * __restrict a = stack[top];
float * __restrict b = input_data[var];    
    __assume_aligned(a,32);
    __assume_aligned(b,32);

	#pragma  ivdep
	#pragma vector always   
#pragma simd
    for (long int i = 0; i < ncases; ++i) 
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
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
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
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
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
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
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
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] = std::abs(b[i]) > eps ? a[i] / b[i]: 0.f;
    }

    top--;
    GPop++;
}


// Tertiary functions, consuming 5 instructions
void eval_plus_plus()    { 
    float * __restrict a = stack[top-2];
    float * __restrict b = stack[top-1];
    float * __restrict c = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] += b[i] + c[i];
    }

    top-=2;
    GPop+=2;
}

void eval_plus_mult()    { 
    float * __restrict a = stack[top-2];
    float * __restrict b = stack[top-1];
    float * __restrict c = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] = (a[i]+b[i]) * c[i];
    }
    top-=2;
    GPop+=2;
}

void eval_plus_minus()    { 
    float * __restrict a = stack[top-2];
    float * __restrict b = stack[top-1];
    float * __restrict c = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] += b[i] - c[i];
    }

    top-=2;
    GPop+=2;
}


void eval_plus_div()    { 
    float * __restrict a = stack[top-2];
    float * __restrict b = stack[top-1];
    float * __restrict c = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] = std::abs(c[i]) > eps ? (a[i] + b[i]) / c[i]: 0.f;
    }

    top-=2;
    GPop+=2;
}


void eval_minus_plus()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] -= b[i] + c[i];
  }
  
  top-=2;
  GPop+=2;
}

void eval_minus_mult()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] = (a[i]-b[i]) * c[i];
  }
  top-=2;
  GPop+=2;
}

void eval_minus_minus()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] -= b[i] - c[i];
  }
  
  top-=2;
  GPop+=2;
}


void eval_minus_div()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] = std::abs(c[i]) > eps ? (a[i] - b[i]) / c[i]: 0.f;
  }
  
  top-=2;
  GPop+=2;
}




void eval_mult_plus()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] = (a[i]*b[i]) + c[i];
  }
  
  top-=2;
  GPop+=2;
}

void eval_mult_mult()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] *= b[i] * c[i];
  }
  top-=2;
  GPop+=2;
}

void eval_mult_minus()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
#pragma ivdep
#pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) {
    a[i] = (a[i]*b[i]) - c[i];
  }
  
  top-=2;
  GPop+=2;
}


void eval_mult_div()    { 
  float * __restrict a = stack[top-2];
  float * __restrict b = stack[top-1];
  float * __restrict c = stack[top];
  __assume_aligned(a,32);
  __assume_aligned(b,32);
  __assume_aligned(c,32);
  
  #pragma ivdep
  #pragma vector always    
#pragma simd
  for (long int i = 0; i < ncases; ++i) 
    a[i] = std::abs(c[i]) > eps ? (a[i] * b[i]) / c[i]: 0.f;
  
  
  top-=2;
  GPop+=2;
}


 void eval_div_plus() {
   float * __restrict a = stack[top-2];
   float * __restrict b = stack[top-1];
   float * __restrict c = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   
	#pragma  ivdep
	#pragma vector always
#pragma simd
	for (unsigned i = 0; i < ncases; ++i) 
	  a[i] = std::abs(b[i]) > eps ? (a[i] / b[i]) + c[i]: c[i];
	
	top-=2;
	GPop+=2;
}
 void eval_div_minus() {
   float * __restrict a = stack[top-2];
   float * __restrict b = stack[top-1];
   float * __restrict c = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   
   #pragma  ivdep
   #pragma vector always
#pragma simd
   for (unsigned i = 0; i < ncases; ++i) 
     a[i] = std::abs(b[i]) > eps ? (a[i] / b[i]) - c[i]: c[i];
   
	top-=2;
	GPop+=2; 
}
 void eval_div_mult() {
   float * __restrict a = stack[top-2];
   float * __restrict b = stack[top-1];
   float * __restrict c = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   
   #pragma  ivdep
   #pragma vector always
#pragma simd
   for (unsigned i = 0; i < ncases; ++i) 
     a[i] = std::abs(b[i]) > eps ? (a[i] / b[i]) * c[i]: c[i];
   
	top-=2;
	GPop+=2; 
}
 void eval_div_div() {
   float * __restrict a = stack[top-2];
   float * __restrict b = stack[top-1];
   float * __restrict c = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   
  #pragma  ivdep
  #pragma vector always
#pragma simd
	for (unsigned i = 0; i < ncases; ++i) {
	  a[i] = std::abs(b[i]) > eps ? (a[i] / b[i]): 0.f;
	  a[i] = std::abs(c[i]) > eps ? (a[i] / c[i]): 0.f;
	}
	top-=2;
	GPop+=2; 
}


// Quartiary functions, consuming 7 instructions
void eval_plus_plus_plus()    { 
    float * __restrict a = stack[top-3];
    float * __restrict b = stack[top-2];
    float * __restrict c = stack[top-1];
    float * __restrict d = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);
    __assume_aligned(d,32);

    #pragma ivdep
    #pragma vector always    
    #pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] += b[i] + c[i] + d[i];
    }

    top-=3;
    GPop+=3;
}

void eval_plus_plus_mult()    { 
    float * __restrict a = stack[top-3];
    float * __restrict b = stack[top-2];
    float * __restrict c = stack[top-1];
    float * __restrict d = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);
    __assume_aligned(d,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] = (a[i]+b[i]) * (c[i]+d[i]);
    }
    top-=3;
    GPop+=3;
}

void eval_plus_plus_minus()    { 
    float * __restrict a = stack[top-3];
    float * __restrict b = stack[top-2];
    float * __restrict c = stack[top-1];
    float * __restrict d = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);
    __assume_aligned(d,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        a[i] += b[i] - c[i] + d[i];
    }

    top-=3;
    GPop+=3;
}


void eval_plus_plus_div()    { 
    float * __restrict a = stack[top-3];
    float * __restrict b = stack[top-2];
    float * __restrict c = stack[top-1];
    float * __restrict d = stack[top];
    __assume_aligned(a,32);
    __assume_aligned(b,32);
    __assume_aligned(c,32);
    __assume_aligned(d,32);

    #pragma ivdep
    #pragma vector always    
#pragma simd
    for (long int i = 0; i < ncases; ++i) {
        c[i] += d[i];
        a[i] = std::abs(c[i]) > eps ? (a[i] + b[i]) / c[i]: 0.f;
    }

    top-=3;
    GPop+=3;
}



 void eval_plus_minus_mult()    { 
   float * __restrict a = stack[top-3];
   float * __restrict b = stack[top-2];
   float * __restrict c = stack[top-1];
   float * __restrict d = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   __assume_aligned(d,32);
   
   #pragma ivdep
   #pragma vector always    
#pragma simd
   for (long int i = 0; i < ncases; ++i) 
     a[i] = (a[i]+b[i]) * (c[i] - d[i]);

	top-=3;  GPop+=3; 
}

 void eval_plus_mult_minus()    { 
   float * __restrict a = stack[top-3];
   float * __restrict b = stack[top-2];
   float * __restrict c = stack[top-1];
   float * __restrict d = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   __assume_aligned(d,32);
   
   #pragma ivdep
   #pragma vector always    
#pragma simd
   for (long int i = 0; i < ncases; ++i) 
     a[i] = (a[i]+b[i]) - (c[i] * d[i]);

	top-=3;  GPop+=3; 
}			 
	 

void eval_plus_mult_div() {
   float * __restrict a = stack[top-3];
   float * __restrict b = stack[top-2];
   float * __restrict c = stack[top-1];
   float * __restrict d = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   __assume_aligned(d,32);
   
   #pragma ivdep
   #pragma vector always    
#pragma simd
   for (long int i = 0; i < ncases; ++i) {
     c[i] *= d[i];
     a[i] = std::abs(c[i]) > eps ? (a[i] + b[i]) / c[i]: 0.f;
   }

	top-=3;
	GPop+=3;
 }


 void eval_plus_div_minus() {
   float * __restrict a = stack[top-3];
   float * __restrict b = stack[top-2];
   float * __restrict c = stack[top-1];
   float * __restrict d = stack[top];
   __assume_aligned(a,32);
   __assume_aligned(b,32);
   __assume_aligned(c,32);
   __assume_aligned(d,32);
   
   #pragma ivdep
   #pragma vector always    
#pragma simd
   for (long int i = 0; i < ncases; ++i) {
     a[i] = std::abs(d[i]) > eps ? (a[i] + b[i]) - (c[i]/d[i]): 0.f;
   }

	top-=3;
	GPop+=3;
}




/* Two more globals, language and functions */
//const int           qtd_tokens = 24;
//*
const int           qtd_tokens = 28;
//const int           qtd_tokens = 4;
const char			tokens[]      = "+*-/!@#$%^&;()_={}[].<>?QWERTY";
const int           arities[qtd_tokens]       = {2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4};
//const int           arities[qtd_tokens]       = {2, 2, 2, 2, 4, 4, 4, 4}; 
const Function		all_functions[qtd_tokens] = { 
                                            eval_plus, eval_mult, eval_minus, eval_div,     // for initializing the jumptable with binary functions
                                            eval_plus_plus, eval_plus_minus, eval_plus_mult, eval_plus_div,     // for initializing the jumptable with ternary functions
                                            eval_minus_plus, eval_minus_minus, eval_minus_mult, eval_minus_div,
                                            eval_mult_plus, eval_mult_minus, eval_mult_mult, eval_mult_div,
                                            eval_div_plus, eval_div_minus, eval_div_mult, eval_div_div,
                                            eval_plus_plus_plus, eval_plus_plus_minus, eval_plus_plus_mult, eval_plus_plus_div,      // for initializing the jumptable with quaternary functions
                                            eval_plus_minus_mult, eval_plus_mult_minus, eval_plus_mult_div, eval_plus_div_minus      // for initializing the jumptable with quaternary functions
                                        };

// */


/*
const int           qtd_tokens = 16;
const char			tokens[]      = "+*-/!@#$%^&;()_={}[].<>?QWERTY";
//const int           arities[qtd_tokens]       = {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};
const int           arities[qtd_tokens]       = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4}; 
const Function		all_functions[qtd_tokens] = { 
                                            eval_plus, eval_mult, eval_minus, eval_div,     // for initializing the jumptable with binary functions
                                            eval_plus_plus, eval_plus_minus, eval_plus_mult, eval_plus_div,     // for initializing the jumptable with ternary functions
//                                            eval_minus_plus, eval_minus_minus, eval_minus_mult,// eval_minus_div,
//                                            eval_mult_plus, eval_mult_minus, eval_mult_mult, //eval_mult_div,
//                                            eval_div_plus, eval_div_minus, eval_div_mult, eval_div_div,
                                            eval_plus_plus_plus, eval_plus_plus_minus, eval_plus_plus_mult, eval_plus_plus_div,      // for initializing the jumptable with quaternary functions
                                            eval_plus_minus_mult, eval_plus_mult_minus, eval_plus_mult_div, eval_plus_div_minus      // for initializing the jumptable with quaternary functions
                                        };


// */
/*
const int           qtd_tokens = 4;
const char			tokens[]      = "+*-/"; // *-/!@#$%^&;()_={}[].<>?QWERTY";
//const int           arities[qtd_tokens]       = {2, 2, 2, 2}; //, 2, 2, 2, 4, 4, 4, 4}; 
//const int           arities[qtd_tokens]       = {3, 3, 3, 3}; //, 2, 2, 2, 4, 4, 4, 4}; 
const int           arities[qtd_tokens]       = {4, 4, 4, 4}; //, 2, 2, 2, 4, 4, 4, 4}; 
const Function		all_functions[qtd_tokens] = { 
//                                            eval_plus,  eval_mult, eval_minus, eval_div,     // for initializing the jumptable with binary functions
//                                            eval_plus_plus, eval_plus_minus, eval_plus_mult, eval_plus_div //,     // for initializing the jumptable with ternary functions
//                                            eval_minus_plus, eval_minus_minus, eval_minus_mult,// eval_minus_div,
//                                            eval_mult_plus, eval_mult_minus, eval_mult_mult, //eval_mult_div,
//                                            eval_div_plus, eval_div_minus, eval_div_mult, eval_div_div,
//                                            eval_plus_plus_plus, eval_plus_plus_minus, eval_plus_plus_mult, eval_plus_plus_div,      // for initializing the jumptable with quaternary functions
                                            eval_plus_plus_plus, eval_plus_plus_plus, eval_plus_plus_plus, eval_plus_plus_plus,      // for initializing the jumptable with quaternary functions
//                                            eval_plus_minus_mult, eval_plus_mult_minus//, eval_plus_mult_div, eval_plus_div_minus      // for initializing the jumptable with quaternary functions
                                        };

//		*/


