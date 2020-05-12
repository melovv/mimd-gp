# A MIMD Interpreter for Genetic Programming

This is the source code used in the paper:

```
https://link.springer.com/chapter/10.1007/978-3-030-43722-0_41
```
 
This is far from production code. It is based on:

```
tinygp, a minimalist genetic programming system in C++ ( (c) M.Keijzer 2004. License: GPL-2 or later)

```


## Abstract

Most Genetic Programming implementations use an interpreter to execute an individual, in order to obtain its outcome. Usually, such interpreter is the main bottleneck of the algorithm, since a single individual may contain thousands of instructions that must be executed on a dataset made of a large number of samples. Although one can use SIMD (Single Instruction Multiple Data) intrinsics to execute a single instruction on a few samples at the same time, multiple passes on the dataset are necessary to calculate the result. To speed up the process, we propose using MIMD (Multiple Instruction Multiple Data) instruction sets. This way, in a single pass one can execute several instructions on the dataset. We employ AVX2 intrinsics to improve the performance even further, reaching a median peak of 7.5 billion genetic programming operations per second in a single CPU core.

## Keywords

Genetic Programming Interpreter Vectorization Multiple Instruction Multiple Data 


