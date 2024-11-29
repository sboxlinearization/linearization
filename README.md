# Algorithms for linearization of S-boxes

This repository contains supporting code for the paper.

## Greedy Extension algorithm

This algorithm is heuristic and was use to produced most of our best approximations.
The code is written in [greedy_extension.py](./greedy_extension.py), it also uses basic classes from other included modules.

Results for monomial S-boxes are given in [results_monomials.txt](./results_monomials.txt). In the form of json `{exp: [approx_for_n4, ..., approx_for_n13]}`, where each approximation is given by the intersection S ∩ A, i.e. the x-coordinates of points of S. The monomial S-boxes were generated using sage's `sage.crypto.sboxes.monomial_function` method.


## Exhaustive Reduction algorithm

This algorithm is used to exhaust the search space and yield guaranteed vectorial linearity (or an upper bound). Its basic idea is to recursively enumerate candidates for the coordinates of the affine approximation, one at a time. At the same time, it maintains the set of satsifying inputs for all the previous guesses. If the set of inputs is lower than the lower bound LB (parameter), then the current subtree is cut and not searched anymore. It also includes an optimization for quadratic maps, by using the derivatives to reduces the search space (typically by a factor 2^n).

The code is written in [ReductionExhaustive.cpp](./ReductionExhaustive.cpp). The [makefile](./makefile) contains options for some hardcoded S-boxes, as well as `RE-custom` target which can process an S-box from an input (up to 8 bits), and `RE-custom-large` (same for up to 16 bits). For example, it can be ran as follows:

```bash
$ make RE-custom
$ cat sboxes/n5_Ascon
4 11 31 20 26 21 9 2 27 5 8 18 29 3 6 28 30 19 7 14 0 13 17 24 16 12 1 25 22 10 15 23
# ./RE-custom <n input bits> <LB> <batch index> <batch count>
$ ./RE-custom 5 12 0 1 <sboxes/n5_Ascon  # Ascon S-box LB=12 - finds nothing
...
$ ./RE-custom 5 11 0 1 <sboxes/n5_Ascon  # Ascon S-box LB=12 - finds 2 solutions
...
S-box is quadratic? 1 1
Computed 32 derivatives (16 unique)
...
SOLUTION depth 5: 11 : 1 2 5 13 16 17 18 20 21 28 30 
SOLUTION depth 5: 11 : 1 5 6 9 16 17 20 21 22 24 26
...

$ ./RE-custom 8 42 0 1 <sboxes/n8_AES
...
# ~3 minutes with LB 42 on a laptop (8 cores)
# much longer with LB 19 (but feasible)
```

The 2 solutions together with translation by constant produce all 32 solutions.

NOTE: We only take first 2 output bits of derivatives which for Ascon produces only 16 (out of 32) unique functions. Perhaps randomizing the S-box with affine maps before running this code could help to use all 32 available translations.
