# Matlab Codes of SWI, DIOM, DQGMRES

## How to cite
If you use swi in your work, please cite using the format given in [CITATION.bib](https://github.com/huangna112/SWI/blob/main/CITATION.bib)

## Content
We provide here some matlab codes of certain methods for solving unsymmetric positive definite linear systems Ax = b.

## Algorithms
- **SWI** Sliding window implementation of SCG with pre-allocated memory (When mk=n, it equals to SCG).
- **SWIWP** Sliding window implementation without pre-allocated memory  (When mk=n, it equals to SCG).
- **DIOM** Direct Incomplete Orthogonalization Method (Algorithms 6.6 and 6.8 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)")
- **DQGMRES** Direct Quasi-GMRES (DQGMRES) (Algorithms 6.6 and 6.13 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)")

## Test
### Convection-diffusion
The linear system arises from the convection diffusion equation. 
- **test_condiff** Runing this file gets the comparing results.


### SuiteSparse Matrix
We select matrix from the SuiteSparse Matrix Collection [https://sparse.tamu.edu/] and b = Ax\* with x\* = (1,1,...,1)'.
- **test_suitesparse.m**  Runing this file gets the comparing results. 
- **testproblem.txt** Storing all the file name of testing matrices.

