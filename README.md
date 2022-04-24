# Matlab Codes of SWI, DIOM, DQGMRES
We provide here some matlab codes of certain methods for solving unsymmetric positive definite linear systems Ax = b.
## Algorithms
- **SWI** Sliding window implementation of SCG with pre-allocated memory (When mk=n, it equals to SCG).
- **SWIWP** Sliding window implementation without pre-allocated memory  (When mk=n, it equals to SCG).
- **DIOM** Direct Incomplete Orthogonalization Method (Algorithms 6.6 and 6.8 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)")
- **DQGMRES** Direct Quasi-GMRES (DQGMRES) (Algorithms 6.6 and 6.13 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)")
## convection-diffusion
The linear system arises from the convection diffusion equation. Running "test_condiff" file can get the results.
## 
We select matrix from the SuiteSparse Matrix Collection [https://sparse.tamu.edu/] and b = Ax\* with x\* = (1,1,...,1)\T
- **test_suitesparse.m**  
- **testproblem.txt** All the file name of testing matrices.
