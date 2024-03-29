# MATLAB code of some Krylov methods for solving unsymmetric positive definite linear systems

## How to cite
If you use swi in your work, please cite using the format given in [CITATION.bib](https://github.com/huangna112/SWI/blob/main/CITATION.bib)

## Content
We provide here some matlab code of certain methods for solving unsymmetric positive definite linear systems Ax = b.

## Algorithms
- **REGMRES** Restarted GMRES (REGMRES) (Algorithm 6.11 in Yousef Saad's [Iterative Methods for Sparse Linear System (2nd Edition)](https://epubs.siam.org/doi/book/10.1137/1.9780898718003)) (When restart=n, it equals to GMRES).
- **REFOM** Restarted Full Orthogonalization Method (Algorithms 6.5 in Yousef Saad's [Iterative Methods for Sparse Linear System (2nd Edition)](https://epubs.siam.org/doi/book/10.1137/1.9780898718003)) (When restart=n, it equals to FOM). 
- **RESCG** Restarted Semi-Conjugate Gradient Method (When restart=n, it equals to SCG).
- **DQGMRES** Direct Quasi-GMRES (DQGMRES) (Algorithms 6.6 and 6.13 in Yousef Saad's [Iterative Methods for Sparse Linear System (2nd Edition)](https://epubs.siam.org/doi/book/10.1137/1.9780898718003)) (When m=n, it equals to GMRES).
- **DIOM** Direct Incomplete Orthogonalization Method (Algorithms 6.6 and 6.8 in Yousef Saad's [Iterative Methods for Sparse Linear System (2nd Edition)](https://epubs.siam.org/doi/book/10.1137/1.9780898718003)) (When m=n, it equals to FOM).
- **SWI** Sliding Window Implementation of SCG with pre-allocated memory (When m=n, it equals to SCG).
- **DYNSWIWP** Sliding Window Implementation of SCG with choosing m dynamically.


## Test
### Convection-diffusion
The linear system arises from the convection diffusion equation. 
- **test_condiff** Runing this file gets the comparing results, for more details, please see [test_condiff.pdf](https://github.com/huangna112/SWI/blob/main/test_condiff.pdf)


### SuiteSparse Matrix
We select matrices from the [SuiteSparse Matrix Collection](https://sparse.tamu.edu/) and b = Ax\* with x\* = (1,1,...,1)'.
- **test_suitesparse.m**  Runing this file gets the comparing results, for more details, please see [test_suitesparse.pdf](https://github.com/huangna112/SWI/blob/main/test_suitesparse.pdf). 
- **testproblem.txt** Storing all the file name of testing matrices. 
