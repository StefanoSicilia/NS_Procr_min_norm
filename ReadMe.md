Given rectangular matrices $X$ and $B$, the NSPDSP problem requires to find the matrix $A$ that minimizes the Frobenius norm of $AX-B$, under the constraint that the symmetric part of A is positive semi-definite, that is, $A+A^\top\succeq 0$. 

This code implements algorithms to solve the non-symmetric positive semidefinite Procrustes (NSPDSP) problem from the paper

[SG24] Sicilia, S. and Gillis, N., Minimum-norm solutions of the non-symmetric semidefinite Procrustes problem, June 2024. 


The four approaches considered in the paper are ANFGM, FGM, MINGD and CARD, which are associated to the following four functions: 
- Procrustes_ANFGM is the method proposed in [BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric semidefinite Procrustes problem, 2021. 
- Procrustes_FGM uses the Fast Gradient Method by Nesterov in [N83] Y. Nesterov, A method of solving a convex programming problem with convergence rate $\mathcal{O}(1/\kappa^2)$, 1983. 
- Procrustes_MINGD implements the main algorithm of our paper [SG24] (Algorithm 2), which guarantees to compute the minimum-norm solution. 
- Procrustes_CARD implements a variant of Algorithm 2 of the paper [SG24], which approximates the minimum-norm solution at a lower computational cost.  

In the repository 'tests' you can find the codes for generating the numerical experiments in the paper.

Before running any code, please run Install.m to have all paths added. 
