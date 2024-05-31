This code implements algorithms to solve the NSPSD Procrustes problem from the paper
(add preprint here) 

Given rectangular matrices X and B, the problem is to find the (non-)symmetric positive semidefinite matrix A that minimizes the Frobenius norm of AX-B. 
A non-symmetric PSD matric is a matrix A such that A+A^T is PSD. 

You can run the files in the folder test to compare the four approaches considered in the paper: ANFGM, FGM, MINGD and CARD.
Procrustes_ANFGM is the method defined in [BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric semidefinite Procrustes problem, 2021. 
Procrustes_FGM is the Fast Gradient Method By Nesterov.
Procrustes_MINGD implements the main algorithm of the paper, with minimum norm of the solution.
Prcorustes_CARD implements the main algorithm of the paper, with approximate minimum norm of the solution.

In the repository 'tests' you can find the codes for generating the numerical experiments in the paper.

Before running any code, please run Install.m to have all paths added. 



