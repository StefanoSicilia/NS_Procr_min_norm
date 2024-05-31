This code implements algorithms to solve the NSPSD Procrustes problem from the paper
(add preprint here) 
[BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric semidefinite Procrustes problem, 2021. 

Given rectangular matrices X and B, the problem is to find the (non-)symmetric positive semidefinite matrix A that minimizes the Frobenius norm of AX-B. 
A non-symmetric PSD matric is a matrix A such that A+A^T is PSD. 

You can run the files in the folder test to compare the four approaches considered in the paper: ANFGM, FGM, MINGD and CARD.

Before running any code, please run Install.m to have all paths added. 



