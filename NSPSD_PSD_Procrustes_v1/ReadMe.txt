This code implements algorithms to solve the PSD Procrustes problem from the paper 
[GP18] N. Gillis and P. Sharma, A semi-analytical approach for the positive  semidefinite procrustes problem, Linear Algebra and its Applications 540  (2018), pp. 112-137,
and for the non-symmetric PSD Procrustes problem from the paper 
[BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric % semidefinite Procrustes problem, 2021. 

Given rectangular matrices X and B, the problem is to find the (non-)symmetric positive semidefinite matrix A that minimizes the Frobenius norm of AX-B. 
A non-symmetric PSD matric is a matrix A such that A+A^T is PSD. 

You can run the files in the folder test to compare our proposed semi-analytical algorithm (procrustes_anly.m) which uses a fast gradient method (Procrustes_FGM.m) and a clever initialization scheme, as done in [BGP21]. 

We have also solved this problem with CVX (Procrustes_cvx.m); see http://cvxr.com/cvx/ to download this toolbox. 

Before running any code, pelase run Install.m to have all paths added. 

Please report any bug to nicolas.gillis@umons.ac.be 


