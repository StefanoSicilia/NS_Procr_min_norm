This code implements algorithms to solve the NSPSD Procrustes problem from the paper
(add preprint here) 

Given rectangular matrices X and B, the problem is to find the (non-)symmetric positive semidefinite matrix A that minimizes the Frobenius norm of AX-B. 
A non-symmetric PSD matric is a matrix A such that A+A^T is PSD. 

The four approaches considered in the paper are ANFGM, FGM, MINGD and CARD, which are associated to the following four functions:
\begin{itemize}
\item Procrustes_ANFGM is the method defined in [BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric semidefinite Procrustes problem, 2021. 
\item Procrustes_FGM is the Fast Gradient Method by Nesterov in [N83] Y. Nesterov, A method of solving a convex programming problem with convergence rate $\mathcal{O}(1/\kappa^2)$, 1983. 
\item Procrustes_MINGD implements the main algorithm of the paper, with minimum norm of the solution.
\item Procrustes_CARD implements the main algorithm of the paper, with approximate minimum norm of the solution.
\end{itemize}

In the repository 'tests' you can find the codes for generating the numerical experiments in the paper.

Before running any code, please run Install.m to have all paths added. 



