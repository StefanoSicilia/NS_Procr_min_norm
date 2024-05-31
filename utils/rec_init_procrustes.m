% Initialization for the PSD Procrustes problem 
% 
% min_{A or (A'+A) PSD} ||AX-B||_F, 
% 
% where X is diagonal with positive diagonal elements. 
% 
% If the problem is well conditioned (cond(X) <= paramcond), then 
% first-order method will converge relatively fast and we use a simple
% diagonal initialization (see init_procrustes.m) and perform a few
% iterations of FGM. 
% 
% If the problem is ill-conditioned (cond(X) > paramcond), then we split X 
% into two diagonal matrices X1 and X2 where X1 and X2 are better 
% conditioned: 
% [ A1 | 0 ] [ X1 | 0 ]       [ AX1 | 0  ]
%   ------    --------    =     --------         approx   B 
% [ 0  | A2] [ 0  | X2 ]      [  0  | AX2] 
% We also impose the corresponding off-diagonal blocks of A to be zero. 
% As long as a submatrix (X1 or X2) is ill-conditioned, we use the same
% procedure recursively in order to eventually partition X into 
% well-conditioned blocks. 
%
% See section 4.2.3 in the paper 
% A semi-analytical approach for the positive semidefinite Procrustes 
% problem, Nicolas Gillis and Punit Sharma, 2016.  

function A0 = rec_init_procrustes(B,X,paramcond,nspsd) 

if nargin <= 2
    paramcond = 100; % Treshold to decide whether a matrix is well-conditioned
end
if nargin <= 3
    nspsd = 1; % Treshold to decide whether a matrix is well-conditioned
end
n = size(X,1); 
% We do not assume the diagonal entries of X are sorted
x = diag(X); 
if n ~= size(X,2) || norm(diag(x) - X) > 1e-6
    error('The input matrix X should be diagonal.');
end
[xs,perms] = sort(x); 
X = X(perms,perms); 
B = B(perms,perms); 
[~,permsinv] = sort(perms); 
% If X is ill-conditioned, partition X: 
if xs(n)/xs(1) > paramcond
    [k,vk] = partition_cond(xs); 
    K1 = 1:k; 
    K2 = k+1:n; 
    A01 = rec_init_procrustes(B(K1,K1),X(K1,K1),paramcond,nspsd); 
    A02 = rec_init_procrustes(B(K2,K2),X(K2,K2),paramcond,nspsd); 
    A0 = [A01 zeros(k,n-k); zeros(n-k,k) A02]; 
% otherwise, run FGM for 100 iterations 
else
    optionsinit.maxiter = 100; 
    optionsinit.delta = 0.01; 
    optionsinit.nspsd = nspsd; 
    A0 = Procrustes_FGM(B,X,optionsinit);  
end 
X = X(permsinv,permsinv); 
A0 = A0(permsinv,permsinv);
B = B(permsinv,permsinv); 