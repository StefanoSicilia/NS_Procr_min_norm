% Solves the PSD and NSPSD Procrustes problem 
% 
% min_{ A } ||AX - B||_F , 
%     such that A+A' is PSD when nspsd=1 (NSPSD Procrustes), 
%     such that A    is PSD otherwise (PSD Procrustes),  
% using CVX 
% 
% Input: Matrices B and X, 
%        options.maxiter = maximum number of iterations of the IPM solver 
%        options.nspsd   = 1 if the NSPSD Procrustes is solved
%        options.complex = 1 if the complex NSPSD Procrustes is solved

function [A,e,t,cvx_status] = Procrustes_cvx(B,X,options)

if nargin <= 2
    options = [];
end
if ~isfield(options,'nspsd')
    options.nspsd = 1; 
end
if ~isfield(options,'complex')
    options.complex = 0; 
end

% Set up the complex problem 
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem; 
    % see Theorem 5 in the paper 
    B = [real(B)  imag(B); -imag(B)  real(B)]; 
    X = [real(X)  imag(X); -imag(X)  real(X)]; 
    n = size(X,1); 
end

t = cputime; 
n = size(X,1); 
cvx_begin sdp quiet 
	if isfield(options,'maxiter')
        cvx_solver_settings('maxit',options.maxiter);
    end
    variable A(n,n);
    minimize( norm( A*X - B , 'fro')  );
    subject to
    if options.nspsd == 1
        (A+A') == hermitian_semidefinite( n ); 
    else
        A == hermitian_semidefinite( n ); 
    end
cvx_end
A = projectNSPSDorPSD(A,options.nspsd); 
t = cputime-t; 
e = norm(A*X-B,'fro'); 

% Put back the complex problem together 
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem; 
    % see Theorem 5 in the paper 
    n = n/2; 
    A = A(1:n,1:n) + sqrt(-1) * A(1:n,n+1:2*n); 
end