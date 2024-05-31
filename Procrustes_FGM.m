% Given n-by-m matrices X and B, this code aims to find, if possible, 
% 
% If options.nspsd = 1 
%  the n-by-n matrix A with A+A' positive semidefinite (PSD) 
% Otherwise 
%   the n-by-n PSD matrix A 
% 
% such that ||AX - B||_F is minimized. 
% 
% ****** Description ******
% This code tries to solve the following optimization problem
% 
% min_{A} || A X - B ||_F such that (A+A') or A is PSD. 
% 
% This problem is referred to as the (NS)PSD Procrustes problem. 
%
% It uses the *fast gradient method* from smooth convex optimization; 
% see p.90 of Yu. Nesterov, Introductory lectures on convex optimization: 
% A basic course, vol. 87, Springer Science & Business Media, 2004.
% 
% See the papers  
% N. Gillis and P. Sharma, A semi-analytical approach for the positive 
% semidefinite procrustes problem, Linear Algebra and its Applications 540, 
% pp. 112-137, 2018, and 
% On non-Hermitian/non-symmetric semidefinite mapping and Procrustes 
% problem, Mohit Kumar Baghel, Nicolas Gillis, Punit Sharma, 2021. 
% 
% ****** Input ******
% X, B     : n-by-m matrices X and B
% 
% -- options ---
% .nspsd    : options.nspsd  = 1 => NSPSD Procrustes solved  
%             otherwise the PSD Procrustes solved (default) 
% .complex  : options.complex  = 1 => complex NSPSD Procrustes solved 
%             (see Theorem 5 in the paper), otherwise the PSD Procrustes is 
%             solved (default) 
% .A0       : n-by-n initial iterate
%            -default: use the initialization from init_procrustes.m
% .maxiter  : the maximum number of iterations performed by the algorithm 
%            -default = 1000. 
% .delta    : stops the algorithm if ||A^k - A^0||<=delta*||A^1 - A^0||
%            where A^k is the kth iterate
%            -default = 0 (inactive)
% .timemax  : time alloted to the algorithm 
%            -default: 60 seconds. 
% .alpha0   : Parameter of the FGM in (0,1) -it can be tuned. 
% 
%
% ****** Output ******
% A        : an n-by-n PSD matrix that minimizes ||AX-B||_F 
% e        : evolution of the approximation error (Frobenius norm)
% t        : cputime to compute e --plot(t,e) displays the error over time 

function [A,e,t] = Procrustes_FGM(B,X,options) 

timeinit = cputime; 
n = size(X,1); 
if nargin <= 2
    options = [];
end
if ~isfield(options,'nspsd')
    options.nspsd = 0; % default 
end
if ~isfield(options,'complex')
    options.complex = 0; % default 
end
if ~isfield(options,'maxiter')
    options.maxiter = 1e4; 
end
if ~isfield(options,'delta')
    options.delta = 1e-6; 
end
if ~isfield(options,'timemax')
    options.timemax = 60; 
end
if ~isfield(options,'alpha0')
    options.alpha0 = 0.1; % Parameter of the FGM in (0,1) -it can be tuned. 
end
% Set up the complex problem 
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem; 
    % see Theorem 5 in the paper 
    B = [real(B)  imag(B); -imag(B)  real(B)]; 
    X = [real(X)  imag(X); -imag(X)  real(X)]; 
    n = size(X,1); 
end
 
if ~isfield(options,'A0')
    A0 = init_procrustes(B,X,2,options.nspsd); 
else
    A0 = projectNSPSDorPSD(options.A0, options.nspsd); 
end
A = A0; 


% Precomputations
% If X is diagonal, it can simplify computations
if size(X,1)==size(X,2) && norm( X - diag(diag(X)) ) < 1e-12
    diagx = 1; 
    x = diag(X); 
    xsqa = x.^2; 
    Lx = max(xsqa);  % Lipschitz constant
    mux = min(xsqa); % strongly convex constant 
    qx = mux/Lx; 
    BXt = B.*repmat(x',n,1); 
else
    diagx = 0; 
    XXt = X*X'; 
    [u,x] = eig(XXt); 
    x = diag(x); 
    Lx = max(x);  % Lipschitz constant
    mux = min(x); 
    qx = mux/Lx; 
    BXt = B*X'; 
end
% Parameters and initalization
alpha(1) = options.alpha0;
Y = A; 
i = 1; 
eps0 = 0; eps = 1;  
t(1) = cputime - timeinit;
tt = cputime; % Do not take into account evaluation of the objective in 
              % the computational cost
if nargout >= 2
    if diagx == 1 
        e(1) = norm( A.*repmat(x',n,1) - B , 'fro' ); 
    else
        e(1) = norm( A*X - B , 'fro');  
    end
    tt = cputime-tt; 
end
while i <= options.maxiter && eps >= options.delta*eps0 ... 
        && t(i) <= options.timemax 
    % Previous iterate
    Ap = A; 
    % FGM Coefficients  
    if alpha(i) == 0
         alpha(i+1) = 0; 
        beta(i) = 0;
    else
        alpha(i+1) = ( sqrt( (alpha(i)^2 - qx)^2 + 4*alpha(i)^2 ) + (qx-alpha(i)^2) ) / (2);  
        beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
    end
    % Projected gradient step from Y
    if diagx == 1
        A = projectNSPSDorPSD( Y - (Y.*repmat(xsqa',n,1) - BXt) / Lx, options.nspsd, options.complex );
    else
        A = projectNSPSDorPSD( Y - (Y*XXt - BXt) / Lx, options.nspsd, options.complex );
    end
    % Linear combination of iterates
    Y = A + beta(i)*(A-Ap);  
    % Error: only computed if required because it takes O(mn^2) operations
    t(i+1) = cputime - timeinit - tt; 
    if nargout >= 2
        tt = cputime; 
        if diagx == 1
            e(i+1) = norm(A.*repmat(x',n,1) - B, 'fro');  
        else
            e(i+1) = norm(A*X - B, 'fro');  
        end
        tt = cputime - tt; 
    end
    % Stopping condition depending on the distance between iterates 
    if i == 1
        eps0 =  norm(A-Ap,'fro'); 
    end
    eps = norm(A-Ap,'fro'); 
    i = i + 1;
end 

% Put back the complex problem together 
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem; 
    % see Theorem 5 in the paper 
    n = n/2; 
    A = A(1:n,1:n) + sqrt(-1) * A(1:n,n+1:2*n); 
end