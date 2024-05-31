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
% This problem is referred to as the (N)SPSD Procrustes problem.
%
% It combines:
% (1) a semi-analytical approach that reduces to the case where X is
% diagonal with positive diagonal elements and of dimension rank(X). This
% subproblem has a unique optimal solution.
% (2) a recursive initialization procedure that divides the problem, if it
% is ill-conditioned, into smaller well-conditioned subproblems; see
% rec_init_procrustes.m.
% (3) the fast gradient method, an optimal first-order method for smooth
% convex optimization; see Procrusted_FGM.m.
%
% See the papers
% [GP18] N. Gillis and P. Sharma, A semi-analytical approach for the positive 
% semidefinite procrustes problem, Linear Algebra and its Applications 540 
% (2018), pp. 112-137, and 
% [BGP21] M.K. Baghel, N. Gillis and P. Sharma, On the non-symmetric 
% semidefinite Procrustes problem, 2021. 
%
% ****** Input ******
% X, B     : n-by-m matrices X and B
% --- options --
% .nspsd    : options.nspsd  = 1 => NSPSD Procrustes solved; see [BGP21] 
%             otherwise the symmetric PSD Procrustes is solved; see [GP18]
% .complex  : options.complex  = 1 => complex NSPSD Procrustes solved 
%             (see Theorem 5 in [BGP21]), otherwise the PSD Procrustes is 
%             solved as in [GP18] (default)
% .maxiter  : the maximum number of iterations performed by the algorithm
%            -default = 1000.
% .delta    : stops the FGM algorithm if ||A^k - A^0||<=delta*||A^1 - A^0||
%            where A^k is the kth iterate
%            -default = 0 (inactive)
% .timemax  : time alloted to the algorithm
%            -default: 60 seconds.
%
% ****** Output ******
% A        : an n-by-n NSPSD matrix (if options.nspsd = 1), an n-by-n PSD
%            matrix otherwise that minimizes ||AX-B||_F
% e        : evolution of the approximation error (Frobenius norm)
% t        : cputime to compute e --plot(t,e) displays the error over time
% status   : indicates whether the PSD Procrustes problem attains its
%            infimum. If it does not, then the algorithms returns a
%            solution whose objective function value is close the infimum.

function [Ahat,es,ts,status] = Procrustes_ANFGM(X,B,options)

tmimec = cputime;

if nargin <= 2
    options = [];
end
if ~isfield(options,'nspsd')
    options.nspsd = 1;
end
if ~isfield(options,'complex')
    options.complex = 0;
end
if ~isfield(options,'maxiter')
    options.maxiter = 1000;
end
if ~isfield(options,'delta')
    options.delta = 0;
end
if ~isfield(options,'timemax')
    options.timemax = 60;
end

% Set up the complex problem
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem;
    % see Theorem 5 in the paper
    B = [real(B)  imag(B); -imag(B)  real(B)];
    X = [real(X)  imag(X); -imag(X)  real(X)];
    n = size(X,1);
end

% Compute the SVD of X, to reduce the original problem
r=rank(X);
[n,m]=size(X);
[U,D,V]=svd(X);
D1=D(1:r,1:r);
U1=U(:,1:r);
U2=U(:,r+1:n);
V1=V(:,1:r);
V2=V(:,r+1:m);

if options.nspsd == 1 % NSPSD PROCRUSTES SOLVED; see [BGP21]
    % This is the case U1'*(B*X'+X*B')*U1 < 0 and
    % the infimum is not attained, but one can
    % compute explicitely a solution whose objective function value is close
    % to the infimum (up to any precision).
    temp = eig(U1'*(B*X'+X*B')*U1);
    if any(temp(:)>0) == 0
        es = sqrt((norm(U1'*B*V1,'fro'))^2 + (norm(B*V2,'fro'))^2);
        alpha = 4*sqrt(n)*norm(D1,'fro')*norm(U1'*B*V1,'fro');
        tol = 1e-6; % Precision of the solution as the infimum is not attained
        if tol/alpha < 1e-9 % This avoids numerical problems
            tol = 1e-9*alpha;
        end
        A11 = zeros(r);
        for j=1:r
            A11(j,j) = tol/alpha;
        end
        K = (U2'*B*V1*inv(D1))*inv(A11)*(U2'*B*V1*inv(D1))';
        Ahat = U1*A11*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U2*K*U2';
        ubound=norm(Ahat*X-B,'fro');
        status = 'inf(||AX-B||_F) is not attained for any NSPSD A because U1''*(B*X''+X*B'')*U1 < 0 ';
        ts = cputime - tmimec;
    else %  the subproblem in A11 has to be solved
        Bprime = U1'*B*V1;
        % Use recursive initialization
        A0rec = rec_init_procrustes(Bprime,D1,100,1);
        % Use FGM to solve the subproblem
        optionsFGM = options; 
        optionsFGM.complex = 0; 
        optionsFGM.nspsd = 1;
        optionsFGM.A0 = A0rec;
        optionsFGM.timemax = options.timemax-(cputime-tmimec);
        [A11,es,ts] = Procrustes_FGM(Bprime,D1,optionsFGM);
        e = es(end);
        es = sqrt(es.^2 + (norm(B*V2,'fro'))^2);
        H11 = (A11+A11')/2;
        S11 = (A11-A11')/2;
%         norm(A11,'fro')
%         cond(H11)
%         cond(D1)
        % depending on H11, we can conclude whether the infimum is
        % attained and compute a optimal solution accordingly.
        % ***CASE 1*** infimum attained if and only if
        % H11 positive definite
        % OR
        % Ker(H11) in ker((1/2)*U2'*B*V1*inv(D1))
        if any(eig(H11)<10^-9) == 0 || norm((1/2)*U2'*B*V1*inv(D1) * null(H11) ,'fro') < 1e-6
            K = (1/4)*(U2'*B*V1*inv(D1))*pinv(H11)*(U2'*B*V1*inv(D1))';
            Ahat = U1*A11*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U2*K*U2';
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is attained by Ahat';
        else % ***CASE 2*** infimim not attained
            H11=(H11+H11')/2;
            [U11,D11] = eig(H11);
            d = diag(D11);
            s = length(d(d<10^-12)); % that means rank(H11) = r-s;
            if norm(A11*D1-U1'*B*V1,'fro') < 1e-9
                beta = 4*sqrt(s)*norm(D1,'fro');
            else
                beta = 4*sqrt(s)*norm(D1,'fro')*norm(A11*D1-U1'*B*V1,'fro');
            end
            tol = 1e-6;  % Precision of the solution if the infimum is not attained
            if tol/beta < 1e-9 % This avoids numerical problems
                tol = 1e-9*beta;
            end
            Heps = D11;
            for j = 1 : s
                Heps(j,j) = tol/beta;
            end
            Heps = U11*Heps*U11';
            Aeps = Heps + S11;
            K = (1/4)*(U2'*B*V1*inv(D1))*inv(Heps)*(U2'*B*V1*inv(D1))';
            Ahat = U1*Aeps*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U2*K*U2';
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is not attained for any NSPSD A because H11 not > 0 and Ker(H11) not in ker((1/2)*U2''*B*V1*inv(D1))';
        end
    end
    
    % Put back the complex problem together
    if options.nspsd == 1 && options.complex == 1
        % Construct the equivalent 2n x 2n real problem;
        % see Theorem 5 in the paper
        n = n/2;
        Ahat = Ahat(1:n,1:n) + sqrt(-1) * Ahat(1:n,n+1:2*n);
    end
    
else % PSD PROCRUSTES SOLVED; see [GP18]
    % This is the case of Theorem 3.2 in [GP18] with U1'*(B*X'+X*B')*U1 < 0.
    % In short, this is a case where the infimum is not attained, but one can
    % compute explicitely a solution whose objective function value is close
    % to the infimum (up to any precision).
    temp=eig(U1'*(B*X'+X*B')*U1);
    if any(temp(:)>0) == 0
        es = sqrt((norm(U1'*B*V1,'fro'))^2 + (norm(B*V2,'fro'))^2);
        alpha = 4*sqrt(n)*norm(D1,'fro')*norm(U1'*B*V1,'fro');
        tol = 1e-6; % Precision of the solution as the infimum is not attained
        if tol/alpha < 1e-9 % This avoids numerical problems
            tol = 1e-9*alpha;
        end
        A11=zeros(r);
        for j=1:r
            A11(j,j)=tol/alpha;
        end
        K = (U2'*B*V1*inv(D1))*inv(A11)*(U2'*B*V1*inv(D1))';
        Ahat = U1*A11*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U1*(U2'*B*V1*inv(D1))'*U2'+U2*K*U2';
        ubound=norm(Ahat*X-B,'fro');
        status = 'inf(||AX-B||_F) is not attained for any PSD A because U1''*(B*X''+X*B'')*U1 < 0 (Theorem 3.2)';
        ts = cputime - tmimec;
    else % Theorem 3.1: the subproblem in A11 has to be solved
        Bprime = U1'*B*V1;
        % Use recursive initialization
        A0rec = rec_init_procrustes(D1,Bprime,100,0);
        % Use FGM to solve the subproblem
        optionsFGM = options; 
        optionsFGM.complex = 0; 
        optionsFGM.A0 = A0rec;
        optionsFGM.timemax = options.timemax-(cputime-tmimec);
        [A11,es,ts] = Procrustes_FGM(D1,Bprime,optionsFGM);
        e = es(end);
        es = sqrt(es.^2 + (norm(B*V2,'fro'))^2);
        % Theorem 3.1: depending on A11, we can conclude whether the infimum is
        % attained and compute a optimal solution accordingly.
        % ***CASE 1*** infimum attained if and only if
        % A11 positive definite
        % OR
        % Ker(A11) in ker(U2'*B*V1*inv(D1))
        if any(eig(A11)<10^-9) == 0 || norm(U2'*B*V1*inv(D1) * null(A11) ,'fro') < 1e-6
            K=(U2'*B*V1*inv(D1))*pinv(A11)*(U2'*B*V1*inv(D1))';
            Ahat=U1*A11*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U1*(U2'*B*V1*inv(D1))'*U2'+U2*K*U2';
            ubound=norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is attained by Ahat';
        else % ***CASE 2*** infimim not attained
            A11 = (A11+A11')/2;
            [U11,D11] = eig(A11);
            d = diag(D11);
            s = length(d(d<10^-12));
            beta = 4*sqrt(r-s)*norm(D1,'fro')*norm(A11*D1-U1'*B*V1,'fro');
            tol = 1e-6;  % Precision of the solution if the infimum is not attained
            if tol/beta < 1e-9 % This avoids numerical problems
                tol = 1e-9*beta;
            end
            Temp = D11;
            for j = 1 : s
                Temp(j,j) = tol/beta;
            end
            Temp = U11*Temp*U11';
            K = (U2'*B*V1*inv(D1))*inv(Temp)*(U2'*B*V1*inv(D1))';
            Ahat = U1*Temp*U1' + U2*(U2'*B*V1*inv(D1))*U1' + U1*(U2'*B*V1*inv(D1))'*U2'+U2*K*U2';
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is not attained for any PSD A because A11 not > 0 and Ker(A11) not in ker(U2''*B*V1*inv(D1))';
        end
    end
    
end

% For highly ill-conditioned problems, it may happen that the analytical
% approach does not lead to an acceptable solution, because of numerical
% instability. This can be detected by comparing es (semi-analytical
% error) with the actual final error of Ahat (ubound).
if abs( ubound - es(end) ) / ( norm(B,'fro')+1e-6 ) > 0.01 % Relative error > 1%
    warning('The problem is highly ill-conditioned and the analytical error does not coincide with the computed solution.');
end
ts = ts + cputime - tmimec - ts(end);
end