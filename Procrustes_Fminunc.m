function [Ahat,es,ts,status] = Procrustes_Fminunc(X,B,options)
%% Procrustes_Fminunc:
% Follows procrustes_anly but the output matrix Ahat is built differently.
% It uses fminunc function to minimize the norm of Ahat in a certain form.
% See procrustes_anyl for further info.

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
if ~isfield(options,'alpha')
    options.alpha = 0.3;
end

% Set up the complex problem
if options.nspsd == 1 && options.complex == 1
    % Construct the equivalent 2n x 2n real problem;
    % see Theorem 5 in the paper
    B = [real(B)  imag(B); -imag(B)  real(B)];
    X = [real(X)  imag(X); -imag(X)  real(X)];
    % n = size(X,1);
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
    if ~any(temp(:)>0)
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
        Ahat=U1*A11*U1'+2*skew(U2*(U2'*B*V1/D1)*U1');
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
        es = sqrt(es.^2 + (norm(B*V2,'fro'))^2);
        H11 = (A11+A11')/2;
        S11 = (A11-A11')/2;
        % depending on H11, we can conclude whether the infimum is
        % attained and compute a optimal solution accordingly.
        % ***CASE 1*** infimum attained if and only if
        % H11 positive definite
        % OR
        % Ker(H11) in ker((1/2)*U2'*B*V1/D1)
        Z=0.5*U2'*B*V1/D1;
        if ~any(eig(H11)<10^-9) || norm(Z*null(H11),'fro')<1e-6
            Ahat=U1*A11*U1'+4*skew(U2*Z*U1');
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is attained by Ahat';
        else % ***CASE 2*** infimum not attained
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
            W1=U11(:,s+1:r);
            Lambda=D11(s+1:r,s+1:r);
            % First alternative
%             f=@(Y) 4*norm(W1*W1'*Y'-Z','fro')^2+norm(Y*W1/Lambda*W1'*Y','fro')^2;
%             Y=fminunc(f,Z);
%             N=Y*W1/Lambda*W1'*Y';
%             M=W1*W1'*Y'-Z';
%             Ahat=U1*Aeps*U1'+2*U2*Z*U1'+2*U1*M*U2'+U2*N*U2';
             % Second alternative
%             Z=U2'*B*V1/D1;
%             f=@(Y) 1/16*norm(Y/Lambda*Y','fro')^2+norm(Y*W1'-Z,'fro')^2;
%             Y=fminunc(f,Z*W1);
%             N=1/4*Y/Lambda*Y';
%             M=W1*Y'-Z';
%             Ahat=U1*Aeps*U1'+U2*Z*U1'+U1*M*U2'+U2*N*U2';
            % Third alternative
            Z=U2'*B*V1/D1;
            Lsq=sqrtm(inv(Lambda));
            f=@(Y) 1/16*norm(Lsq*Y'*Y*Lsq,'fro')^2+norm(Y*W1'-Z,'fro')^2;
            Y=fminunc(f,Z*W1);
            N=1/4*Y/Lambda*Y';
            M=W1*Y'-Z';
            Ahat=U1*Aeps*U1'+U2*Z*U1'+U1*M*U2'+U2*N*U2';
            %
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is not attained for any NSPSD A because H11 not > 0 and Ker(H11) not in ker((1/2)*U2''*B*V1/D1)';
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
    if ~any(temp(:)>0)
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
        K = ((U2'*B*V1/D1)/A11)*(U2'*B*V1/D1)';
        Ahat = U1*A11*U1' + U2*(U2'*B*V1/D1)*U1' + U1*(U2'*B*V1/D1)'*U2'+U2*K*U2';
        ubound=norm(Ahat*X-B,'fro');
        status = 'inf(||AX-B||_F) is not attained for any PSD A because U1''*(B*X''+X*B'')*U1 < 0 (Theorem 3.2)';
        ts = cputime - tmimec;
    else % Theorem 3.1: the subproblem in A11 has to be solved
        Bprime = U1'*B*V1;
        % Use recursive initialization
        A0rec = rec_init_procrustes(Bprime,D1,100,0);
        % Use FGM to solve the subproblem
        optionsFGM = options; 
        optionsFGM.complex = 0; 
        optionsFGM.A0 = A0rec;
        optionsFGM.timemax = options.timemax-(cputime-tmimec);
        [A11,es,ts] = Procrustes_FGM(Bprime,D1,optionsFGM);
        es = sqrt(es.^2 + (norm(B*V2,'fro'))^2);
        % Theorem 3.1: depending on A11, we can conclude whether the infimum is
        % attained and compute a optimal solution accordingly.
        % ***CASE 1*** infimum attained if and only if
        % A11 positive definite
        % OR
        % Ker(A11) in ker(U2'*B*V1/D1)
        Z=U2'*B*V1/D1;
        if ~any(eig(A11)<10^-9) || norm(Z*null(A11),'fro')<1e-6
            K=Z*pinv(A11)*Z';
            Ahat=U1*A11*U1' + U2*Z*U1' + U1*Z'*U2'+U2*K*U2';
            ubound=norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is attained by Ahat';
        else % ***CASE 2*** infimum not attained
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
            K = ((U2'*B*V1/D1)/Temp)*(U2'*B*V1/D1)';
            Ahat = U1*Temp*U1' + U2*(U2'*B*V1/D1)*U1' + U1*(U2'*B*V1/D1)'*U2'+U2*K*U2';
            ubound = norm(Ahat*X-B,'fro');
            status = 'inf(||AX-B||_F) is not attained for any PSD A because A11 not > 0 and Ker(A11) not in ker(U2''*B*V1/D1)';
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