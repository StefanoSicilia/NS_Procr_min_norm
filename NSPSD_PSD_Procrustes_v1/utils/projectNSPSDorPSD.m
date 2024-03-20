% Project the matrix Q onto the 
% NSPSD cone if nspsd = 1, 
% PSD cone otherwise. 

function Qp = projectNSPSDorPSD(Q,nspsd,iscomplex) 

if nargin <= 2 || nspsd == 0 || iscomplex == 0
    if nspsd == 1
        Qp = projectNSPSD(Q);
    else
        Qp = projectPSD(Q);
    end
else
    n = size(Q,1)/2; 
    Q1 = Q(1:n,1:n); 
    Q1H = (Q1+Q1')/2; 
    Q1S = (Q1-Q1')/2; 
    Q2 = Q(1:n,n+1:2*n); 
    Q2H = (Q2+Q2')/2;
    Q2S = (Q2-Q2')/2; 
    Q3 = Q(n+1:2*n,1:n); 
    Q3H = (Q3+Q3')/2; 
    Q3S = (Q3-Q3')/2; 
    Q4 = Q(n+1:2*n,n+1:2*n); 
    Q4H = (Q4+Q4')/2;
    Q4S = (Q4-Q4')/2;
  % X1X2 = projectNSPSDorPSD( 0.5*(Q1H+Q4H) + 0.5*sqrt(-1)*(Q2S-Q3S) ,1);
    X1X2 = projectPSD( 0.5*(Q1H+Q4H) + 0.5*sqrt(-1)*(Q2S-Q3S));
    X1H = real(X1X2); 
    X2S = imag(X1X2); 
    X1S = (Q1S+Q4S)/2; 
    X2H = (Q2H-Q3H)/2; 
    Qp = [ X1H+X1S X2S+X2H; -X2S-X2H X1H+X1S ]; 
end
end

% Project the matrix Q onto the PSD cone 
%
% This requires the eigendecomposition of (Q+Q')/2 and then setting the 
% negative eigenvalues to zero. 
function Qp = projectPSD(Q) 
Q = (Q+Q')/2; 
if max(max(isnan(Q))) == 1 || max(max(isinf(Q))) == 1
    error('Input matrix has infinite or NaN entries');
end
[V,e] = eig(Q); 
Qp = V*diag(max(diag(e),0))*V'; 
end

% Project the matrix Q onto the set of NSPSD matrices  
% 
% This requires the eigendecomposition of (Q+Q')/2 and then setting the 
% negative eigenvalues to zero. 
function Qp = projectNSPSD(Q) 
Qh = (Q+Q')/2; 
Qs = (Q-Q')/2;
if max(max(isnan(Q))) == 1 || max(max(isinf(Q))) == 1
    error('Input matrix has infinite or NaN entries');
end
[V,e] = eig(Qh); 
Qp = V*diag(max(diag(e),0))*V' + Qs; 
end