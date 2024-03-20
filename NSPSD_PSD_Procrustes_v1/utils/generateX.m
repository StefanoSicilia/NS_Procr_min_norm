% Generate the n-by-m matrix X in different ways: 
% typeR = 1: well-conditioned  - X = randn(n,m) 
% typeR = 2: ill-conditioned   - cond(X) = 10^6
% typeR = 3: rank-deficient    - rank(X) = min(m,n)/2

function X = generateX(n,m,typeR); 

% Well-conditioned
if typeR == 1
    X = randn(n,m); 
% Ill-conditioned
elseif typeR == 2
    X = randn(n,m); 
    [u,s,v] = svd(X,'econ'); 
    k = min(m,n); 
    alpha = 10^(6/(k-1)); 
    s = alpha.^[0:(k-1)]; 
    X = u*diag(s)*v';
% Rank-deficient
elseif typeR == 3
    X = randn(n,m); 
    [u,s,v] = svd(X,'econ'); 
    r = ceil( min(m,n)/2 ); 
    s = diag(s); 
    s(r+1:end) = 0; 
    X = u*diag(s)*v';
else
    error('typeR should be 1, 2 or 3.');
end