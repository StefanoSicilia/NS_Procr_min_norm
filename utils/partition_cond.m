% Given a positive vector x of size n, sorted in nondecreasing order, 
% it returns an index set kb such that vb=max(x(kb)/x(1),x(n)/x(kb+1)) is 
% minimized. The vector f contains the values of that criterion for all k.

function [kb,vb,f] = partition_cond(x) 

n = length(x); 

for k = 1 : n
    if k == 1
        f(k) = x(n)/x(2); 
    elseif k == n
        f(k) = x(n-1)/x(1); 
    else
        f(k) = max( x(k)/x(1)  , x(n)/x(k+1) );
    end
end

[vb,kb] = min(f);  