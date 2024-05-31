% Initialization for the PSD Procrustes problem 
% 
% min_{A or (A'+A) PSD} ||AX-B||_F. 
% 
% It compares the following two initializations: 
% (1) The projection onto the PSD cone of the optimal solution of the
% unconstrained problem, and 
% (2) The diagonal optimal solution, 
% and keeps the best among the two. 
% 
% If the input choicei is specified, 
% choicei == 1 => initialization (1).
% choicei == 2 => initialization (2). 
% 
% See section 4.2 in the paper 
% A semi-analytical approach for the positive semidefinite Procrustes 
% problem, Nicolas Gillis and Punit Sharma, 2016. 

function A0 = init_procrustes(B,X,choicei,nspsd);  
 
if nargin <= 2
    choicei = 0;
end
if nargin <= 3
    nspsd = 1;
end
n = size(X,1); 
% Projection of the optimal solution of the unconstrained problem
if choicei == 1 || choicei == 0
    A1 = projectNSPSDorPSD((X'\B')',nspsd) ;
    e1 = norm(A1*X-B,'fro'); 
    if choicei == 1
        A0 = A1; 
        return; 
    end
end
% Optimal diagonal solution
if choicei == 2 || choicei == 0
    A2 = zeros(n); 
    for i = 1 : n
        A2(i,i) = max( 0 , (X(i,:)*B(i,:)') / (norm( X(i,:) )^2+1e-6) );
    end
    e2 = norm(A2*X-B,'fro'); 
    if choicei == 2
        A0 = A2; 
        return; 
    end
end
if choicei == 0
    if e1 < e2 
        A0 = A1;
    else
        A0 = A2;
    end
end