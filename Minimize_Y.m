function Y=Minimize_Y(D,W,Z)
%% Minimize_Y: 
% Solves, through gradient descent method the optimization problem
% min_Y norm(Y*D^(-1)*Y','fro')^2+norm(Y*W'-Z,'fro')^2,
% by using the gradient G=4*Y*D^(-1)*Y'*Y*D^(-1)+Y-Z*W.

    rng(1)
    %Yold=zeros(size(Z,1),size(D,1));
    %Yold=randn(size(Z,1),size(D,1));
    Yold=Z*W;
    rho=norm(Yold,'fro');
    lambda=norm(inv(D),'fro');
    L=(0.75*rho^2*lambda^2+2);
    maxit=500000;
    tol=1e-8;
    j=0;
    err=tol+1;
    while j<maxit && err>=tol
        G=0.25*Yold/D*(Yold'*Yold)/D+2*(Yold-Z*W);
        Ynew=Yold-1/L*G;
        Ynew=min(rho/norm(Ynew,'fro'),1)*Ynew;
        j=j+1;
        err=norm(Ynew-Yold,'fro');%/norm(Ynew,'fro');
        Yold=Ynew;
    end
    Y=Ynew;
    disp(['The number of iterations is ',num2str(j),'.'])
end