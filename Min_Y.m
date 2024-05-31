function Y=Min_Y(D,W,Z)
%% Min_Y: 
% Solves, through gradient descent method the optimization problem
% min_Y 1/16*norm(Y*D^(-1)*Y','fro')^2+norm(Y*W'-Z,'fro')^2,
% by using the gradient G=1/4*Y*D^(-1)*Y'*Y*D^(-1)+2*(Y-Z*W).
% It implements the line search back-tracking method.

    Yold=Z*W;
    sqD=sqrtm(inv(D));
    YD=Yold*sqD;
    normYDY=norm(YD'*YD,'fro');
    p=8*norm(Yold,'fro')^2/normYDY^2;
    alpha=Cardano(p,-p);
    Yold=alpha*Yold;
    maxit=200;
    tol=1e-8;
    theta1=2;
    theta2=1.1;
    
    % First iteration
    fold=norm(Yold*W'-Z,'fro')^2+1/16*normYDY^2;
    Gold=0.25*Yold/D*(Yold'*Yold)/D+2*(Yold-Z*W);
    gamma=(norm(Yold,'fro')/norm(Gold,'fro'))^2;
    fnew=fold+1;
    while fnew>=fold && gamma>1e-16
        Ynew=Yold-gamma*Gold;
        YD=Ynew*sqD;
        fnew=norm(Ynew*W'-Z,'fro')^2+1/16*norm(YD'*YD,'fro')^2;
        gamma=gamma/theta1;
    end
    gamma=gamma*theta2;
    
    % Main Iterations
    j=1;
    err=tol+1;
    while j<maxit && err>=tol
        fnew=fold+1;
        G=0.25*Yold/D*(Yold'*Yold)/D+2*(Yold-Z*W);
        k=1;
        while fnew>=fold && k<maxit
            Ynew=Yold-gamma*G;
            YD=Ynew*sqD;
            fnew=norm(Ynew*W'-Z,'fro')^2+1/16*norm(YD'*YD,'fro')^2;
            gamma=gamma/theta1;
            k=k+1;
        end
        gamma=gamma*theta2;
        j=j+1;
        err=norm(Ynew-Yold,'fro')/norm(Ynew,'fro');
        Yold=Ynew;
    end
    Y=Ynew;
    %disp(['The number of iterations is ',num2str(j),'.'])
end