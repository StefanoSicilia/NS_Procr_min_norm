function Y=Min_Y(Lambda,W,Z)
%% Min_Y: 
% Solves, through gradient descent method the optimization problem
% min_Y 1/16*norm(Y*Lambda^(-1)*Y','fro')^2+norm(Y*W'-Z,'fro')^2,
% by using the gradient G=1/4*Y*Lambda^(-1)*Y'*Y*Lambda^(-1)+2*(Y-Z*W).
% It implements the line search back-tracking method of Algorithm 2, by
% verifying the Armijo condition.

    Yold=Z*W;
    D=diag(diag(ones(size(Lambda,1))./diag(Lambda)));
    sqinvD=sqrtm(D);
    YD=Yold*sqinvD;
    normYDY=norm(YD'*YD,'fro');
    p=8*norm(Yold,'fro')^2/normYDY^2;
    alpha=Cardano(p,-p);
    Yold=alpha*Yold;
    YD=Yold*sqinvD;
    normYDY=norm(YD'*YD,'fro');
    maxit=500;
    tol=1e-8;
    theta1=1.5;
    c=1e-4;
    jmax=50;
    
    % First values
    fold=norm(Yold*W'-Z,'fro')^2+(1/16)*normYDY^2;
    G=0.25*Yold*D*(Yold'*Yold)*D+2*(Yold-Z*W);
    g=norm(G,'fro');
    relerr=1;
    
    % Main Iterations
    i=1;
    while i<maxit && g>=tol && relerr>=1e-15
        gamma=1e-1*norm(Yold,'fro')/g;
        for j=1:jmax 
            Ynew=Yold-gamma*G;
            YDnew=Ynew*sqinvD;
            fnew=norm(Ynew*W'-Z,'fro')^2+(1/16)*norm(YDnew'*YDnew,'fro')^2;
            if fnew<=fold-c*gamma*g^2
                break
            else
                gamma=gamma/theta1;
            end
        end
        if j==jmax
            disp('Armijo failed.')
            i=maxit;
        end
        i=i+1;
        G=0.25*Ynew*D*(Ynew'*Ynew)*D+2*(Ynew-Z*W);
        g=norm(G,'fro');
        relerr=abs(fold-fnew)/fold;
        Yold=Ynew;
        fold=fnew;
    end
    Y=Yold;
end
