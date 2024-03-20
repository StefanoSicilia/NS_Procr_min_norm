function x=Cardano(p,q)
%% Cardano:
% Solves the equation x^3+px+q=0 in the case that q^2/4+p^3/27>0.

    Delta=q^2/4+p^3/27;
    if Delta<0
        error('The Cardano formula does not work in this case.')
    end
    x=nthroot(-q/2+sqrt(Delta),3)+nthroot(-q/2-sqrt(Delta),3);

end