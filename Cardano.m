function x=Cardano(p,q)
%% Cardano:
% Solves the equation x^3+px+q=0 in the case that Delta=q^2/4+p^3/27>0.

    sqDelta=sqrt(q^2/4+p^3/27);
    x=nthroot(-q/2+sqDelta,3)+nthroot(-q/2-sqDelta,3);

end