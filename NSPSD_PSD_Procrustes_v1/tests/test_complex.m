% Compare CVX, FGM and AN-FGM on the complex NSPSD problem from the paper 
% J. Kiskiras and G. Halikias, A note on the complex semi-definite matrix 
% procrustes problem, Numer. Linear Algebra Appl., 14 (2007), pp. 485-502

clear all; clc; 

datacomplexNLA; 

disp('Comparing CVX, AN-FGM and FGM on the complex NSPSD problem from');  
disp('the paper of Kiskiras and Halikias (NLA, 2007).');

% CVX
optionscvx.nspsd = 1;
optionscvx.complex = 1;
[Acvx,ecvx,tcvx,statcvx] = Procrustes_cvx(B,X,optionscvx);
ecvx = 100*norm( Acvx*X - B  )/norm(B,'fro'); 
fprintf('CVX    : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',ecvx,tcvx)
% AN-FGM
optionsANFGM.nspsd = 1;
optionsANFGM.timemax = tcvx;
optionsANFGM.delta = 1e-6;
optionsANFGM.complex = 1; 
optionsANFGM.maxiter = +Inf; 
[Aanfgm,eanfgm,tanfgm,stataf] = procrustes_anly(X,B,optionsANFGM);
eanfgm = 100*norm( Aanfgm*X - B  )/norm(B,'fro'); 
tanfgm(end); 
fprintf('AN-FGM : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',eanfgm,tanfgm(end))
% FGM
optionsFGM.nspsd = 1;
optionsFGM.timemax = tcvx;
optionsFGM.delta = 1e-6;
optionsFGM.complex = 1; 
[Afgm,efgm,tfgm] = Procrustes_FGM(B,X,optionsFGM);
efgm = 100*norm( Afgm*X - B  )/norm(B,'fro'); 
fprintf('FGM    : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',eanfgm,tfgm(end))