% Compare CVX, FGM and AN-FGM on the data set from the paper 
% N. Krislock, J. Lang, J. Varah, D. K. Pai, and H.-P. Seidel, 
% Local compliance estimation via positive semidefinite constrained least 
% squares, IEEE Transactions on Robotics, 20 (2004), pp. 1007-1011 

clear all; clc; 
load jochen_data2
disp('Comparing CVX, AN-FGM and FGM on the NSPSD problem from');  
disp('the paper of Krislock et al. (IEEE Robotics, 2004).');
% CVX
optionscvx.nspsd = 1;
[Acvx,ecvx,tcvx,statcvx] = Procrustes_cvx(B,X,optionscvx);
ecvx = 100*norm( Acvx*X - B  )/norm(B,'fro'); 
fprintf('CVX    : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',ecvx,tcvx)
% AN-FGM
optionsANFGM.nspsd = 1;
optionsANFGM.timemax = tcvx;
optionsANFGM.delta = 1e-6;
optionsANFGM.maxiter = +Inf; 
[Aanfgm,eanfgm,tanfgm,stataf] = procrustes_anly(X,B,optionsANFGM);
eanfgm = 100*norm( Aanfgm*X - B  )/norm(B,'fro'); 
tanfgm(end); 
fprintf('AN-FGM : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',eanfgm,tanfgm(end))
% FGM
optionsFGM.nspsd = 1;
optionsFGM.timemax = tcvx;
optionsFGM.delta = 1e-6;
[Afgm,efgm,tfgm] = Procrustes_FGM(B,X,optionsFGM);
efgm = 100*norm( Afgm*X - B  )/norm(B,'fro'); 
fprintf('FGM    : ||AX-B||_F/||B||_F = %2.4f%%, time = %2.2fs.\n',eanfgm,tfgm(end))