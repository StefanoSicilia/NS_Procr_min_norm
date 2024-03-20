% This file allows you to run exactly the same experiment on synthetic data
% sets as in 
% On the non-symmetric semidefinite Procrustes problem, 
% by Mohit Kumar Baghel, Nicolas Gillis and Punit Sharma, 2021. 
% It compares CVX/SDPT3, FGM and AN-FGM.   

clc; clear all; close all; 

% Niter is the number of generated matrices for each case
% In the paper, Niter = 20. 
Niter = 2; 
% parameter for the size (m,n): max(m,n)=sizemn
% In the paper, sizemn = 60. 
sizemn = 60; % originally 20 
disp('*****************************************************************'); 
disp('Running the experiment from Section 6.1:'); 
disp('Comparaison of CVX/SDPT3, FGM and AN-FGM on synthetic data set.'); 
disp('*****************************************************************'); 
disp('Parameters:'); 
fprintf('Number of generated matrices for each case: %2.0f \n', Niter); 
fprintf('Size of the matrices: max(m,n) = %2.0f \n', sizemn); 
% Save results (errors and timings) in the following matrices: 
ecvxm = zeros(3,3,Niter); 
eanfgmm = zeros(3,3,Niter); 
tanfgmm = zeros(3,3,Niter); 
efgmm = zeros(3,3,Niter); 
tfgmm = zeros(3,3,Niter); 
tcvxm = zeros(3,3,Niter); 
statuscvx = zeros(3,3,3); 
% NSPSD PROCRUSTES SOLVED : 
optionscvx.nspsd = 1;  
optionsFGM.nspsd = 1;  
optionsANFGM.nspsd = 1;
% Main loop 
for typeR = 1 : 3
    for dimi = 1 : 3
        if dimi == 1
            n = sizemn; 
            m = n;
        elseif dimi == 2
            n = sizemn/2; 
            m = 2*n;
        elseif dimi == 3
            n = sizemn; 
            m = n/2;
        end
        disp('*****************************************************************'); 
        fprintf('Sarting case typeR = %1.0f, n = %2.0f, m = %2.0f.\n', typeR, n, m);
        for nit = 1 : Niter
            % Generate X and B
            B = randn(n,m);
            X = generateX(n,m,typeR);
            
            % CVX 
            [Acvx,ecvx,tcvx,statcvx] = Procrustes_cvx(B,X,optionscvx); 
            ecvx = 100*norm( Acvx*X - B  )/norm(B,'fro'); 
            ecvxm(typeR,dimi,nit) = ecvx; 
            tcvxm(typeR,dimi,nit) = tcvx; 
            % status for cvx: 
            if strcmp(statcvx , 'Solved')
                statuscvx(typeR,dimi,1) = statuscvx(typeR,dimi,1)+1; 
            elseif strcmp(statcvx , 'Inaccurate/Solved')
                statuscvx(typeR,dimi,2) = statuscvx(typeR,dimi,2)+1; 
            elseif strcmp(statcvx , 'Failed')
                statuscvx(typeR,dimi,3) = statuscvx(typeR,dimi,3)+1; 
            end
            % AN-FGM
            optionsANFGM.timemax = tcvx; 
            optionsANFGM.delta = 1e-6; 
            [Aanfgm,eanfgm,tanfgm,stataf] = procrustes_anly(X,B,optionsANFGM); 
            eanfgmm(typeR,dimi,nit) =  100*norm( Aanfgm*X - B  )/norm(B,'fro');
            tanfgmm(typeR,dimi,nit) = tanfgm(end); 
            % FGM
            optionsFGM.timemax = tcvx;  
            optionsFGM.delta = 1e-6; 
            [Afgm,efgm,tfgm] = Procrustes_FGM(B,X,optionsFGM);  
            efgmm(typeR,dimi,nit) = 100*norm( Afgm*X - B  )/norm(B,'fro');
            tfgmm(typeR,dimi,nit) = tfgm(end); 
            fprintf('%2.0f...',nit); 
            if mod(nit,10) == 0
                fprintf('\n');       
            end
        end
        fprintf('\n');
        fprintf('Status of cvx: solved (%2.0f), solved/innacurate (%2.0f), failed (%2.0f)', statuscvx(typeR,dimi,1) , statuscvx(typeR,dimi,2) , statuscvx(typeR,dimi,3));
        fprintf('\n');
    end
end

fprintf('\n'); 
disp('******************* DONE *************************'); 

TEST_paper_synth_display