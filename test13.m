%% Script for the figure in the paper, just for FGM and for n=2^9=512

    %% INIZIALIZATION
    %Install
    addpath('./NSPSD_PSD_Procrustes_v1')
    addpath('./NSPSD_PSD_Procrustes_v1\utils') 
    
    % Parameters for reduced approach
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    options.rank=1;
    ranktol=1e-8;
    
    %% OUTPUT STRUCTURES
    j=9;     % dimension n=2^j
    kmax=10; % number of the sample of random matrices
    A={[kmax,1]};
    CPUtime=zeros(kmax,1);
    
    %% FGM
    fprintf('FGM... ')
    for k=1:kmax
        n=2^j;
        m=n;
        r=n/2;
        rng(k)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        disp(['Number ',num2str(k),'.'])
        tic;
        A{k}=Procrustes_FGM(B,X,options); 
        CPUtime(k)=toc;
    end
    fprintf('Done!\n')
    
    MeanTime=mean(CPUtime);
    
    
    