%% Script with the problem of large dimension (RANDOM) 

    %% INIZIALIZATION
    %Install
    addpath('./NSPSD_PSD_Procrustes_v1')
    addpath('./NSPSD_PSD_Procrustes_v1\utils')
    
    % Example
    n=2048; 
    m=n;
    r=ceil(n/2);
    rng(1)
    X=randn(n,r)*randn(r,m);
    B=randn(n,m);
    
    % Parameters for reduced approach
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    options.rank=1;
    ranktol=1e-8;
    
    %% OUTPUT STRUCTURES
    nex=2;
    nsample=1;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: FMINGD
    fprintf('1: FMINGD... ')
    i=1;
    tic;
    A{i}=Procrustes_Min_GD(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 2: CARDANO
    fprintf('2: CARDANO... ')
    i=2;
    tic;
    A{i}=Procrustes_Cardano(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% TABLE OF RESULTS
    RelErr=Functional/norm(B,'fro');
    rownames={'MINGD','CARDANO'};    
    close all
    columnnames={'Fun','Norm_sol','Time','Rk(sym(A))','Rk(skew(A))'};
    figure   
    Res={[5,1]};
    Res{1}=RelErr;
    Res{2}=Norms;
    Res{3}=CPUtime;
    Res{4}=RanksSym;
    Res{5}=RanksSkew;
    T1=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames); 
    UT1=uitable('Data',T1{:,:},'ColumnName',columnnames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

    % Use table generator website to produce latex output
    
    