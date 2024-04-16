%% Script with the problem from Krislock et al. paper 

    %% INIZIALIZATION
    %Install
    addpath('./NSPSD_PSD_Procrustes_v1')
    addpath('./NSPSD_PSD_Procrustes_v1\utils')
    
    % Example
    n=3; 
    m=12;
    r=3;  
    X=[-0.32 -0.33 -0.36 -0.30 -0.32 -0.34 -0.24 -0.21 -0.33 -0.25 -0.22 -0.31;
        0.03 -0.02  0.08  0.03 -0.00  0.07  0.07 -0.01  0.16  0.09  0.00  0.15;
        0.06  0.06  0.06  0.05  0.07  0.05  0.05  0.02  0.10  0.06  0.03  0.09];
    B=[-1.43 -1.40 -1.38 -1.43 -1.40 -1.37 -1.43 -1.40 -1.38 -1.43 -1.40 -1.37;
        0.15 -0.31  0.44  0.14 -0.31  0.43  0.16 -0.32  0.42  0.15 -0.33  0.42;
       -0.44 -0.42 -0.42 -0.44 -0.42 -0.42 -0.43 -0.42 -0.43 -0.44 -0.42 -0.44];
    
    % Parameters for reduced approach
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    options.rank=1;
    ranktol=1e-8;
    
    %% OUTPUT STRUCTURES
    nex=4;
    nsample=1;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: GILLIS ET AL. METHOD
    fprintf('1: ANFGM... ')
    i=1;
    tic;
    A{i}=procrustes_anly(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 2: FGM
    fprintf('2: FGM... ')
    i=2;
    tic;
    A{i}=Procrustes_FGM(B,X,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 3: FMINGD
    fprintf('3: FMINGD... ')
    i=3;
    tic;
    A{i}=Procrustes_Min_GD(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 4: CARDANO
    fprintf('4: CARDANO... ')
    i=4;
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
    rownames={'ANFGM','FGM','MINGD','CARDANO'};    
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
    
    