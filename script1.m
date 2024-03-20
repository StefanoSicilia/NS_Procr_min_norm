%% Script to compare all the approaches to solve NSPSDP with minimal norm 
% Single matrix

    %% INIZIALIZATION
    %Install
    addpath('./NSPSD_PSD_Procrustes_v1')
    
    % Example
    n=9; 
    m=8;
    r=4;
    rng(1)
    X=randn(n,r)*randn(r,m);
    B=rand(n,m);    
    
    % Parameters for reduced approach
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    ranktol=1e-8;
    
    %% OUTPUT STRUCTURES
    nex=6;
    A={[nex,1]};
    Functional=zeros(nex,1);
    CPUtime=zeros(nex,1);
    Norms=zeros(nex,1);
    RanksSym=zeros(nex,1);
    RanksSkew=zeros(nex,1);
    
    %% 1: GILLIS ET AL. METHOD
    fprintf('1: ANFGM... ')
    i=1;
    tic;
    A{i}=procrustes_anly(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 2: PLR
    fprintf('2: PLR... ')
    i=2;
    tic;
    A{i}=Procrustes_LowRank(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 3: FGM
    fprintf('3: FGM... ')
    i=3;
    tic;
    options.timemax = 1;
    A{i}=Procrustes_FGM(B,X,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 4: FMINUNC
    fprintf('4: FMINUNC... ')
    i=4;
    tic;
    A{i}=Procrustes_Fminunc(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 5: FMINGD
    fprintf('5: FMINGD... ')
    i=5;
    tic;
    A{i}=Procrustes_Min_GD(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 6: CARDANO
    fprintf('6: CARDANO... ')
    i=6;
    tic;
    A{i}=Procrustes_Cardano(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro')^2;
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% TABLE
    rownames={'ANFGM','PLR','FGM','FMINUNC','MINGD','CARDANO'};    
    close all
    columnnames={'Fun_fr','Norm_A_fr','Time_fr','Rk(sym(A))_fr',...
        'Rk(skew(A))_fr','Soldiff'};
    figure    
    Res={[5,1]};
    Res{1}=Functional;
    Res{2}=Norms;
    Res{3}=CPUtime;
    Res{4}=RanksSym;
    Res{5}=RanksSkew;
    T1=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames); 
    UT1=uitable('Data',T1{:,:},'ColumnName',columnnames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);
    title('Results')     
    
    