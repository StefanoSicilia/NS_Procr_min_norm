%% Script to compare all the approaches to solve NSPSDP with minimal norm 
% Perturbed Low-rank X

    %% INIZIALIZATION
    %Install
    addpath('./NSPSD_PSD_Procrustes_v1')
    addpath('./NSPSD_PSD_Procrustes_v1\utils')
    
    % Example
    n=200; 
    m=200;
    r=50;  
    
    % Parameters for reduced approach
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    options.rank=0;
    ranktol=1e-8;
    eta=1e-8;
    
    %% OUTPUT STRUCTURES
    nex=6;
    nsample=20;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: GILLIS ET AL. METHOD
    fprintf('1: ANFGM... ')
    i=1;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=procrustes_anly(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 2: FGM
    fprintf('2: FGM... ')
    i=2;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_FGM(B,X,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 3: FMINGD
    fprintf('3: FMINGD... ')
    i=3;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Min_GD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 4: CARDANO
    fprintf('4: CARDANO... ')
    i=4;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Cardano(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    % Low rank correction option
    options.rank=1;
    
    %% 5: MINGD LR
    fprintf('5: MINGD LR... ')
    i=5;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Min_GD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 6: CARDANO LR
    fprintf('6: CARDANO LR... ')
    i=6;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m)+eta*randn(n,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Cardano(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro')^2;
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% MEANS AND STD OF DATA
    RelErr=sqrt(Functional)/norm(B,'fro');
    
    MeanRelErr=mean(RelErr,2);
    MeanNorms=mean(Norms,2);
    MeanCPUtime=mean(CPUtime,2);
    MeanRanksSym=mean(RanksSym,2);
    MeanRanksSkew=mean(RanksSkew,2);
    
    StdRelErr=std(RelErr,0,2);
    StdNorms=std(Norms,0,2);
    StdCPUtime=std(CPUtime,0,2);
    StdRanksSym=std(RanksSym,0,2);
    StdRanksSkew=std(RanksSkew,0,2);
    
    %% TABLES
    close all
    rownames={'ANFGM','FGM','MINGD','CARD','MINGD_LR','CARD_LR'};    
    columnnames={'Fun','Norm_sol','Time','Rk(sym(A))','Rk(skew(A))'};
    
    % Means
    figure   
    Res={[5,1]};
    Res{1}=MeanRelErr;
    Res{2}=MeanNorms;
    Res{3}=MeanCPUtime;
    Res{4}=MeanRanksSym;
    Res{5}=MeanRanksSkew;
    T1=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames); 
    UT1=uitable('Data',T1{:,:},'ColumnName',columnnames,...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);
    title('Means') 

    % Standard Deviations
    figure   
    Res={[5,1]};
    Res{1}=StdRelErr;
    Res{2}=StdNorms;
    Res{3}=StdCPUtime;
    Res{4}=StdRanksSym;
    Res{5}=StdRanksSkew;
    T2=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames); 
    UT2=uitable('Data',T2{:,:},'ColumnName',columnnames,...
    'RowName',T2.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);    
    title('Stds')     
    
    % All data
    MeanRelErr=round(MeanRelErr,4);
    MeanNorms=round(MeanNorms,4);
    MeanCPUtime=round(MeanCPUtime,2);
    MeanRanksSym=round(MeanRanksSym,2);
    MeanRanksSkew=round(MeanRanksSkew,2);
    StdRelErr=round(StdRelErr,4);
    StdNorms=round(StdNorms,4);
    StdCPUtime=round(StdCPUtime,2);
    StdRanksSym=round(StdRanksSym,2);
    StdRanksSkew=round(StdRanksSkew,2);
    
    figure
    AllRelErr=strcat('$',num2str(MeanRelErr),' \pm ',num2str(StdRelErr),'$');
    Res{1}=cellstr(AllRelErr);
    AllNorms=strcat('$',num2str(MeanNorms),' \pm ',num2str(StdNorms),'$');
    Res{2}=cellstr(AllNorms);
    AllCPUtime=strcat('$',num2str(MeanCPUtime),' \pm ',num2str(StdCPUtime),'$');
    Res{3}=cellstr(AllCPUtime);
    AllRanksSym=strcat('$',num2str(MeanRanksSym),' \pm ',num2str(StdRanksSym),'$');
    Res{4}=cellstr(AllRanksSym);
    AllRanksSkew=strcat('$',num2str(MeanRanksSkew),' \pm ',num2str(StdRanksSkew),'$');
    Res{5}=cellstr(AllRanksSkew);
    T3=table(Res{1},Res{2},Res{3},Res{4},Res{5},'RowNames',rownames);
    UT3=uitable('Data',T3{:,:},'ColumnName',columnnames,...
    'RowName',T3.Properties.RowNames,'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);    
    title('All')
    
    % Use table generator website to produce latex output
    
    