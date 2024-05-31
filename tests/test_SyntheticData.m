%% Comparison of ANFGM, FGM, MINGD and CARD in solving NSPSD
% Synthetic data in Section 5.1: random low-rank X and random B. 
    
    %% Example dimensions
    type=1; % refers to the Table number in the paper
    switch type
        case 1
            n=50;
            m=70;
            r=20; 
        case 2
            n=100;
            m=100;
            r=40;
        case 3
            n=200;
            m=200;
            r=50;
    end
    
    %% Inizialization
    % Parameters for Fast Gradient Method
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    
    % Parameters for the ranks
    options.rank=1;
    ranktol=1e-8;    
    
    % Output structures
    nex=4;
    nsample=20;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: ANFGM
    fprintf('1: ANFGM... ')
    i=1;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_ANFGM(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
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
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_FGM(B,X,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 3: MINGD
    fprintf('3: MINGD... ')
    i=3; 
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_MINGD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 4: CARD
    fprintf('4: CARD... ')
    i=4;
    for j=1:nsample
        rng(j)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_CARD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% Means and standard deviations of data
    RelErr=Functional/norm(B,'fro');
    
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
    
    %% Tables  
    close all
    rownames={'ANFGM','FGM','MINGD','CARDANO'};
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
    title('MEANs') 

    % Standard deviations
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
    title('STDs')     
    