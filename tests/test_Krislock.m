%% Script with the problem from Krislock et al. paper 

    %% Definition of the example
    n=3; 
    m=12;
    r=3;  
    load krislock_plush
    [U,D,V]=svd(X);
    X2=U(:,1:2)*D(1:2,1:2)*V(:,1:2)';
    
    %% Inizialization
    % Parameters for Fast Gradient Method
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    
    % Parameters for the ranks
    options.rank=1;
    ranktol=1e-8;
    
    % Output structures
    nex=8;
    nsample=1;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: ANFGM
    fprintf('1: ANFGM... ')
    i=1;
    tic;
    A{i}=Procrustes_ANFGM(X,B,options); 
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
    
    %% 3: MINGD
    fprintf('3: MINGD... ')
    i=3;
    tic;
    A{i}=Procrustes_MINGD(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 4: CARD
    fprintf('4: CARD... ')
    i=4;
    tic;
    A{i}=Procrustes_CARD(X,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')

    %% 5: ANFGM
    fprintf('5: ANFGM... ')
    i=5;
    tic;
    A{i}=Procrustes_ANFGM(X2,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 6: FGM
    fprintf('6: FGM... ')
    i=6;
    tic;
    A{i}=Procrustes_FGM(B,X2,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 7: MINGD
    fprintf('7: MINGD... ')
    i=7;
    tic;
    A{i}=Procrustes_MINGD(X2,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% 8: CARD
    fprintf('8: CARD... ')
    i=8;
    tic;
    A{i}=Procrustes_CARD(X2,B,options); 
    CPUtime(i)=toc;
    Functional(i)=norm(A{i}*X-B,'fro');
    Norms(i)=norm(A{i},'fro');
    RanksSym(i)=rank(symm(A{i}),ranktol);
    RanksSkew(i)=rank(skew(A{i}),ranktol);
    fprintf('Done!\n')
    
    %% Table of results
    RelErr=Functional/norm(B,'fro');
    rownames={'ANFGM','FGM','MINGD','CARD',...
        'ANFGMLR','FGMLR','MINGDLR','CARDLR'};    
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
    
    
    