%% Script to compare the scale of computational time for random examples

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
    nex=2;
    nsample=11;
    A={[nex,nsample]};
    Functional=zeros(nex,nsample);
    CPUtime=zeros(nex,nsample);
    Norms=zeros(nex,nsample);
    RanksSym=zeros(nex,nsample);
    RanksSkew=zeros(nex,nsample);
    
    %% 1: FMINGD
    fprintf('1: FMINGD... ')
    i=1;
    for j=1:nsample
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Min_GD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% 2: CARDANO
    fprintf('2: CARDANO... ')
    i=2;
    for j=1:nsample
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        X=randn(n,r)*randn(r,m);
        B=randn(n,m);
        tic;
        A{i,j}=Procrustes_Cardano(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% SCALING OF COMPUTATIONAL TIME
    close all
    vecj=2:nsample+1;
    vecn=2.^vecj;
    figure(1)
    semilogy(vecj,CPUtime(1,:),'b-o','LineWidth',1.5)
    hold on
    semilogy(vecj,CPUtime(2,:),'r-s','LineWidth',1.5)
    legend('MINGD','CARD','Location','northwest')
    
    figure(2)
    plot(vecn,CPUtime(1,:),'b-o','LineWidth',1.5)
    hold on
    plot(vecn,CPUtime(2,:),'r-s','LineWidth',1.5)
    legend('MINGD','CARD','Location','northwest')
    