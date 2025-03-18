%% Scaling of ANFGM, FGM, MINGD and CARD in solving NSPSD
% Example in Section 5.5: X=J_nJ_n' and B=H.

    %% Inizialization    
    % Parameters for Fast Gradient method
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6; 
    
    % Parameters for the ranks
    options.rank=1;
    ranktol=1e-8;
    
    % Output structures
    nex=4;
    jmax=11; 
    jmaxFGM=8; 
    A={[nex,jmax]};
    Functional=zeros(nex,jmax);
    CPUtime=zeros(nex,jmax);
    Norms=zeros(nex,jmax);
    RanksSym=zeros(nex,jmax);
    RanksSkew=zeros(nex,jmax);
    
    %% 1: ANFGM
    fprintf('1: ANFGM... ')
    i=1;
    for j=1:jmax
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        J=triu(ones(n));
        J=J(:,1:r);
        H=triu(ones(m));
        H=H(:,1:r);
        X=J*H';
        B=toeplitz(1:n,[1;zeros(m-1,1)]);
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
    for j=1:jmaxFGM
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        J=triu(ones(n));
        J=J(:,1:r);
        H=triu(ones(m));
        H=H(:,1:r);
        X=J*H';
        B=toeplitz(1:n,[1;zeros(m-1,1)]);
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
    for j=1:jmax
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        J=triu(ones(n));
        J=J(:,1:r);
        H=triu(ones(m));
        H=H(:,1:r);
        X=J*H';
        B=toeplitz(1:n,[1;zeros(m-1,1)]);
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
    for j=1:jmax
        n=2^(j+1);
        m=n;
        r=n/2;
        rng(j)
        J=triu(ones(n));
        J=J(:,1:r);
        H=triu(ones(m));
        H=H(:,1:r);
        X=J*H';
        B=toeplitz(1:n,[1;zeros(m-1,1)]);
        tic;
        A{i,j}=Procrustes_CARD(X,B,options); 
        CPUtime(i,j)=toc;
        Functional(i,j)=norm(A{i,j}*X-B,'fro');
        Norms(i,j)=norm(A{i,j},'fro');
        RanksSym(i,j)=rank(symm(A{i,j}),ranktol);
        RanksSkew(i,j)=rank(skew(A{i,j}),ranktol);
    end
    fprintf('Done!\n')
    
    %% SCALING OF COMPUTATIONAL TIME
    close all
    figure(1)
    shift=1;
    vecj=(1+shift):(jmax+shift);
    vecjFGM=vecj(1:jmaxFGM);
    l=length(vecjFGM);
    figure(1)
    semilogy(vecj,CPUtime(1,:),'b--o','LineWidth',1.5)
    hold on
    semilogy(vecjFGM,CPUtime(2,1:l),'r-.s','LineWidth',1.5)
    hold on
    semilogy(vecj,CPUtime(3,:),'m-^','LineWidth',1.5)
    hold on
    semilogy(vecj,CPUtime(4,:),'c:d','LineWidth',1.5)
    xlabel('$j$','Interpreter','Latex')
    ylabel('time (s.)','Interpreter','Latex')
    legend('ANFGM','FGM','MINGD','CARD','Location','northwest')
    fontsize(14,'points')
    grid on

    figure(2)
    plot(vecj,CPUtime(1,:),'b--o','LineWidth',1.5)
    hold on
    plot(vecjFGM,CPUtime(2,1:l),'r-.s','LineWidth',1.5)
    hold on
    plot(vecj,CPUtime(3,:),'m-^','LineWidth',1.5)
    hold on
    plot(vecj,CPUtime(4,:),'c:d','LineWidth',1.5)
    xlabel('$j$','Interpreter','Latex')
    ylabel('time (s.)','Interpreter','Latex')
    legend('ANFGM','FGM','MINGD','CARD','Location','northwest')
    fontsize(14,'points')