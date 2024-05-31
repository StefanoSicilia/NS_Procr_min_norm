%% Scaling of ANFGM, FGM, MINGD and CARD in solving NSPSD
% Example in Section 5.4: average of a sample of 10 random X and B.

    %% Inizialization 
    % Parameters for Fast Gradient Method
    options.nspsd=1;
    options.timemax=10; 
    options.delta=1e-6;
    
    % Parameters for the ranks
    options.rank=1;
    
    %% OUTPUT STRUCTURES
    shift=1; % shift for j: the values used are from 1+shift to jmax+shift
    imax=4;  % number of methods
    jmax=11;  % number of maximum dimension n=2^(jmax+shift) 
    jmaxFGM=8; % number of maximum dimension for FGM 
    kmax=10; % number of the sample of random matrices
    A={[imax,jmax,kmax]};
    CPUtime=zeros(imax,jmax,kmax);
    
    %% 1: GILLIS ET AL. METHOD
    fprintf('1: ANFGM... ')
    i=1;
    for j=1:jmax
        for k=1:kmax
            n=2^(j+shift);
            m=n;
            r=n/2;
            rng(k)
            X=randn(n,r)*randn(r,m);
            B=randn(n,m);
            tic;
            A{i,j,k}=Procrustes_ANFGM(X,B,options); 
            CPUtime(i,j,k)=toc;
        end
    end
    fprintf('Done!\n')
    
    %% 2:FGM
    fprintf('2: FGM... ')
    i=2;
    for j=1:jmaxFGM
        for k=1:kmax
            n=2^(j+shift);
            m=n;
            r=n/2;
            rng(k)
            X=randn(n,r)*randn(r,m);
            B=randn(n,m);
            tic;
            A{i,j,k}=Procrustes_FGM(B,X,options); 
            CPUtime(i,j,k)=toc;
        end
    end
    fprintf('Done!\n')
    
    %% 3: MINGD
    fprintf('3: MINGD... ')
    i=3;
    for j=1:jmax
        for k=1:kmax
            n=2^(j+shift);
            m=n;
            r=n/2;
            rng(k)
            X=randn(n,r)*randn(r,m);
            B=randn(n,m);
            tic;
            A{i,j,k}=Procrustes_MINGD(X,B,options); 
            CPUtime(i,j,k)=toc;
        end
    end
    fprintf('Done!\n')
    
    %% 4: CARD
    fprintf('4: CARD... ')
    i=4;
    for j=1:jmax
        for k=1:kmax
            n=2^(j+shift);
            m=n;
            r=n/2;
            rng(k)
            X=randn(n,r)*randn(r,m);
            B=randn(n,m);
            tic;
            A{i,j,k}=Procrustes_CARD(X,B,options); 
            CPUtime(i,j,k)=toc;
        end
    end
    fprintf('Done!\n')

    MeanCPUtime=mean(CPUtime,3);
    
    %% Scaling on cpu time
    close all
    vecj=(1+shift):(jmax+shift);
    vecjFGM=vecj(1:jmaxFGM);
    l=length(vecjFGM);
    figure(1)
    semilogy(vecj,MeanCPUtime(1,:),'b-o','LineWidth',1.5)
    hold on
    semilogy(vecjFGM,MeanCPUtime(2,1:l),'r-s','LineWidth',1.5)
    hold on
    semilogy(vecj,MeanCPUtime(3,:),'m-^','LineWidth',1.5)
    hold on
    semilogy(vecj,MeanCPUtime(4,:),'c-d','LineWidth',1.5)
    xlabel('$j$','Interpreter','Latex')
    ylabel('Time (sec.)')
    legend('ANFGM','FGM','MINGD','CARD','Location','northwest')
    
    figure(2)
    plot(vecj,MeanCPUtime(1,:),'b-o','LineWidth',1.5)
    hold on
    plot(vecjFGM,MeanCPUtime(2,1:l),'r-s','LineWidth',1.5)
    hold on
    plot(vecj,MeanCPUtime(3,:),'m-^','LineWidth',1.5)
    hold on
    plot(vecj,MeanCPUtime(4,:),'c-d','LineWidth',1.5)
    xlabel('$j$','Interpreter','Latex')
    ylabel('Time (sec.)')
    legend('ANFGM','FGM','MINGD','CARD','Location','northwest')
    