clc;clear all;close all;


K = 0;
format long
Rho = [0.1,10,20,50];
F_t = [];
R_t = [];
for n  = 1:length(Rho)
    tic
    disp('========================================================')
   
    file = 'case30';
    % file = ['matpower5.1\case',num2str(n),'.m'];
%     if ~exist(file,'file')
%         disp('File not exists')
%         continue;
%     end
   
    try
        eval(['opf=',file,';']);
    catch ME
        disp('Fail to find OPF object')
        continue;
    end
    
    try
        Ybus = makeYbus(opf);
    catch ME
        disp('Forming Y bus failed')
        continue;
    end
    [N,temp] = size(Ybus);

    K = K+1; % number of tested cases
    Y = full(Ybus);
    
    [Phi,Psi,M,T] = Matrix_SDP(Y);
    
    bus = opf.bus;
    %[N,temp] = size(bus);
    
    gen = opf.gen;
    genBus = gen(:,1); % number of buses with generatros;
    
    
    baseMVA = 100;
    
    %% constraints
    DemandP = 1*bus(:,3);
    DemandQ = 1*bus(:,4);
    GenPmax = zeros(N,1);
    GenPmax(genBus) = gen(:,9);
    GenPmin = zeros(N,1);
    GenPmin(genBus) = gen(:,10);
    
    GenQmax = zeros(N,1);
    GenQmax(genBus) = gen(:,4);
    GenQmin = zeros(N,1);
    GenQmin(genBus) = gen(:,5);
    
    % Pmin Pmax
    Pmax = (GenPmax-DemandP)/baseMVA;
    Pmin = (GenPmin-DemandP)/baseMVA;
    Qmax = (GenQmax-DemandQ)/baseMVA;
    Qmin = (GenQmin-DemandQ)/baseMVA;
    
    
    
    % Vmin Vmax
    Vmin = bus(:,13);
    Vmax = bus(:,12);
    
    
    
    %% Cost functions
    % cost function
    %C = eye(N); %min \|V\|_2
    %C = (Y+Y')/2; % min power loss
    %C = sum_k ck*Phik;  %min production costs.
    MinCost = 0;
    genCost = opf.gencost;
    C = 0;
    for k = 1:length(genBus)
        BusNum = genBus(k);
        if genCost(k,1) == 2
            c = genCost(k,end-1);
%             for n = 1:genCost(k,4)
%                c = c+genCost(k,end-1); 
%             end
            C = C+c*Phi(:,(BusNum-1)*N+1:BusNum*N);
            MinCost = MinCost + c*DemandP(BusNum);
        else
            disp('Type 1 Found; the cost function needs to be fixed')
            
            return;
        end
        
        
    end
    
    %% SDP formulation
    
    %disp('Matrix preparation ready')
    %pause
    cvx_begin sdp quiet
    variable W(N,N) complex semidefinite
    minimize(trace(W*C))
    subject to
    for k = 1:N
        trace(W*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
        trace(W*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
        trace(W*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
        trace(W*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
        W(k,k) >= Vmin(k).^2
        W(k,k) <= Vmax(k).^2
    end
    cvx_end
    cvx_clear
    if strcmp(cvx_status,'Infeasible')
        disp('SDP relaxation is infeasible')
        continue;
    end
    if rank(W) == 1
        disp('SDP is exact')
        continue;
    else
        [X,f_t,r_t] = RLRA4SDP(W,Phi,Psi,Pmax,Pmin,Qmax,Qmin,Vmax,Vmin,300,Rho(n),file);
        F_t = [F_t,f_t];
        R_t = [R_t,r_t];
        [s,v,d] = svd(W);
        v = diag(v);
        disp('optimal of relaxed problem: ')
        rank(W);
        disp(['objective value: ',num2str(trace(W*C)*baseMVA+MinCost)])
        perc = 1-v(1)/sum(v);
        if perc <= 10^(-7)
            disp('feasible')
        else
            disp('infeasible')
        end
        disp(['spetrum of the solution: ',num2str(v(1)/sum(v))])
        
        
        disp('Solution of RLRA: ')
        [s,v,d] = svd(X);
        v = diag(v);
        perc = 1-v(1)/sum(v);
        if perc <= 10^(-7)
            disp('feasible')
        else
            disp('infeasible')
        end
        disp(['spetrum of the solution: ',num2str(v(1)/sum(v))])
        disp(['new objective value: ',num2str(trace(X*C)*baseMVA+MinCost)])
        
        disp(['distance: ',num2str(norm(X-W,'fro')/norm(W,'fro'))])
        
    end
    
    toc
    
    
end
save('Bus30.mat')
format short