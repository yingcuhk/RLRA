clc;clear all;close all;

inexact = {};
notsolved = {};
files = dir('matpower5.1\case*.m');
files = {files.name}';
OKfiles = {};
[N,temp] = size(files);
K = 0;
for n  = 1:N
    file = files{n};
    file = file(1:end-2)
    % file = ['matpower5.1\case',num2str(n),'.m'];
    
    if ~exist(file,'file')
        continue;
    end
    
    try
        eval(['OPF=',file,';']);
    catch ME
        continue;
    end
    % OPF = case30;
    try
        Ybus = makeYbus(OPF);
    catch ME
        continue;
    end
    [N,temp] = size(Ybus);
    if N > 500
        continue;
    end
    K = K+1; % number of tested cases
    Y = full(Ybus);
    
    [Phi,Psi,M,T] = Matrix_SDP(Y);
    
    bus = OPF.bus;
    %[N,temp] = size(bus);
    
    gen = OPF.gen;
    genBus = gen(:,1); % number of buses with generatros;
    
    
    baseMVA = OPF.baseMVA;
    
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
    C = eye(N); %min \|V\|_2
    %C = (Y+Y')/2; % min power loss
    %C = sum_k ck*Phik;  %min production costs.
    
    
    %% SDP formulation
    
    disp('Matrix preparation ready')
    %pause
    cvx_begin sdp quiet
    variable W(N,N) complex semidefinite
    minimize(trace(W*C))
    subject to
    % V constraints
    %W-diag(Vmin.^2) >= 0
    %W-diag(Vmax.^2) <= 0
    for k = 1:N
        trace(W*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
        trace(W*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
        trace(W*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
        trace(W*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
        W(k,k) >= Vmin(k).^2
        W(k,k) <= Vmax(k).^2
    end
    cvx_end
    if strcmp(cvx_status,'Solved')
        if rank(W) > 1
            rank(W)
            [s,v,d] = svd(W);
            v = diag(v);
            ratio = v(1)/sum(v);
            if ratio < 0.99
                ratio
                inexact = {inexact,file};
            end
        end
        OKfiles = {OKfiles,file};
    else
        notsolved = {notsolved,file};
    end
    cvx_clear
end



