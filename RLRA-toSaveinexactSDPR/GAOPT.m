


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
    return;
end

try
    Ybus = makeYbus(opf);
catch ME
    disp('Forming Y bus failed')
    return
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
options = gaoptimset;

options = gaoptimset('PopulationSize', 50,'Display','iter','Generations', 200,'StallGenLimit',100,'TimeLimit', 300,'StallTimeLimit', 300,'PlotFcns',  {@gaplotbestf,@gaplotbestindiv});


[x ff]=ga(@opf1,5,options);