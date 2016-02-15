%%
clc;
close all;
clear;

% min(X-D), s.t. 
%% input and initialize
%W = randn(20,30);
% M = 100;
% N = 80;
% W = randn(M,N);
load data.mat
[M,N] = size(W);
W = W;
K = 20;  % upper bound of rank

rho = 20; 
U = 0;
X = max(W,0);
Y = W;
%Y = max(W,0);
X_AD = W;
Y_AD = max(W,0);
maxiter = 30;  % max update 
Err = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);
%Dual = zeros(maxiter,1);
Err_AP = zeros(maxiter,1);
%% Alternating Direction Optimization
% P_AD = rand(M,K);
% L_AD = (X\P_AD)';
[P_AD,L_AD] = nnmf(X,K);
Err_PL = zeros(maxiter,1);
%% NMF

Lambda = zeros(M,K);
Pi = zeros(K,N);
alpha = 1;
beta= 1;
gamma = 0.01;
for t = 1:maxiter
    % Y update
    %Y =  max(X+U,0);
    % X update
    
    T = (W+rho/2*(Y-U))/(1+rho/2);
    X = max(T,0);
    
    T = X+U;
    
    [u,l,v] = svd(T);
    Y = 0;
    l = diag(l);
    for k = 1:K
        Y = Y+l(k)*u(:,k)*v(:,k)'; 
    end
    
    
    % U update
    U = U+X-Y;
    
    %[u,l,v] = svd(Y);
    %l = diag(l);
    % statistics
    %Dual(t) = norm(X-Y,'from');
    Err(t) = norm(W-Y,'fro');
    Resi(t) = sum(l(1:K))/sum(l);
    Gap(t) = norm(X-Y,'fro');
%     
    %% by Lift and Projection, X_AD,Y_AD
    [u,l,v] = svd(Y_AD); % projection
    X_AD = 0;
    l = diag(l);
    for k = 1:K
        X_AD = X_AD+l(k)*u(:,k)*v(:,k)'; 
    end
    
    Y_AD = max(X_AD,0); % lift
    Err_AD(t) = norm(W-Y_AD,'fro'); 
    %% 
    P_AD = minP(L_AD,W);
    L_AD = minL(P_AD,W);
    Err_PL(t) = norm(W-P_AD*L_AD,'fro');

end
opt = statset('maxiter',maxiter,'display','iter');
[A,B] = nnmf(W,K,'options',opt);
matlab_obj = norm(A*B-W,'fro');
% hold on; plot(Err,'r');
figure; hold on; 
plot(Err,'-sr',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','r')

plot(Err_AD,'--dk',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','k')

plot(Err_PL,':^g',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','g')

% plot(Err_MF,':^b',...
%     'LineWidth',2,...
%     'MarkerSize',8,...
%     'MarkerEdgeColor','b')
set(gca,'FontSize',18)
xlim([0 maxiter])
grid on
box on
xlabel('Number of iterations')
ylabel('Objective values')
h = legend('ADMM-RLRA','Lift-and-Project','AD-Optimization');
set(h, 'Box', 'On', 'FontSize', 18, 'Orientation', 'vertical', 'Location', 'NorthEast','FontName', 'Arial');
%plot(Err_AD,'blue','LineWidth',2);
%plot(Err_MF,'black','LineWidth',2);
%hold on;plot(Gap,'r');%plot(Gap_AD); title('norm(X-Y)')
% figure; hold on; plot(Resi,'r');plot(Resi_AD); title('Spectrum Residual')
% figure; hold on; plot(Err,'r');plot(Err_AD);title('objective value')