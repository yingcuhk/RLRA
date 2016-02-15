clc;clear;close all;
M = 20;N = 20;
W = randn(1,M+N-1);
W = hankel(W(1:M),W(M:end));
% by RLRA ADMM
rho = 5; % rho = 5 is good for exp 1,2
K = 5;
U = 0;
X = Xupdate(W);
[u,l,v] = svd(W);
Y = 0;
l = diag(l);
for k = 1:K
    Y = Y+l(k)*u(:,k)*v(:,k)';
end

P = ones(M,K);

Y_PL = Xupdate(W);
X_PL = zeros(size(Y_PL));

[P_AD,L_AD] = nnmf(Y_PL,K);
%P_AD = ones(M,K);
maxiter = 30;  % max update
Err = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);

Err_AD = zeros(maxiter,1);
Resi_AD = zeros(maxiter,1);
Gap_AD = zeros(maxiter,1);
%% ADMM

for t = 1:maxiter
    
    % X update

    t
    T = (W+rho/2*(Y-U))/(1+rho/2);
    
    X =  Xupdate(T);
    % Y update
    T = X+U;
    [u,l,v] = svd(T);
    Y = 0;
    l = diag(l);
    for k = 1:K
        Y = Y+l(k)*u(:,k)*v(:,k)';
    end
    % U update
    U = U+X-Y;
    
    %     [u,l,v] = svd(Y);
    %     l = diag(l);
    % statistics
    Err(t) = norm(W-X,'fro');
    %     Resi(t) = sum(l(1:K))/sum(l);
    Gap(t) = norm(X-Y,'fro');
    
    [u,l,v] = svd(Y_PL);
    X_PL = 0;
    l = diag(l);
    for k = 1:K
        X_PL = X_PL+l(k)*u(:,k)*v(:,k)';
    end
    
    Y_PL = Yupdate(X_PL);
    Err_PL(t) = norm(W-Y_PL,'fro');
    
    L_AD = minL(P_AD,W);
    P_AD = minP(L_AD,W);
    
    Err_AD(t) = norm(P_AD*L_AD-W,'fro');
    
    
end

figure; hold on;
plot(Err,'-sr',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','r')

plot(Err_PL,'--dk',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','k')

plot(Err_AD,':^g',...
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
