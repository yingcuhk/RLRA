%%
clc;
%close all;
clear all;

% min(X-D), s.t.
%% input and initialize
%W = randn(20,30);
%M = 100;
%N = 80;
%W = randn(M,N);
load data.mat
[M,N] = size(W);

% a = imread('kaola.jpg');
% A = rgb2gray(a);
% A = double(A)/255;
% W = A;
% [M,N] = size(W);

K = 20;  % upper bound of rank

rho = 30; % rho = 5 is good for exp 1,2
U = 0;
X = max(W,0);
Y = W;

maxiter = 50;  % max update
Err = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);
Dual = zeros(maxiter,1);

Err_AD = zeros(maxiter,1);
X_diff = zeros(maxiter,1);
Err_MF = zeros(maxiter,1);
Resi_MF = zeros(maxiter,1);
Gap_MF = zeros(maxiter,1);
%% ADMM

%% NMF
X_MF = zeros(M,K);
Y_MF = rand(K,N);
U_MF = zeros(M,K);
V_MF = zeros(K,N);

Lambda = zeros(M,K);
Pi = zeros(K,N);
alpha = 1;
beta= 1;
gamma = 0.01;
Rho = [1,20,70,100];
for k = 1:4
    rho = Rho(k);
    U = 0;
    X = max(W,0);
    Y = W;
    for t = 1:maxiter
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
        
    end
    % hold on; plot(Err,'r');
    hold on;
    plot(Gap,...
        'LineWidth',3,...
        'MarkerSize',8)
end
set(gca,'FontSize',18)
legend('\rho=5','\rho=20','\rho=35','\rho=50')
%xlim([0 maxiter])
%grid on
%box on
xlabel('Number of iterations')
ylabel('Objective values')