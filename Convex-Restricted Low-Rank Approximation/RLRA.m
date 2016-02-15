%%
clc;close all;clear all;

%% input and initialize
W = randn(20,30);
W = W*W';



K = 10;  % upper bound of rank

rho = 20; % rho = 5 is good for exp 1,2
U = 0;
X = W;
Y = W;
X_AD = W;
Y_AD = W;
maxiter = 25;  % max update 
Err = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);

Err_AD = zeros(maxiter,1);
Resi_AD = zeros(maxiter,1);
Gap_AD = zeros(maxiter,1);
%% ADMM

for t = 1:maxiter
    
    % X update

    T = (W+rho/2*(Y-U)/(1+rho/2));
    [u,l,v] = svd(T);
    X = 0;
    l = diag(l);
    for k = 1:K
        X = X+l(k)*u(:,k)*v(:,k)'; 
    end
    
    % Y update
    Y =  Yupdate(X+U,W);
    
    
    % U update
    U = U+X-Y;
    
    [u,l,v] = svd(Y);
    l = diag(l);
    % statistics
    Err(t) = norm(W-Y,'fro');
    Resi(t) = sum(l(1:K))/sum(l);
    Gap(t) = norm(X-Y,'fro');
    
    T_AD = Y_AD;
    [u,l,v] = svd(T_AD);
    X_AD = 0;
    l = diag(l);
    for k = 1:K
        X_AD = X_AD+l(k)*u(:,k)*v(:,k)'; 
    end
       
    
    Y_AD = Yupdate(X_AD,W);
  
    [u,l,v] = svd(Y_AD);
    l = diag(l);
    % statistics
    Err_AD(t) = norm(Y_AD-W,'fro');
    Resi_AD(t) = sum(l(1:K))/sum(l);
    Gap_AD(t) = norm(X_AD-Y_AD,'fro');
end
figure; hold on; plot(Gap,'r');plot(Gap_AD); title('norm(X-Y)')
figure; hold on; plot(Resi,'r');plot(Resi_AD); title('Spectrum Residual')
figure; hold on; plot(Err,'r');plot(Err_AD);title('objective value')