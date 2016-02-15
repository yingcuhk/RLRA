%%
clc;
close all;
clear all;

% min(X-D), s.t. 
%% input and initialize
%W = randn(20,30);
%M = 100;
%N = 80;
%W = randn(M,N);
%load data.mat
%[M,N] = size(W);

%a = imread('mit.jpg');
%a = imread('chinese.png');
% A = [];
% for k =1:3
%     A = [A,double(a(:,:,k))/255];
% end
a = double(imread('ie.jpg'))/255;
A = rgb2gray(a);
%A = [];
% for k =1:3
%     A = [A,double(a(:,:,k))];
% end
%A = a;
%A = double(A)/255;
%Edge = edge(A);
[M,N] = size(A);
W = A+randn(M,N)*2;
%imshow(A);
b = A*ones(N,1);
c = A'*ones(M,1);
mm = rand(M,N);
MM = mm > 0.95;
Ob = A.*MM;
%figure;imshow(W);
% a noise 
K = 20;  % upper bound of rank
TSVD = zeros(M,N);

W_Ob = W.*(1-MM)+Ob;
[u,l,v] = svd(W_Ob);
l = diag(l);
for k = 1:K
    TSVD = TSVD+l(k)*u(:,k)*v(:,k)'; 
end
%figure;imshow(TSVD);
%pause(3);
close all;


rho = 15;
U = 0;
X = W_Ob;
Y = max(W,0);

maxiter = 30;  % max update 
Err = zeros(maxiter,1);
SNR = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);
tic
for t = 1:maxiter
    
    T = (W_Ob+rho/2*(Y-U))/(1+rho/2);
    X = T.*(1-MM)+Ob;
    
    T = X+U;
    
    [u,l,v] = svd(T);
    Y = 0;
    l = diag(l);
    for k = 1:K
        Y = Y+l(k)*u(:,k)*v(:,k)'; 
    end
    
    % U update
    U = U+X-Y;
    
    [u,l,v] = svd(Y);
    l = diag(l);
    % statistics
    Err(t) = norm(A-Y,'fro');
    Resi(t) = sum(l(1:K))/sum(l);
    Gap(t) = norm(X-Y,'fro');
    [SNR(t),~] = psnr(X,A);
%     
end
toc

%figure;plot(SNR)
% B = zeros(M,N/3,3);
% for k = 1:3
%    B(:,:,k) = W(:,(k-1)*N/3+1:k*N/3);
% end

figure;imshow(W)
% for k = 1:3
%    B(:,:,k) = TSVD(:,(k-1)*N/3+1:k*N/3);
% end
% [Px,SNRx] = psnr(B,a)
figure;imshow(TSVD)
% for k = 1:3
%    B(:,:,k) = X(:,(k-1)*N/3+1:k*N/3);
% end
% 
% [Px,SNRx] = psnr(B,a)
figure;imshow(X)

%[Px,SNRx] = psnr(X,A)
%figure;imshow(TSVD)
%[Pt,SNRt] = psnr(TSVD,A)
%figure;hold on; plot(Err,'r');
%figure; hold on; plot(Err,'r');plot(Err_AD,'LineWidth',2);plot(Err_MF);title('objective value')
%hold on;plot(Gap,'r');%plot(Gap_AD); title('norm(X-Y)')
% figure; hold on; plot(Resi,'r');plot(Resi_AD); title('Spectrum Residual')
% figure; hold on; plot(Err,'r');plot(Err_AD);title('objective value')