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
a = imread('chinese.png');
% A = [];
% for k =1:3
%     A = [A,double(a(:,:,k))/255];
% end
%a = imread('ie.jpg');
A = rgb2gray(a);
A = double(A)/255;
Edge = edge(A);
[M,N] = size(A);
Noise = randn(M,N)*1.5;

Noise(:,1:N*0.5) = 0;
W = A+Noise;
%imshow(A);
b = A*ones(N,1);
c = A'*ones(M,1);
mm = rand(M,N);
MM = zeros(M,N);
MM(:,1:N*0.5) = 1;
Ob = A.*MM;
%figure;imshow(W);
% a noise 
K = 5;  % upper bound of rank
TSVD = zeros(M,N);
[u,l,v] = svd(W);
l = diag(l);
for k = 1:K
    TSVD = TSVD+l(k)*u(:,k)*v(:,k)'; 
end
%figure;imshow(TSVD);
%pause(3);
close all;


rho = 10;
U = 0;
X = W;
Y = max(W,0);

maxiter = 50;  % max update 
Err = zeros(maxiter,1);
SNR = zeros(maxiter,1);
Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);
tic
for t = 1:maxiter
    
    % X update

    T = (W+rho/2*(Y-U))/(1+rho/2);
    [u,l,v] = svd(T);
    X = 0;
    l = diag(l);
    for k = 1:K
        X = X+l(k)*u(:,k)*v(:,k)'; 
    end
    
    % Y update
    %Y =  max(X+U,0);
    %Y = Yupdate(X+U,Edge);
    %Y = Yupdate(X+U,MM,Ob)
    Y = X+U;
    Y = Y.*(1-MM)+Ob;
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
% B = zeros(M,N,3);
% for k = 1:3
%    B(:,:,k) = X(:,(k-1)*N/3+1:k*N/3);
% end
figure;imshow(A)
figure;imshow(X)
[Px,SNRx] = psnr(X,A)
figure;imshow(TSVD)
[Pt,SNRt] = psnr(TSVD,A)
%figure;hold on; plot(Err,'r');
%figure; hold on; plot(Err,'r');plot(Err_AD,'LineWidth',2);plot(Err_MF);title('objective value')
%hold on;plot(Gap,'r');%plot(Gap_AD); title('norm(X-Y)')
% figure; hold on; plot(Resi,'r');plot(Resi_AD); title('Spectrum Residual')
% figure; hold on; plot(Err,'r');plot(Err_AD);title('objective value')