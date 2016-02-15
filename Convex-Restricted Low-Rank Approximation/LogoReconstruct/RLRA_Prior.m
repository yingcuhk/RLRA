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


%a = imread('chinese.png');
% A = [];
% for k =1:3
%     A = [A,double(a(:,:,k))/255];
% end
%a = imread('ie.jpg');
a = imread('mit.jpg');
A = rgb2gray(a);
A = double(A)/255;
Edge = edge(A);
[M,N] = size(A);
W = A+randn(M,N)*2;
%imshow(A);
b = A*ones(N,1);
c = A'*ones(M,1);
mm = rand(M,N);

SNR = [];
thre = 0.9;
for thre = 0.95:-0.05:0.5
MM = mm > thre;
Ob = A.*MM;

W_Ob = W.*(1-MM)+Ob;
%figure;imshow(W);
% a noise 
K = 5;  % upper bound of rank
TSVD = zeros(M,N);
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
X = W;
X = X.*(1-MM)+Ob;
Y = TSVD;

maxiter = 50;  % max update 
Err = zeros(maxiter,1);

Resi = zeros(maxiter,1);
Gap = zeros(maxiter,1);
k = 1;
tic
for t = 1:maxiter
    
    T = (W+rho/2*(Y-U))/(1+rho/2);
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
    
    [u,l,v] = svd(X);
    l = diag(l);
    % statistics
    Err(t) = norm(A-Y,'fro');
    Resi(t) = sum(l(1:K))/sum(l);
    Gap(t) = norm(X-Y,'fro');
    %[SNR(t),~] = psnr(X,A);
%     
end
SNR = [SNR, psnr(X,A)];

end
toc
X = 0.95:-0.05:0.5;
X = 1-X;
plot(X,SNR,'-sr',...
    'LineWidth',3,...
    'MarkerSize',8,...
    'MarkerEdgeColor','r')

grid on
box on
xlabel('Percentage of unpolluted pixels')
ylabel('PSNR')
set(gca,'FontSize',18)
%figure;plot(SNR)
% B = zeros(M,N,3);
% for k = 1:3
%    B(:,:,k) = X(:,(k-1)*N/3+1:k*N/3);
% end
% figure;imshow(A)
% figure;imshow(Y)
% [Px,SNRx] = psnr(Y,A)
% figure;imshow(TSVD)
% [Pt,SNRt] = psnr(TSVD,A)
%figure;hold on; plot(Err,'r');
%figure; hold on; plot(Err,'r');plot(Err_AD,'LineWidth',2);plot(Err_MF);title('objective value')
%hold on;plot(Gap,'r');%plot(Gap_AD); title('norm(X-Y)')
% figure; hold on; plot(Resi,'r');plot(Resi_AD); title('Spectrum Residual')
% figure; hold on; plot(Err,'r');plot(Err_AD);title('objective value')