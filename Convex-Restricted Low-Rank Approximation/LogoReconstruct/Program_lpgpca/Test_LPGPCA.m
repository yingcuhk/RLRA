%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% Please cite the paper when you use th code:
%
% Lei Zhang, Weisheng Dong, David Zhang, Guangming Shi 
% Two-stage Image Denoising by Principal Component Analysis with Local
% Pixel Grouping, Pattern Recognition, vol. 43, issue 4, pp. 1531-1549,
% Apr. 2010
%
% AUTHORS:
%     Weisheng Dong, Lei Zhang, email:
%     wsdong@mail.xidian.edu.cn / lzhang@comp.polyu.edu.hk
%
%  Latest modify at Mar. 17, 2011
%
%  Two profiles for this code: 'fast' and  'normal'
%  'fast'---fastest speed but slight performance decreases;
%  'normal'---higher performance but computation time increases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
addpath('Code');

profile  =   'fast';   
v        =   20;
%oI       =   double( imread('Images\house.tif') );

oI       =   double( imread('ie.jpg') );

seed   =  0;
randn( 'state', seed );
%noise      =   randn(size( oI ));
%noise      =   noise/sqrt(mean2(noise.^2));
%nI         =   oI + v*noise;
a = imread('ie.jpg');
A = rgb2gray(a);
A = double(A)/255;
Edge = edge(A);
[M,N] = size(A);
W = A+randn(M,N)*2;
K          =   0;    % The width of the excluded boundaries, set to 20 to get the results in our paper

[ d_im psnr1 ssim1 psnr2 ssim2 ]   =  LPGPCA_denoising( W,A, v, profile, K );
imwrite( d_im/255,'lena_lpgpca.tif','tif' );
psnr1
ssim1
psnr2
ssim2

