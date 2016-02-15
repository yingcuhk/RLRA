function [d_im psnr_1 ssim_1 psnr_2 ssim_2]  = LPGPCA( nI, oI, v )

% load data\paint\paint10.mat;


[n,m]=size(oI);

figure(1),clf;
imshow(blow(oI), [0 255]);

figure(2),clf;
imshow(nI, [0 255]);

%%%%%%%%%%%%%% denoising %%%%%%%%%%%%%%%%%%%
%%%initial denoising
k=20;       %training block (2k+1)*(2k+1)
b=3;        %variable block (2b+1)*(2b+1)
K=20;       %K>=k

dI=nI;
V=0;
Vk=0;

dI       =  zeros( n, m );
cnt_im   =  zeros( n, m );

for i=K+1:n-K
   for j=K+1:m-K
      Block=nI(i-k:i+k,j-k:j+k);
      CB=nI(i-b:i+b,j-b:j+b);
      [B,nr]=lpgpca_1(Block,CB,2*b+1,v,180,25);
      [dI cnt_im ]  =  Write2Im( B, b, dI, cnt_im, i, j, n, m );
      V=V+nr;
      Vk=Vk+1;
   end
end
dI(K-1:n-K+2,K-1:m-K+2)  =  dI(K-1:n-K+2,K-1:m-K+2)./cnt_im(K-1:n-K+2,K-1:m-K+2);


psnr_1    = csnr(oI,dI,K,K);
ssim_1   = cal_ssim( oI, dI, K, K );

figure(3),clf;
imshow(blow(dI(K+1:n-K,K+1:m-K)), [0 255]);


dif=dI(K+1:n-K,K+1:m-K)-nI(K+1:n-K,K+1:m-K);
vd=v^2-(mean(mean(dif.^2)));
v1=sqrt(abs(vd));

dI2       =  zeros( n, m );
cnt_im    =  zeros( n, m );


%-----------------------------------------------------------------
% The second stage: refinement
%
%-----------------------------------------------------------------

for i=K+1:n-K
   for j=K+1:m-K
       
      Block2  = dI(i-k:i+k,j-k:j+k);
      CB2     = dI(i-b:i+b,j-b:j+b);          
      [B,nr]  = lpgpca_1(Block2,CB2,2*b+1,0.35*v1, 150,17);
      [dI2 cnt_im ]  =  Write2Im( B, b, dI2, cnt_im, i, j, n, m );
   end
end

dI2(K-1:n-K+2,K-1:m-K+2)   =  dI2(K-1:n-K+2,K-1:m-K+2)./cnt_im(K-1:n-K+2,K-1:m-K+2);

psnr_2    =  csnr(oI,dI2,K,K);
ssim_2   =  cal_ssim( oI, dI2, K, K );

figure(4),clf;
imshow(blow(dI2(K+1:n-K,K+1:m-K)), [0 255]);

d_im  =  dI2(K+1:n-K,K+1:m-K);

