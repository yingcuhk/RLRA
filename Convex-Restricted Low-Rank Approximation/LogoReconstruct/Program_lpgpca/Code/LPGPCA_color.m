function [d_im psnr_1 ssim_1 psnr_2 ssim_2]  = LPGPCA_color( nI, oI, v )

[n m ch]=size(oI);


%%%%%%%%%%%%%% denoising %%%%%%%%%%%%%%%%%%%
%%%initial denoising
k=20;       %training block (2k+1)*(2k+1)
b=2;        %variable block (2b+1)*(2b+1)
K=20;       %K>=k

V=0;
Vk=0;

dI       =  zeros( n, m, ch );
cnt_im   =  zeros( n, m, ch );

for i=K+1:n-K
   for j=K+1:m-K
      Block  =  nI(i-k:i+k,j-k:j+k, :);
      CB     =  nI(i-b:i+b,j-b:j+b, :);
      [B(:,:,:),nr1] =  lpgpca_new(Block(:,:,:),CB(:,:,:),2*b+1,v,320,25);
      
      [dI(:,:,1) cnt_im(:,:,1) ]  =  Write2Im( B(:,:,1), b, dI(:,:,1), cnt_im(:,:,1), i, j, n, m );
      [dI(:,:,2) cnt_im(:,:,2) ]  =  Write2Im( B(:,:,2), b, dI(:,:,2), cnt_im(:,:,2), i, j, n, m );
      [dI(:,:,3) cnt_im(:,:,3) ]  =  Write2Im( B(:,:,3), b, dI(:,:,3), cnt_im(:,:,3), i, j, n, m );      

      V=V+nr1;
      Vk=Vk+1;
   end
end
dI(K+1-b:n-K+b,K+1-b:m-K+b,:)  =  dI(K+1-b:n-K+b,K+1-b:m-K+b,:)./cnt_im(K+1-b:n-K+b,K+1-b:m-K+b,:);

imwrite(dI(K+1:n-K,K+1:m-K,:)/255,'Results\basic_lpgpca.tif','tif');

psnr_1   = csnr_color(oI,dI,K,K);
ssim_1   = cal_ssim( oI(:,:,1), dI(:,:,1), K, K );
ssim_2   = cal_ssim( oI(:,:,2), dI(:,:,2), K, K );
ssim_3   = cal_ssim( oI(:,:,3), dI(:,:,3), K, K );
ssim_1   = (ssim_1+ssim_2+ssim_3)/3;



dif_1=dI(K+1:n-K,K+1:m-K,1)-nI(K+1:n-K,K+1:m-K,1);
dif_2=dI(K+1:n-K,K+1:m-K,2)-nI(K+1:n-K,K+1:m-K,2);
dif_3=dI(K+1:n-K,K+1:m-K,3)-nI(K+1:n-K,K+1:m-K,3);
vd_1=v^2-(mean(mean(dif_1.^2)));
vd_2=v^2-(mean(mean(dif_2.^2)));
vd_3=v^2-(mean(mean(dif_3.^2)));
vd  = (vd_1+vd_2+vd_3)/3;

v1=sqrt(abs(vd));



dI2       =  zeros( n, m, ch );
cnt_im    =  zeros( n, m, ch );
%-----------------------------------------------------------------
% The second stage: refinement
%
%-----------------------------------------------------------------

c  =  0.35;

for i=K+1:n-K
   for j=K+1:m-K
      Block2  = dI(i-k:i+k,j-k:j+k,:);
      CB2     = dI(i-b:i+b,j-b:j+b,:);      
      [B(:,:,:),nr] =  lpgpca_new(Block2(:,:,:),CB2(:,:,:),2*b+1,c*v1,300,17);
      
%       [B(:,:,1),nr] =  lpgpca_1(Block2(:,:,1),CB2(:,:,1),2*b+1,c*v1,150,17);
%       [B(:,:,2),nr] =  lpgpca_1(Block2(:,:,2),CB2(:,:,2),2*b+1,c*v1,150,17);      
%       [B(:,:,3),nr] =  lpgpca_1(Block2(:,:,3),CB2(:,:,3),2*b+1,c*v1,150,17);            


      [dI2(:,:,1) cnt_im(:,:,1)]  =  Write2Im( B(:,:,1), b, dI2(:,:,1), cnt_im(:,:,1), i, j, n, m );
      [dI2(:,:,2) cnt_im(:,:,2)]  =  Write2Im( B(:,:,2), b, dI2(:,:,2), cnt_im(:,:,2), i, j, n, m );
      [dI2(:,:,3) cnt_im(:,:,3)]  =  Write2Im( B(:,:,3), b, dI2(:,:,3), cnt_im(:,:,3), i, j, n, m );      
      
   end
end

dI2(K+1-b:n-K+b,K+1-b:m-K+b,:)   =  dI2(K+1-b:n-K+b,K+1-b:m-K+b,:)./cnt_im(K+1-b:n-K+b,K+1-b:m-K+b,:);

psnr_2    =  csnr_color(oI,dI2,K,K);
ssim_21   = cal_ssim( oI(:,:,1), dI2(:,:,1), K, K );
ssim_22   = cal_ssim( oI(:,:,2), dI2(:,:,2), K, K );
ssim_23   = cal_ssim( oI(:,:,3), dI2(:,:,3), K, K );
ssim_2   = (ssim_21+ssim_22+ssim_23)/3;

d_im  =  dI2(K+1:n-K,K+1:m-K,:);

