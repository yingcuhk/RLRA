function [d_im psnr_1 ssim_1 psnr_2 ssim_2]  = LPGPCA_color_denoising( nI, oI, v, profile, K )
s      =   2;
if strcmp(profile, 'fast')
    s   =  2;
elseif strcmp(profile, 'normal')
    s   =  2;
end
[n,m,ch]=size(oI);

figure(1),clf;
imshow(oI);

figure(2),clf;
imshow(nI);
[h w ch]    =   size(nI);

%%%%%%%%%%%%%% denoising %%%%%%%%%%%%%%%%%%%
%%%initial denoising
v2     =   v^2;
S      =   20;       %training block (2k+1)*(2k+1)
t      =   3;        %variable block (2t+1)*(2t+1)
nblk   =   250;
b      =   2*t+1;
b2     =   b*b;

k     =  0;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b*ch,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  nI(i:h-b+i,j:w-b+j, 1);
        blk  =  blk(:);
        X(k,:) =  blk';            
        blk  =  nI(i:h-b+i,j:w-b+j, 2);
        blk  =  blk(:);
        X(k+b2,:) =  blk';            
        blk  =  nI(i:h-b+i,j:w-b+j, 3);
        blk  =  blk(:);
        X(k+b2*2,:) =  blk';                            
    end
end
if strcmp(profile, 'fast')
    X1       =   dim_reduction(X);
    XT       =   X1';
else
    XT       =   X';
end
I        =   (1:L);
I        =   reshape(I, N, M);
N1       =   length(r);
M1       =   length(c);
L        =   N1*M1;
Y        =   zeros( b2*ch, L );

for  i  =  1 : N1
    for  j  =  1 : M1
        
        row      =   r(i);
        col      =   c(j);
        off      =   (col-1)*N + row;
        off1     =   (j-1)*N1 + i;        
        
        indc              =   LPG_new( XT, row, col, off, nblk, S, I );    
        [coe, P, V, mX]   =   getpca( X(:, indc) );
        py                =   mean(coe.^2, 2);
        px                =   max(0, py-v2);
        wei               =   px./py;
        Y(:,off1)         =   P'*(coe(:,1).*wei) + mX(:,1);
    
    end
end

% Output the processed image
dI       =  zeros(h,w,ch);
im_wei   =  zeros(h,w,ch);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        dI(r-1+i,c-1+j,1)       =  dI(r-1+i,c-1+j, 1) + reshape( Y(k,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,1)   =  im_wei(r-1+i,c-1+j,1) + 1;
        dI(r-1+i,c-1+j,2)       =  dI(r-1+i,c-1+j, 2) + reshape( Y(k+b2,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,2)   =  im_wei(r-1+i,c-1+j,2) + 1;
        dI(r-1+i,c-1+j,3)       =  dI(r-1+i,c-1+j, 3) + reshape( Y(k+b2*2,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,3)   =  im_wei(r-1+i,c-1+j,3) + 1;
    end
end
dI        =   dI./(im_wei+eps);
figure(3),clf;
imshow(dI(K+1:n-K,K+1:m-K,:));
imwrite(dI(K+1:n-K,K+1:m-K,:),'Results\basic_lpgpca.tif','tif');

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
%-----------------------------------------------------------------
% The second stage: refinement
%-----------------------------------------------------------------
nI    =  dI;
v     =  v1*0.37;
v2    =  v^2;

k     =   0;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b*ch,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  nI(i:h-b+i,j:w-b+j, 1);
        blk  =  blk(:);
        X(k,:) =  blk';            
        blk  =  nI(i:h-b+i,j:w-b+j, 2);
        blk  =  blk(:);
        X(k+b2,:) =  blk';            
        blk  =  nI(i:h-b+i,j:w-b+j, 3);
        blk  =  blk(:);
        X(k+b2*2,:) =  blk';       
    end
end
if strcmp(profile, 'fast')
    X1       =   dim_reduction(X);
    XT       =   X1';
else
    XT       =   X';
end
I        =   (1:L);
I        =   reshape(I, N, M);
N1       =   length(r);
M1       =   length(c);
L        =   N1*M1;
Y        =   zeros( b2*ch, L );

for  i  =  1 : N1
    for  j  =  1 : M1        
        row     =   r(i);
        col     =   c(j);
        off     =  (col-1)*N + row;
        off1    =  (j-1)*N1 + i;        
        
        indc              =   LPG_new( XT, row, col, off, nblk, S, I );    
        [coe, P, V, mX]   =   getpca( X(:, indc) );
        py                =   mean(coe.^2, 2);
        px                =   max(0, py-v2);
        wei               =   px./py;
        Y(:,off1)         =   P'*(coe(:,1).*wei) + mX(:,1);    
    end
end

% Output the processed image
dI       =  zeros(h,w,ch);
im_wei   =  zeros(h,w,ch);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        dI(r-1+i,c-1+j,1)       =  dI(r-1+i,c-1+j, 1) + reshape( Y(k,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,1)   =  im_wei(r-1+i,c-1+j,1) + 1;
        dI(r-1+i,c-1+j,2)       =  dI(r-1+i,c-1+j, 2) + reshape( Y(k+b2,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,2)   =  im_wei(r-1+i,c-1+j,2) + 1;
        dI(r-1+i,c-1+j,3)       =  dI(r-1+i,c-1+j, 3) + reshape( Y(k+b2*2,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j,3)   =  im_wei(r-1+i,c-1+j,3) + 1;
    end
end
dI        =   dI./(im_wei+eps);

psnr_2    =  csnr_color(oI,dI,K,K);
ssim_21   = cal_ssim( oI(:,:,1), dI(:,:,1), K, K );
ssim_22   = cal_ssim( oI(:,:,2), dI(:,:,2), K, K );
ssim_23   = cal_ssim( oI(:,:,3), dI(:,:,3), K, K );
ssim_2    = (ssim_21+ssim_22+ssim_23)/3;

figure(4),clf;
imshow(dI(K+1:n-K,K+1:m-K,:));

d_im  =  dI(K+1:n-K,K+1:m-K,:);
return;


function Y    =  dim_reduction(X)
n      =   size(X,1);
n      =   floor(n*0.4);
[coe, P, V, mX]   =   getpca( X );
Y      =   P(1:n,:)*X;
return;

