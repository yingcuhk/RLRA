function [d_im psnr_1 ssim_1 psnr_2 ssim_2]  = LPGPCA_denoising( nI, oI, v, profile, K )
s      =   2;
if strcmp(profile, 'fast')
    s   =  2;
elseif strcmp(profile, 'normal')
    s   =  1;
end
[n,m]=size(oI);

figure(1),clf;
imshow(oI);

figure(2),clf;
imshow(nI);
[h w]      =   size(nI);

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
X     =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  nI(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
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
Y        =   zeros( b2, L );

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
dI       =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        dI(r-1+i,c-1+j)      =  dI(r-1+i,c-1+j) + reshape( Y(k,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;
    end
end
dI        =   dI./(im_wei+eps);
imwrite(dI(K+1:n-K,K+1:m-K,:),'Results\basic_lpgpca.tif','tif');

psnr_1    =   csnr(oI,dI,K,K);
ssim_1    =   cal_ssim( oI, dI, K, K );
figure(3),clf;
imshow((dI(K+1:n-K,K+1:m-K)));

dif=dI(K+1:n-K,K+1:m-K)-nI(K+1:n-K,K+1:m-K);
vd=v^2-(mean(mean(dif.^2)));
v1=sqrt(abs(vd));

%-----------------------------------------------------------------
% The second stage: refinement
%-----------------------------------------------------------------
nI    =  dI;
v     =  v1*0.36;
v2    =  v^2;

k     =   0;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  nI(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
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
Y        =   zeros( b2, L );

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
dI       =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        dI(r-1+i,c-1+j)      =  dI(r-1+i,c-1+j) + reshape( Y(k,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
    end
end
dI        =   dI./(im_wei+eps);


psnr_2    =  csnr(oI,dI,K,K);
ssim_2    =  cal_ssim( oI, dI, K, K );

figure(4),clf;
imshow(dI(K+1:n-K,K+1:m-K));

d_im  =  dI(K+1:n-K,K+1:m-K);
return;

function Y    =  dim_reduction(X)
n      =   size(X,1);
n      =   floor(n*0.4);
[coe, P, V, mX]   =   getpca( X );
Y      =   P(1:n,:)*X;
return;

