function [I,nr]=lpgpca_new(I,CB,w,v,num,th)

%%pca-based denoising, w*w vector block
%%I:   input noisy block n*m

[n,m,ch]=size(I);
N=n-w+1;
M=m-w+1;
L=N*M;

X=zeros(w*w*ch,L);

CB2(:,:,1)  =  CB(:,:,1)';
CB2(:,:,2)  =  CB(:,:,2)';
CB2(:,:,3)  =  CB(:,:,3)';
CB          =  CB2;
w2          =  w*w;
%%%
% CB=CB';
CB=CB(:);
CB=repmat(CB,1,L);

k=0;
for i=1:w
   for j=1:w
      k=k+1;
      T=I(i:n-w+i,j:m-w+j,:);
      T1  =  T(:,:,1);
      T2  =  T(:,:,2);
      T3  =  T(:,:,3);
      
      T1  =  T1(:);
      X(k,:)  = T1';

      T2  =  T2(:);
      X(k+w2,:)  = T2';

      T3  =  T3(:);
      X(k+w2*2,:)  = T3';
      
%       T=T(:);
%       X(k,:)=T';
   end
end

%%%grouping
E=abs(X-CB).^2;
mE=mean(E);
[val,ind]=sort(mE);

%%%%%%%%
ind1=find(val<th+2*v^2);
n1=length(ind1);
if n1>500
   X=X(:,ind(1:500));
elseif n1<num
   X=X(:,ind(1:num));
else
   X=X(:,ind(1:n1));
end

v2=v^2;
w2=w^2*ch;

% [X J]  =  scale_data( X, 0.75 );
[Y, P, V, mX]=getpca(X);
nr=0;


for i=1:w2
   y=Y(i,:);
   py=mean(y.^2)+0.01;
   mpv=max(0,py-v2);
   c=mpv/py;
   Y(i,:)=c*Y(i,:);
   nr=nr+mpv*v2/(mpv+v2);
end

% num  =  nn;
% for i=1:w2
%    y=Y(i,1:num);
%    mm  =  mean(y,2);
%    py=mean((y-mm).^2)+0.0001;
%    mpv=max(0,py-v2);
%    c=mpv/py;
%    Y(i,1)=c*(Y(i,1)-mm)+mm;
%    nr=nr+mpv*v2/(mpv+v2);
% end


nr=nr/w2;

B=(P'*Y+mX);

I=B(:,1);
I=reshape(I,w,w,ch);
% I=I';
I1(:,:,1)  = I(:,:,1)';
I1(:,:,2)  = I(:,:,2)';
I1(:,:,3)  = I(:,:,3)';
I   =  I1;

return;


function [X_out J]  =  scale_data( X, r )

[h w]  =  size( X );
s      =  h/3;

m1     =  repmat( mean( X(1:s, :), 1 ), s, 1);
m2     =  repmat( mean( X(s+1:s*2, :), 1 ), s, 1);
m3     =  repmat( mean( X(2*s+1:3*s, :), 1 ), s, 1);
J      =  [m1;m2;m3].*r;

X_out  =  X + J;