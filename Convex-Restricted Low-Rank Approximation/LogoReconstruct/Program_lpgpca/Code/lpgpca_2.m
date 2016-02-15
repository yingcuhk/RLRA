function [I]=lpgpca_2( I, CB, v_n, I2, CB2, w, v, num, th )

%%pca-based denoising, w*w vector block
%%I:   input noisy block n*m

[n,m]=size(I);
N=n-w+1;
M=m-w+1;
L=N*M;

X   = zeros(w*w,L);
X2  = zeros(w*w,L);

%%%
CB   = CB';
CB   = CB(:);
CB   = repmat(CB,1,L);

CB2  = CB2';
CB2  = CB2(:);
CB2  = repmat(CB2,1,L);

k=0;
for i=1:w
   for j=1:w
      k  = k+1;
      T  = I(i:n-w+i,j:m-w+j);
      T  = T(:);
      X(k,:)=T';

      T2 = I2(i:n-w+i,j:m-w+j);
      T2 = T2(:);
      X2(k,:)=T2';

   end
end

%%%grouping
E=abs(X2-CB2).^2;
mE=mean(E);
[val,ind]=sort(mE);

%%%%%%%%
ind1=find(val<th+2*v^2);
n1=length(ind1);
if n1>500
   X2 = X2(:,ind(1:500));
   X  = X(:,ind(1:500));   
elseif n1<num
   X2 = X2(:,ind(1:num));
   X  = X(:,ind(1:num));   
else
   X2 = X2(:,ind(1:n1));
   X  = X(:,ind(1:n1));   
end


[Y, P, V, mX]=getpca(X);
nr=0;

% X2  = X2 - mX;
% Y2  = P*X2;

v2 = v^2;
vn = v_n^2;
w2=w^2;

for i=1:w2
   y     = Y(i,:);   
   py    = mean(y.^2)+0.01;   % compute signal var   
   mpv   = max(0,py-vn);   
   c     = mpv/py;
   Y(i,:)= c*Y(i,:);
end

B=(P'*Y+mX);

I=B(:,1);
I=reshape(I,w,w);
I=I';
return;