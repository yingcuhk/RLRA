  function [F]=opf(x)
%  x=[50 50 50 50 50;40 50 50 50 50];
%   x=[67.4162    2.6417    7.4733    2.1379   12.8046;
%       50 50 50 50 50;50 50 50 50 50];
[mm1 mm2]=size(x)
    for ii=1:mm1
    y1=x(ii,:);
basemva=100;
n=30;
accuracy=0.0001;
maxiter=5;
busdata=[1	1	1.06	0	0	0	48	0	0	50	0
2	2	1.043	0	21.7	12.7	40	0	-40	50	0
3	0	1	0	2.4	1.2	 0	0	0	0	0
4	0	1	0	7.6	1.6 0	0	0	30	0
5	2	1.01 0	94.2	19.0	0	0	 -40 40	0
6	0	1	0	0	0	0	0	0	30	0
7	0	1	0	22.8	10.9	0	0	0	0	0
8	2	1.01	0	30 30	0	0	-10	40	0
9	0	1	0	0	0	0	0	0	0	0
10	0	1	0	5.8	2	0	0	0	0	19
11	2	1.082	0	0	0	0	0	-6	24	0
12	0	1	0	11.2	7.5	0	0	0	0	0
13	2	1.071	0	0	0	0	0	-6	24	0
14	0	1	0	6.2	1.6	0	0	0	0	0
15	0	1	0	8.2	2.5	0	0	0	0	0
16	0	1	0	3.5	1.8	0	0	0	0	0
17	0	1	0	9	5.8	0	0	0	0	0
18	0	1	0	3.2	0.9	0	0	0	0	0
19	0	1	0	9.5	3.4	0	0	0	0	0
20	0	1	0	2.2	0.7	0	0	0	0	0
21	0	1	0	17.5	11.2	0	0	0	0	0
22	0	1	0	0	0	0	0	0	0	0
23	0	1	0	3.2	1.6	0	0	0	0	0
24	0	1	0	8.7	6.7	0	0	0	0	4.3
25	0	1	0	0	0	0	0	0	0	0
26	0	1	0	3.5	2.3	0	0	0	0	0
27	0	1	0	0	0	0	0	0	0	0
28	0	1	0	0	0	0	0	0	0	0
29	0	1	0	2.4	0.9	0	0	0	0	0
30	0	1	0	10.6	1.9	0	0	0	0	0
];     
linedata=[ 1	2	0.0192	0.0575	0.0264	1
1	3	0.0452	0.1852	0.0204	1
2	4	0.057	0.1737	0.0184	1
3	4	0.0132	0.0379	0.0042	1
2	5	0.0472	0.1983	0.0209	1
2	6	0.0581	0.1783	0.0187	1
4	6	0.0119	0.0414	0.0045	1
5	7	0.046	0.116	0.0102	1
6	7	0.0267	0.082	0.085	1
6	8	0.012	0.042	0.0045	1
6	9	0	0.208	0	.978
6	10	0	0.556	0	.969
9	11	0	0.208	0	1
9	10	0	0.11	0	1
4	12	0	0.256	0	.932
12	13	0	0.14	0	1
12	14	0.1231	0.2559	0	1
12	15	0.0662	0.1304	0	1
12	16	0.0945	0.1987	0	1
14	15	0.09	0.1997	0	1
16	17	0.0824	0.1932	0	1
15	18	0.107	0.2185	0	1
18	19	0.0639	0.1292	0	1
19	20	0.034	0.068	0	1
10	20	0.0936	0.209	0	1
10	17	0.0324	0.0845	0	1
10	21	0.0348	0.0749	0	1
28	22	0.0727	0.1499	0	1
21	22	0.0116	0.0236	0	1
15	23	0.1	0.202	0	1
22	24	0.115	0.179	0	1
23	24	0.132	0.27	0	1
24	25	0.1885	0.3292	0	1
25	26	0.2544	0.38	0	1
25	27	0.1093	0.2087	0	1
27	28	0	0.396	0	.968
27	29	0.2198	0.4153	0	1
27	30	0.3202	0.6027	0	1
29	30	0.2399	0.4533	0	1
8	28	0.0636	0.2	0.0214	1
6	28	0.0169	0.0599	0.065	1];
  
% formation of  Y bus
j=sqrt(-1); 
nl = linedata(:,1); nr = linedata(:,2); R = linedata(:,3);
X = linedata(:,4); Bc = j*linedata(:,5); a = linedata(:, 6);
nbr=length(linedata(:,1)); nbus = max(max(nl), max(nr));
Z = R + j*X; y= ones(nbr,1)./Z;        %branch admittance
for n = 1:nbr
if a(n) <= 0  a(n) = 1; else end
Ybus=zeros(nbus,nbus);     % initialize Ybus to zero
% formation of the off diagonal elements
for k=1:nbr;
       Ybus(nl(k),nr(k))=Ybus(nl(k),nr(k))-y(k)/a(k);
       Ybus(nr(k),nl(k))=Ybus(nl(k),nr(k));
    end
end
% formation of the diagonal elements
for  n=1:nbus
     for k=1:nbr
         if nl(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k)/(a(k)^2) + Bc(k);
         elseif nr(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k) +Bc(k);
         else, end
     end
end
gencost = [1 	0.02	2	0  0 80;
		2 0.0175	1.75	0  0 80;
		5 0.0625	1	    0  0 50;
		8 0.00834	3.25	0  0 55;
		11 0.025	3	    0  0 30;
		13 0.025	3	0  0 40];
    nn1=length(gencost(:,1));
    for i=2:nn1;
       xx=gencost(i,1);
       busdata(xx,7)=y1(i-1);
    end
    
ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));
for k=1:nbus
n=busdata(k,1);
kb(n)=busdata(k,2); Vm(n)=busdata(k,3); delta(n)=busdata(k, 4);
Pd(n)=busdata(k,5); Qd(n)=busdata(k,6); Pg(n)=busdata(k,7); Qg(n) = busdata(k,8);
Qmin(n)=busdata(k, 9); Qmax(n)=busdata(k, 10);
Qsh(n)=busdata(k, 11);
    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
end
for k=1:nbus
if kb(k) == 1, ns = ns+1; else, end
if kb(k) == 2 ng = ng+1; else, end
ngs(k) = ng;
nss(k) = ns;
end
Ym=abs(Ybus); t = angle(Ybus);
m=2*nbus-ng-2*ns;
maxerror = 1; converge=1;
iter = 0;
% Start of iterations
clear A  DC   J  DX
while maxerror >= accuracy & iter <= maxiter % Test for max. power mismatch
for i=1:m
for k=1:m
   A(i,k)=0;      %Initializing Jacobian matrix
end, end
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;
J11=0; J22=0; J33=0; J44=0;
   for i=1:nbr
     if nl(i) == n | nr(i) == n
        if nl(i) == n,  l = nr(i); end
        if nr(i) == n,  l = nl(i); end
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        if kb(n)~=1
        J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        else, end
        if kb(n) ~= 1  & kb(l) ~=1
        lk = nbus+l-ngs(l)-nss(l)-ns;
        ll = l -nss(l);
      % off diagonalelements of J1
        A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
              if kb(l) == 0  % off diagonal elements of J2
              A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));end
              if kb(n) == 0  % off diagonal elements of J3
              A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); end
              if kb(n) == 0 & kb(l) == 0  % off diagonal elements of  J4
              A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));end
        else end
     else , end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if kb(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P
     if kb(n) == 2  Q(n)=Qk;
         if Qmax(n) ~= 0
           Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
           if iter <= 7                  % Between the 2th & 6th iterations
              if iter > 2                % the Mvar of generator buses are
                if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
                Vm(n) = Vm(n) + 0.01;    % is changed in steps of 0.01 pu to
                elseif Qgc  > Qmax(n),   % bring the generator Mvar within
                Vm(n) = Vm(n) - 0.01;end % the specified limits.
              else, end
           else,end
         else,end
     end
   if kb(n) ~= 1
     A(nn,nn) = J11;  %diagonal elements of J1
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2
     A(lm,nn)= J33;        %diagonal elements of J3
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4
     DC(lm) = Q(n)-Qk;
   end
end
DX=A\DC';
for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 1
    delta(n) = delta(n)+DX(nn); end
    if kb(n) == 0
    Vm(n)=Vm(n)+DX(lm); end
 end
  maxerror=max(abs(DC));
%      if iter == maxiter & maxerror > accuracy 
%    fprintf('\nWARNING: Iterative solution did not converged after ')
%    fprintf('%g', iter), fprintf(' iterations.\n\n')
%    fprintf('Press Enter to terminate the iterations and print the results \n')
%    converge = 0; pause, else, end
%    
end

% if converge ~= 1
%    tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
%    tech=('                   Power Flow Solution by Newton-Raphson Method');
% end   
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
i=sqrt(-1);
k=0;
for n = 1:nbus
     if kb(n) == 1
     k=k+1;
     S(n)= P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);     %june 97
     elseif  kb(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  % June 1997
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';
Pgt = sum(Pg);  Qgt = sum(Qg); Pdt = sum(Pd); Qdt = sum(Qd); Qsht = sum(Qsh);
if Pgg(1)>gencost(1,6);
    Pgg(1)=gencost(1,6);
else
end
a1=gencost(:,2);b1=gencost(:,3);
c1=gencost(:,4);
 F1=(Pgg.*Pgg)*a1+Pgg*b1+sum(c1);
 F(ii,1)=F1;
 clear F1 Pgg
    end
    