function [Phi,Psi,M,T] = Matrix_SDP(Y)

    [N,~] = size(Y);
    Phi = [];
    Psi = [];
    M = [];
    T = [];
    for k = 1:N
        phi = zeros(N);
        psi = zeros(N);
        m = zeros(N);
        t = zeros(N);
%         for i = 1:N
%             for j = 1:N
%                 if k == i && k == i
%                     phi(i,j) = (Y(i,j)+conj(Y(i,j)))/2;
%                 elseif k == i
%                     phi(i,j) = 1/2*Y(i,j);
%                 elseif k == j
%                     phi(i,j) = 1/2*conj(Y(j,i));
%                 end
%             end
%         end
        e = zeros(N,1);
        e(k) = 1;
        Yk = e*e'*Y;
        phi = (Yk+Yk')/2;
%         if sum(abs(tphi-phi))>0
%             disp('WARNING, Phi');
%             pause;
%         end
        Phi = [Phi,phi];
%         for i = 1:N
%             for j = 1:N
%                 if k == i && k == i
%                     psi(i,j) = (Y(i,j)+conj(Y(i,j)))/2;
%                 elseif k == i
%                     psi(i,j) = -1/2i*Y(i,j);
%                 elseif k == j
%                     psi(i,j) = 1/2i*conj(Y(j,i));
%                 end
%             end
%         end
        psi = (Yk'-Yk)/2i;
        %psi
        %pause
%         if sum(abs(tpsi-psi))>0
%             disp('WARNING,Psi');
%             pause;
%         end
        Psi = [Psi,psi];
                 
    end
end