function [X,f_t,r_t] = RLRA4SDP(W,Phi,Psi,Pmax,Pmin,Qmax,Qmin,Vmax,Vmin,maxiter,rho,casename)

[M,N] = size(W);


K = 1;  % upper bound of rank


Y = 0;
[u,l,v] = svd(W);

% l = diag(l);
% for k = 1:K
%     Y = Y+l(k)*u(:,k)*v(:,k)';
% end
%rho = 200;
U = 0;
f_t = [];
r_t = [];
for t = 1:maxiter
    
    % X update
    X = Xupdate(Y,W,U,Phi,Psi,Pmax,Pmin,Qmax,Qmin,Vmax,Vmin,rho);
    
    
    % Y update
    T = X+U;
    [u,l,v] = svd(T);
    Y = 0;
    l = diag(l);
    for k = 1:K
        Y = Y+l(k)*u(:,k)*v(:,k)';
    end
    
    
    % U update
    U = U+X-Y;
    f_t = [f_t;norm(W-X,'fro')];
    r_t = [r_t;norm(X-Y,'fro')];
end


% title(casename)
% subplot(2,1,1);hold on;
% plot(f_t)
% title(['rho = ',num2str(rho),':distance, norm(X-W,''fro'')'])
% subplot(2,1,2);hold on
% plot(r_t)
% title('residual, norm(X-Y,''fro'')')
end