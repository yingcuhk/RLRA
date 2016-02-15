function New_X = Xupdate(Y,W,U,Phi,Psi,Pmax,Pmin,Qmax,Qmin,Vmax,Vmin,rho)


D = (W+rho/2*(Y-U))/(1+rho/2);
[N,temp] = size(W);
cvx_begin sdp quiet
variable X(N,N) complex semidefinite
minimize(norm(X-D,'fro'))
subject to
% V constraints
%W-diag(Vmin.^2) >= 0
%W-diag(Vmax.^2) <= 0
for k = 1:N
    trace(X*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
    trace(X*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
    trace(X*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
    trace(X*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
    X(k,k) >= Vmin(k).^2
    X(k,k) <= Vmax(k).^2
end
cvx_end
cvx_clear
New_X = X;
end

% 
% for k = 1:N
%     trace(X*Phi(:,(k-1)*N+1:k*N)) - Pmax(k)
%     -trace(X*Phi(:,(k-1)*N+1:k*N)) + Pmin(k)
%     trace(X*Psi(:,(k-1)*N+1:k*N)) - Qmax(k)
%     -trace(X*Psi(:,(k-1)*N+1:k*N)) + Qmin(k)
%     -X(k,k) + Vmin(k).^2
%     X(k,k) - Vmax(k).^2
% end