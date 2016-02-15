function X = ConvexProjecton(W,Y,U,rho,Phi,Psi,Vmin,Vmax)

P = (W+rho*(Y-U)/(2+rho)); % the original point

% the purpose of this function is trying to projecting the point P onto a
% convex set defined by the constraints of Optimal Power Flow



%% SDP formulation
disp('Matrix preparation ready')
%pause
cvx_begin sdp
    variable W(N,N) complex semidefinite
    minimize(nrom(W-P),'fro')
    subject to
    % V constraints
        %W-diag(Vmin.^2) >= 0
        %W-diag(Vmax.^2) <= 0
        %W == complex semidefinite(N)
        for k = 1:N
             trace(W*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
             trace(W*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
             trace(W*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
             trace(W*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
             W(k,k) >= Vmin(k).^2
             W(k,k) <= Vmax(k).^2
        end
cvx_end

X = W;
emd