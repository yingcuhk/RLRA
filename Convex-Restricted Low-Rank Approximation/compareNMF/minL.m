function L = minL(P,W)
    % min |PL-W|, PL>=0
    [M,N] = size(W);
    [~,K] = size(P);
    cvx_begin 
    variable L(K,N) 
    minimize(norm(P*L-W,'fro'))
    subject to
        P*L >= 0
%         for k = 1:N
%              trace(W*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
%              trace(W*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
%              trace(W*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
%              trace(W*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
%              W(k,k) >= Vmin(k).^2
%              W(k,k) <= Vmax(k).^2
%         end
    cvx_end
    
    cvx_clear
end