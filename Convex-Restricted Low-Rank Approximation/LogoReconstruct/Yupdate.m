function Y = Yupdate(X,Edge)
    
    % min norm(Y-X,'fro'),
    % s.t. g(Y,D) <= 0
    %W*1=b
    [M,N] = size(X);
    
    cvx_begin sdp quiet
    variable W(M,N)
    minimize(norm(W-X,'fro'))
    subject to
        edge(X) == Edge
        %W.*MM == Ob
        %W*ones(N,1) == b
        %W'*ones(M,1) == c
        %norm(W,'fro') <= norm(D,'fro')/M/N
        %diag(W) == diag(D) 
        %trace(W) <= trace(D)
        %min(min(W))>=0
    % W constraints
        % g(X)<=0
        %W >= 0 
%         for k = 1:N
%              trace(W*Phi(:,(k-1)*N+1:k*N)) <= Pmax(k)
%              trace(W*Phi(:,(k-1)*N+1:k*N)) >= Pmin(k)
%              trace(W*Psi(:,(k-1)*N+1:k*N)) <= Qmax(k)
%              trace(W*Psi(:,(k-1)*N+1:k*N)) >= Qmin(k)
%              W(k,k) >= Vmin(k).^2
%              W(k,k) <= Vmax(k).^2
%         end
    cvx_end
    
    Y = W;
    cvx_clear
end