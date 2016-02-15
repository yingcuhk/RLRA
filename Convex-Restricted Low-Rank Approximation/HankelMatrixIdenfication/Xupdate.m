function X = Xupdate(T)
    
    % min norm(Y-X,'fro'),
    % s.t. g(Y) <= 0
    %Y = X;
    %Y = Y.*(Y>0);
    [M,N] = size(T);
    
    cvx_begin quiet
    variable W(M,M) semidefinite
    minimize(norm(W-T,'fro'))
    subject to

    cvx_end
    
    X = W;
    cvx_clear
end