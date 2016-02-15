function X = Xupdate(W,nuclearnorm)
    
    % min norm(Y-X,'fro'),
    % s.t. g(Y,D) <= 0
    %Y = X;
    %Y = Y.*(Y>0);
    [M,N] = size(W);
    
    cvx_begin quiet
    variable Y(M,N) 
    minimize(norm(Y-W,'fro'))
    subject to
        norm(Y,'fro') <= nuclearnorm
    cvx_end
    
    X = Y;
    
    cvx_clear
end