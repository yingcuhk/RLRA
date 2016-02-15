function P_AD = minP(L,W,nuclearnorm)


    % min norm(P*L-W,'fro'),
    % s.t. g(Y,D) <= 0
    %Y = X;
    %Y = Y.*(Y>0);
    [M,~] = size(W);
    [K,~] = size(L);
    cvx_begin quiet
    variable P(M,K) 
    minimize(norm(P*L-W,'fro'))
    subject to
        norm(P*L,'fro') <= nuclearnorm
    cvx_end
    
    P_AD = P;
    cvx_clear
end