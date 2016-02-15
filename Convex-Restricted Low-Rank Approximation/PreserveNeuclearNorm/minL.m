function L_AD = minL(P,W,nuclearnorm)


    % min norm(P*L-W,'fro'),
    % s.t. g(Y,D) <= 0
    %Y = X;
    %Y = Y.*(Y>0);
    [~,N] = size(W);
    [~,K] = size(P);
    cvx_begin quiet
    variable L(K,N) 
    minimize(norm(P*L-W,'fro'))
    subject to
        norm(P*L,'fro') <= nuclearnorm
    cvx_end
    
    L_AD = L;
    cvx_clear
end