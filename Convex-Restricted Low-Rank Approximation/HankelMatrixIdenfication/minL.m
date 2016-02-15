function L_AD = minL(P,W)


    % min norm(P*L-W,'fro'),
    % s.t. g(Y,D) <= 0
    %Y = X;
    %Y = Y.*(Y>0);
    [M,N] = size(W);
    [~,K] = size(P);
    cvx_begin quiet
    variable L(K,N)  
    minimize(norm(P*L-W,'fro'))
    subject to
       P*L <In> semidefinite(M)
    cvx_end
    cvx_clear
    
    L_AD = L;
end