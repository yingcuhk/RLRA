function H = MosiacHankelMatrix(P,Mv,Nv)
    
   
    H = zeros(sum(Mv),sum(Nv));
    marker = 1;
    for m = 1:length(Mv)
        for n = 1:length(Nv)
            
            M = Mv(m);
            N = Nv(n);
            vec = P(marker:marker+M+N-1-1);
            %whos vec
            marker = marker+M+N-1;
            h = hankel(vec(1:M),vec(M:end));
            m_pos = sum(Mv(1:m))-Mv(m)+1;
            n_pos =sum(Nv(1:n))-Nv(n)+1;
            %m_pos,n_pos
            H(m_pos:m_pos+M-1,n_pos:n_pos+N-1) = h;
        end
    end
    
    


end