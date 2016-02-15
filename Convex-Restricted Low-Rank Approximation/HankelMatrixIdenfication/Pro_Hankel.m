function H_bar = Pro_Hankel(H)
    
    [M,N] = size(H);
    
    for k = 2:M+N
        total = 0;
        count = 0;
        for m = 1:M
            n = k-m;
            if n <= 0 
                continue;
            end
            if n <= N
                total = total + H(m,n);
                count = count + 1;
            end
        end
        ave = total/count;
        for m = 1:M
            n = k-m;
            if n <= 0 
                continue;
            end
            if n <= N
                H(m,n) = ave;
            end
        end
        
    end
    H_bar = H;
end