function Lie = LieFind(X1, dimsF, dimsX1, Lie, i, TrX1)
    
    % dimsX1 => n_v + 1 (Eq. 7)
    % dimsX2 => n_w + 1 (Eq. 7)
    dimsX2 = dimsX1 + dimsF - 1; % (n_v + 1) + (n_f + 1) -1
                                 % in general dims are vectors eg. [3,3]+[3,3]-1=[5,5]
   
                                 %see solveSDPAlternative.m for reasoning
                                 %of the following lines
    d = length(dimsX2);                                 
    bb = max(dimsX2); 
    B = cell(1,d);
    [B{:}] = ndgrid(-bb:bb);
    B = cellfun(@(M) M(:), B, 'uniform', 0);
    S = [B{end:-1:1}]';

        for l = 1 : length(dimsF) % for each dimension of F, theta1, ...

            for j = S(:,:) % for each coefficient
                
                % Accumulates the lie derivative (RHS of equation (17) in "https://www.preprints.org/manuscript/202409.1361/v1")
                Lie = Lie - TrFind(X1, dimsX1, TrX1, j) * F(l,i - j) * (sqrt(-1)) * j(l) + TrFind(X1, dimsX1, TrX1, (i - j)) * (sqrt(-1)) * j(l) * F(l,j);
            
            end
            
        end
  
    

    

    
end