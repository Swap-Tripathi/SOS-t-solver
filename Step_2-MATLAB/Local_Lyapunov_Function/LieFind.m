function Lie = LieFind(X1, dimsF, dimsX1, i)
    
    

    dimsX2 = dimsX1 + dimsF - 1; 

    d = length(dimsX2); 
    bb = max(dimsX2); 
    B = cell(1,d);
    [B{:}] = ndgrid(-bb+1:bb-1);
    B = cellfun(@(M) M(:), B, 'uniform', 0);
    S = [B{end:-1:1}]';
    Lie = 0;

    preDefineJ = sqrt(-1);
    
        for l = 1 : length(dimsF)

            for j = S(:,:)
    
                Lie = Lie - TrFind(X1, dimsX1, j) * F(l,i - j) * (preDefineJ) * j(l) + TrFind(X1, dimsX1, (i - j)) * (preDefineJ) * j(l) * F(l,j);
            
            end
            
        end
    
end
