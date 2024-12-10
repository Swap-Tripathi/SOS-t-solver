function Lie = LieFind(X1, dimsF, dimsX1, i)
    
    % dims = [2,3,5], metinde m_v = [1,2,4]
    
    % Örnek kronecker product of [x1, x2] * [y1 y2 y3] * [z1 z2 z3 z4 z5]^T
    % [x1y1z1 x1y1z2 x1y1z3 x1y1z4 x1y1z5 x1y2z1 x1y2z2 x1y2z3 (i = 7)

    % Sırf upper-diagonal trace'leri düşünerek trace elde edilebilir mi?

    dimsX2 = dimsX1 + dimsF - 1; 

    d = length(dimsX2); 
    bb = max(dimsX2); 
    B = cell(1,d);
    [B{:}] = ndgrid(-bb+1:bb-1);% ALKIM CHECK
    B = cellfun(@(M) M(:), B, 'uniform', 0);
    S = [B{end:-1:1}]';
    Lie = 0;

    preDefineJ = sqrt(-1);
    
        for l = 1 : length(dimsF)

            for j = S(:,:)
    
                % Lie = Lie - TrFind(X1, dimsX1, TrX1, (i - j)) * F(l,i - j) * (sqrt(-1)) * j(l) + TrFind(X1, dimsX1, TrX1, j) * (sqrt(-1)) * j(l) * F(l,j);
                Lie = Lie - TrFind(X1, dimsX1, j) * F(l,i - j) * (preDefineJ) * j(l) + TrFind(X1, dimsX1, (i - j)) * (preDefineJ) * j(l) * F(l,j);
            
            end
            
        end
    
end