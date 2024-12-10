function [status, V, W, X, Y] = solveSDPAlternative(dimsF, dimsGSV,r)
    
    dimsGSW = dimsGSV + dimsF - 1; % GSV is sq matrix of size (n_s^v+1)^d 
                                   % dimsGSV=[n_s^v+1,...,n_s^v+1]
                                   % GF is sq matrix of size (n_f+1)^d
                                   % dimsF=[n_f+1,....,n_f+1]
                                
                                   % GSW is sq matrix of size (n_s^v+n_f+1)^d
                                   % dimsGSW=[n_s^v+n_f+1,....,n_s^v+n_f+1]
    
    cvx_begin sdp quiet
    cvx_precision low

    variable GS0W(prod(dimsGSW), prod(dimsGSW)) hermitian
    variable GS1W(prod(dimsGSW), prod(dimsGSW)) hermitian
    variable GS0V(prod(dimsGSV), prod(dimsGSV)) hermitian
    variable GS1V(prod(dimsGSV), prod(dimsGSV)) hermitian
    
    d = length(dimsGSW);
    bb = max(dimsGSW)+r;              % max(dimsGSW)=n_s^v+n_f+1+r
    B = cell(1,d);
    [B{:}] = ndgrid(-bb+1:bb-1); 
    B = cellfun(@(M) M(:), B, 'uniform', 0);
    S = [B{end:-1:1}]';                             % All the steps in total produce all possible index vectors for the trig poly corresponding to S_l^w as columns
    

    
        

    minimize(0);
    subject to
 
        (GS0W)>=0;
        (GS1W)>=0;
        (GS0V)>=0;
        (GS1V)>=0;

        %% LOCAL PROBLEM

        for k = S(:,1:(end+1)/2)
           lhs = 0;
           
           % Calculate the first term: Tr[T^{nf+nsv}_k G_W]
           % when d=2, take p=1 and q=2 in (26) and m comes from the set
           % {[r,-r],[-r,r]}
           lhs = lhs + (TrFind(GS0W, dimsGSW, k))+ 0.5 *(TrFind(GS1W, dimsGSW, k-[-r,r])+TrFind(GS1W, dimsGSW, k-[r,-r]));
           
           % Define the right-hand side as per equation (26)
           rhs = 0;
           
           for l = 1:d
               for j = S(:,:)
                   % Compute the RHS of eq 26 as per the summation
                   rhs = rhs - 1i * j(l)* F(l,k-j) * ( TrFind(GS0V, dimsGSV, j) + 0.5 * (TrFind(GS1V,dimsGSV,j-[r,-r]) + TrFind(GS1V,dimsGSV,j-[-r,r]) ) );
               end
           end
       
           % Now add the constraint to the CVX problem
           lhs == rhs;
        end
        %% LOCAL PROBLEM RELATED PART ENDS

        trace(GS0V) >= 2;
        trace(GS1V) >= 2;
        trace(GS0W) >= 2;
        trace(GS1W) >= 2;
        sum(GS0V(:)) == 0; %S0v(theta)=0 at origin
        sum(GS1V(:)) == 0; %S1v(theta)=0 at origin

    cvx_end

    status = cvx_status;
    V = GS0V;
    W = GS1V;
    X = GS0W;
    Y = GS1W;
end