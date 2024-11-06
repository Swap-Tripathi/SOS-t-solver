

% compute the trace of a product of a matrix x 
% with a tensor product of Toeplitz matrices. 

function TrX = TrFind(x, n, TrX, k)
    %producing the required tensor product of Toeplitz matrices
    y = 1;
    for i = 1 : length(k)
        if(abs(k(i)) > n(i)) % (offset from the diagonal) specified by k(i) is too large for the size of the matrix n(i) x n(i).
            y = kron(y, zeros(n(i), n(i)));
        else
            y = kron(y,diag(ones(n(i) - abs(k(i)), 1), k(i))); % toeplitz matrix is created, T_j
        end
    end
    TrX = TrX + trace(y' * x); % trace(T_k * R) in paper (https://www.preprints.org/manuscript/202409.1361/v1).

end