function [status, V, W] = solveSDPAlternative(dimsF, dimsX1)
    

    % dimsF => number of columns represents the dimension of the vector
    % field F(theta), and each entry of dimsF represents the number of
    % harmonics of corresponding dimension.

    % dimsX1 => number of columns represents the dimension of the lyapunov
    % density (or lyapunov function), V(theta), ann each entry of the dimsX1
    % represents the number of harmonics of V(theta).


    % thanks to symmetry,
    % dimsF => n_f + 1 ( normally -nf, -nf +1, .... 0, 1, .. nf -1, nf)
    % (Eq.7)

    % dimsX1 => n_v + 1 (Eq. 7)
    % dimsX2 => n_w + 1 (Eq. 7)
    dimsX2 = dimsX1 + dimsF - 1; % (n_v + 1) + (n_f + 1) -1
                                 % in general dims are vectors eg. [3,3]+[3,3]-1=[5,5]
    
    cvx_clear % to clear the cash memory inculding previous iteration info.
    cvx_begin sdp quiet % initialize the cvx as sdp, "quite" for no debug output
    cvx_precision low

    % define cvx variables to represent sdp vars.
    % cvx variables cannot be assigned to a non-cvx variable (e.g. double)
    % output of each cvx variable including function are definde as cvx
    % var.

    variable X1(prod(dimsX1), prod(dimsX1)) hermitian
    variable X2(prod(dimsX2), prod(dimsX2)) hermitian   

    % to get each combination of Fourier coefficient indexes of "W"
    % e.g. [-5,-5,-5], [-5,-5,-4], ...., [-5,-5,5], ..., [5, 5, 5]
    % but, due to the symmetry, [0,0,0] is assumed as last coefficient index.
    d = length(dimsX2);                                                     % number of Harmonics in W
    bb = max(dimsX2)-1;                                                     % largest degree of a Harmonic in W minus 1
    B = cell(1,d);                                                          % B is a 1*d cell array of empty matrices
    [B{:}] = ndgrid(-bb:bb);                                                % writing -bb,...0,....,bb as a column vector B
                                                                            % B is a 1*d cell array which each cell as a (2*bb+1) x ...d times.. x (2*bb+1) hypercube matrix 
    B = cellfun(@(M) M(:), B, 'uniform', 0);                                % B is now written as 1*d cell array with each cell being a vector of length (2*bb+1)^d
                                                                            % cellfun(func,C) applies the function func to the contents of each cell of cell array C, one cell at a time. 
                                                                            % @(M) M(:) returns the hypercube matrix as vectors as stated above
                                                                            % 'uniform',0 is required whenever function outputs a non-scalar value (here it inputs a matrix and outputs a vector)
    S = [B{end:-1:1}]';                                                     % Cell has indices 1 to end. This command reverses the order of the entries
                                                                            % All the steps in total produce all possible index vectors for the trig poly corresponding to X2 (or V)
    

    minimize(0);
    subject to
 
        (X1)>=0; % positivity of V, sos(V)
        (X2)>=0; % positivity of W, sos(W)
        for i = S(:,1:(end+1)/2) % a for loop holding a vector              % Each column of S is an index for subscript in a specific order. There are (2*bb+1)^d columns (odd). 
                                                                            % Since only half of the indices are relevant we take columns till the middle part of the matrix S. And all rows are chosen using ":" in first component.
            % define Lie derivative, and traceVectors holding coefficients.
            expression Lie
            expression TrX1
            expression TrX2
            LieFind(X1, dimsF, dimsX1, Lie, i, TrX1) == ( TrFind(X2, dimsX2, TrX2, i) );
        end
        % positivity of constant coefficients of V(theta) and W(theta)
        trace(X1) >= 2;
        trace(X2) >= 2;
        sum(X1(:)) == 0; % lyapunov function should be zero at the origin.

    cvx_end

    status = cvx_status;
    V = X1;
    W = X2;
end