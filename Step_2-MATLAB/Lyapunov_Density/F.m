

function y = F(l,i)

    % It should be symmetric with respect to reflection about the middle
    % point of the matrix.
    % -J = J*

    %     F1(theta1, theta2)
%                  c-2          c-1                 c0          
    % f1 = (      [0,           0,                  1j*a,       0,                  0.5j*a,...       
    %              0,           0,                  -1j*cos(b), -0.5j*exp(1j*b),    0,...
    %              0.5*1j*a,    -0.5j*exp(-1j*b),   0,          0.5j*exp(1j*b),     -0.5j*a,...
    %              0,           0.5j*exp(-1j*b),    1j*cos(b),  0,                  0,...
    %              -0.5j*a,     0,                  -1j*a,      0 ,                 0         ]);

global a;
global b;

    
    F2 = [0.5j*a;...
          -0.5j*a;
          -0.5j*exp(-1j*b);
          0.5j*exp(-1j*b);
                1j*a;
                -1j*cos(b)];
    
    FF   = [-2  0 1j*a;...
            -2  2 0.5*1j*a;
            -1  0 -1j*cos(b);
            -1  1 -0.5j*exp(1j*b);
             0  -2 0.5*1j*a;
             0  -1 -0.5j*exp(-1j*b)];
    
    FF = [FF, F2];                                      % Adding F2 as a columns in previous line FF which only contained information of non zero elements of F1
                                                        % FF can also be obtained using mathematica (STEP 1)

    % In FF, columns from (1) to (n, dimension of F) holds the non-zero
    % coefficient indexes,
    % columns form (n+1) to (n+n) represents the non-zero coefficients


    
    % we indexed the non-zero coefficients in FF.
    % if " i th " coefficient is non-zero, return it in y.
    % return zero for zero coefficients.
    
    % When A is a matrix and B is a row vector, ismember(A,B,"rows") compares each row of A with B and returns a vector
    % with 1 in the position where the match is made. Otherwise, gives a zero matrix.

    % Find of the above ismember vector returns the indices having non-zero
    % entries.
    
    % FF(:,1:2) creates a submatrix from FF containing all indices of first
    % dimension and indices 1,2 for second dimension. Thus, the submatrix
    % containing position of nonzero elements

    if (isempty(find(ismember(FF(:,1:2), i', 'rows'))) && isempty(find(ismember(FF(:,1:2), -i', 'rows'))))  % if both subscript vectors i and -i are not in subscript part of FF, then the coefficient is 0 vector
        y = 0;
    elseif ~isempty(find(ismember(FF(:,1:2), i', 'rows')))                                                  % if subscript vector i is an entry in FF, then return the coefficient
        y = FF(find(ismember(FF(:,1:2), i', 'rows')),2+l); 
    else
        y = conj(FF(find(ismember(FF(:,1:2), -i', 'rows')),2+l));                                           % if negative of subscript vector i is an entry in FF, then return the conjugate of the coefficient
                                                                                                            % ALL THE 2s IN THE IF LOOP CHANGE TO N when F IS N-DIMENSIONAL
    end
    
end