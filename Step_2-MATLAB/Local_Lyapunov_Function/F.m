

function y = F(l,i)

   
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

    FF = [FF, F2];

    % In FF, columns from (1) to (n, dimension of F) holds the non-zero
    % coefficient indexes,
    % columns form (n+1) to (n+n) represents the non-zero coefficients


    
    % we indexed the non-zero coefficienttics in FF.
    % if " i th " coefficient is non-zero, return it in y.
    % return zero for zero coefficients.
    if (isempty(find(ismember(FF(:,1:2), i', 'rows'))) && isempty(find(ismember(FF(:,1:2), -i', 'rows'))))
        y = 0; % zero coefficients,
    elseif ~isempty(find(ismember(FF(:,1:2), i', 'rows')))
        % return non zero coeffs
        y = FF(find(ismember(FF(:,1:2), i', 'rows')),2+l);
    else
        % return conjugates of the coefficients for non-existing indexes.
        y = conj(FF(find(ismember(FF(:,1:2), -i', 'rows')),2+l));

    end
   
   

   



end