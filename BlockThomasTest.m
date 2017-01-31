%% Block Thomas Test Script %%
%% Verify that the BlockThomas and BlockThomasLU functions work by comparing them to the verified Gauss function
% Generate the input parameters for the function "BlockThomas"

    n = 5; % dimensions of the square block matrices
    m = 3; % number of B blocks in the full matrix

    a = randi(10, 1, n*(m-1)); % Generate the "A" diagonal

    Ba = randi(10, 1, m*n-1); % Generate the "Ba" sub-diagonal
    Bb = randi(10, 1, n*m);     % Generate the "Bb" diagonal
    Bc = randi(10, 1, m*n-1); % Generate the "Bc" super-diagonal

    c = randi(10, 1, n*(m-1)); % Generate the "C" diagonal
    
    G = randi(10, n*m, 1);

    % Create the blocks to be used to construct the main matrix

    A = zeros(n,n*(m-1));
    B = zeros(n,n*m);
    C = zeros(n,n*(m-1));

for i = 1 : m % Loop through each block
    fac = (i-1)*n; % Create a factor to allow manipulation of block "m"

    % Populate the blocks A and C
    if i ~= m
        C(:, fac+1 : fac+n) = diag(c(fac+1 : fac+n));
        A(:, fac+1 : fac+n) = diag(a(fac+1 : fac+n));
    end

    % Populate block B
    B(:, fac+1 : fac+n) = diag(Bb(fac+1 : fac+n)); % Populate the main diagonal of B

    % Populate the super and sub-diagonals of B
    for k = 1 : n % Loop through the columns of block m of B
        for j = 1 : n % Loop through the rows of block m of B
            if k - j == 1 % If true, insert a super-diagonal element
                B(j, fac+k) = Bc(fac+j);
            elseif j - k == 1 % If true, insert a sub-diagonal element
                B(j, fac+k) = Ba(fac+k);
            end
        end
    end
end

    Main = zeros(m*n,m*n);

for i = 1 : m
    fac = (i-1)*n;
    Main(fac+1 : fac+n, fac+1 : fac+n) = B(:, fac+1 : fac+n);
    
    if i ~= m
    Main(fac+1 : fac+n, fac+n+1 : fac+2*n) = C(:, fac+1 : fac+n);
    Main(fac+n+1 : fac+2*n, fac+1 : fac+n) = A(:, fac+1 : fac+n);
    end
end


% Use the function "BlockThomas" to find the solution set X
[ X, M, B ] = BlockThomas(a, Ba, Bb, Bc, c, G, n, m);

% Use the function "Gauss" to find the solution set X and store it in a vector Y
[Y,dummy] = Gauss(Main,G,m*n);

% Use the function "BlockThomasLU" to find the solution set X and store it in a vector Y
[ Z ] = BlockThomasLU(M, B, C, G, m, n);

% Since my version of Matlab only has a round function to round to the nearest integer...
% I will round X-Y to the nearest integer and if that is equal to 0 I will say they are equal.

% Check to see the solutions are the same
if abs(X-Y) < 1e-9
    fprintf( 'The function "BlockThomas" works! :)\n' );
else
    fprintf( 'The function "BlockThomas" doesn''t work :(\n' );
end

if abs(Z-Y) < 1e-9
    fprintf( 'The function "BlockThomasLU" works! :)\n' );
else
    fprintf( 'The function "BlockThomasLU" doesn''t work :(\n' );
end