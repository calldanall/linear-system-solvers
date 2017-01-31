function [ G, A, B ] = BlockThomas(a, Ba, Bb, Bc, c, G, n, m)
% This function solves Main * X = G for X using the Thomas algorithm , 
% where Main = tridiag(A,B,C). The solution X is returned on exit.

% n is the dimension of the square blocks and m is the number of blocks

%% Create and Populate all blocks A_m, B_m, C_m %%

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
    B(:, fac+1 : fac+n) = diag(Bb(fac+1 : fac+n)); % Populate the main 
                                                   % diagonal of B

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

%% FORWARD SWEEP %%
for i = 1 : m - 1 % Loop through each block
    fac = (i-1)*n; % Create a factor to allow manipulation of block "m"

    A_trans = A(:, fac+1 : fac+n)'; % Creates A^T
    B_trans = B(:, fac+1 : fac+n)'; % Creates B^T

    [M_trans, dummy] = Gauss(B_trans, A_trans, n); % Call function "Gauss" 
                                                   % from NR

    M = M_trans';
    A(:, fac+1 : fac+n) = M;

    B(:, fac+n+1 : fac+2*n) = B(:, fac+n+1 : fac+2*n) - M * C(:, fac+1 : fac+n); % Modify B_m+1
    G(fac+n+1 : fac+2*n) = G(fac+n+1 : fac+2*n) - M * G(fac+1 : fac+n); % Modify G_m+1
end

%% BACK SUBSTITUTE %%  

for i = m : -1: 1; % Loop through the blocks "m" in reverse order
    fac = (i-1)*n; % Create a factor to allow manipulation of block "m"
    
    if i == m
        [G(fac+1 : fac+n), dummy] = Gauss(B(:, fac+1 : fac+n), G(fac+1 : fac+n), n); % Solve for B_m
    else
        G(fac+1 : fac+n) = G(fac+1 : fac+n) - C(:, fac+1 : fac+n) * G(fac+1+n : fac+2*n); % Subtract C_m * X_m+1 from both sides
        [G(fac+1 : fac+n), dummy] = Gauss(B(:, fac+1 : fac+n), G(fac+1 : fac+n), n);
    end
end
end