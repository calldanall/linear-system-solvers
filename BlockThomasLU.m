function [ G ] = BlockThomasLU(M, B, C, G, m, n)

	%% FORWARD SUBSTITUTE %%  

	%Solve the equation L * X = G
	for i = 2 : m; % Loop through the blocks "m"
	    fac = (i-1)*n; % Create a factor to allow manipulation of block "m"
        
	    G(fac+1 : fac+n) = G(fac+1 : fac+n) - M(:, fac+1-n : fac) * G(fac+1-n : fac); % Subtract M_m-1 * X_m-1 from both sides
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