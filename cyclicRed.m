function [ xSol ] = cyclicRed( M, xInd, b, xSol )
    %% Use cyclic reduction to recursively solve M * x = b by splitting it into smaller matrices
	% M is the matrix to be solved for x in M * x = b
	% xInd is a vector to track the positions of the solutions "x" when they are permutated
	% b is the vector containing the result of the equation M * x = b
	% x is the vector containing the solution to the equation M * x = b

	[n, n] = size(M);
    
	P = eye(n);

	odd = 1 : 2 : n;
	even = 2 : 2 : n;

	P = [ P(:, odd), P(:, even)];
   
	% Permutate the matrix M, indX, and bSol
	M = P' * M * P;
	xInd = P' * xInd;
	b = P' * b;

	% Split the matrix M into four blocks
	k = n/2;

	A = [ M(1:k, 1:k)];
	B = [ M(1:k, k+1:n)];
	C = [ M(k+1:n, 1:k)];
	D = [ M(k+1:n, k+1:n)];

	% Assemble the Schur complements for the odd and even components of M
	Modd = A - B * inv(D) * C; 
	Meven = D - C * inv(A) * B;

	% Update bSol
    bUpdate = [ eye(k,k), -B * inv(D); -C * inv(A), eye(k, k) ];
	b = bUpdate * b;
    
	if n == 2 % base case
		xSol(xInd(1:k)) = inv(Modd) * b(1:k);
        xSol(xInd(k+1:n)) = inv(Meven) * b(k+1:n);
	else
		xSol = cyclicRed( Modd, xInd(1:k), b(1:k) , xSol );
		xSol = cyclicRed( Meven, xInd(k+1:n), b(k+1:n), xSol );
	end
end