function [A, b, G, D, K, dx, X, E, DX_u, DX_q] = assemble1d(num_cells, f, k)
% assemble 1D matrices. Boundary is NOT part of the mesh
%
%	A stiffnes Matrix
%	b rhs
%	G gradient Operator
%	D divergence Operator
%	K ks s.t. -(D*K*G)*u = b
%	dx grid size
%	X cell centers of inner cells
%	E boundary positions of inner cells without domain boundary
%	DX_u central difference op for u. This is equivalent to forward diff
% of the linear over the cell boundaries interpolated u
%	DX_q forward diff of q on the inner cells. This is the div operator!

N = num_cells;
dx = 1 / (N + 1);
% grid
X = linspace(dx, 1-dx, N);
% edge grid
E = linspace(dx/2, 1 - dx/2, N + 1);


% G grad op
G = spdiags(1/dx * ones(N, 1)*[-1 1], [-1, 0], N + 1, N);

D = -transpose(G);
% TODO fancy formeln von Erlend/Kaja?
%K = spdiags(transpose(k(E)), 0, N + 1, N + 1);

ls = k([0, X]);
rs = k([X, 1]);
K = spdiags(transpose(2 ./ (1 ./ ls + 1 ./ rs)), 0, N + 1, N + 1);

A = -D * K * G;
b = transpose(f(X));

DX_u = spdiags(1/(2*dx) * ones(N, 1)*[-1 1], [-1, 1], N, N);

DX_q = D;

end
