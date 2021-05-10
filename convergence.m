function [err_conv, grid_size] = convergence(nx)
% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
k = @(x,y) 1+(x-1).*(y-1);

% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);

%RHS function
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
%f = @(x,y) 1;

%Determine number of cells in each direction
ny = nx;
dx = 1/(nx-1);
dy = 1/(ny-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);

% Solve linear system
u = A\b;

u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
err_conv = norm(u-u_fabricated_vect,2)*sqrt(dx*dy);
grid_size = dx
