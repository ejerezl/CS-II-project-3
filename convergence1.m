%%%%% Convergence of 1d method for f=1 and k=1 %%%%%

function [err_conv, grid_size, err_calc, err_pot, err_v] = convergence1(nx)

% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
k1 = @(x) 1;
k_const = true;


%RHS function
f = @(x) 1;
f_const = true;


%Determine number of cells
dx = 1/(nx-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices1d(nx,f,k1);

% Solve linear system
u = A\b;

%Compute flux vector
q = flux1d(k1,u,dx,nx, k_const, f_const);

[err_true,err_calc, err_v, err_pot] = energy_error_norm1(u, q, nx, edges, k1);

err_conv = err_true;
grid_size = dx;
end