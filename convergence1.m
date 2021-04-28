%Convergence of 1d method for f=1

function [err_conv, grid_size] = convergence1(nx)

% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
%k1 = @(x) 1+(x-1).*(1/4-1);
k1 = @(x) 1;


%RHS function
%f = @(x) 4*pi*pi*sin(2*pi*x).*k1(x); %-(2*pi*cos(2*x*pi).*(1/4-1));
f = @(x) 1;


%Determine number of cells
dx = 1/(nx-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices1d(nx,f,k1);

% Solve linear system
u = A\b;

%Compute flux vector
X = zeros(nx+1,1);
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
q = flux1d(k1,u,dx,nx);

[err_true,err_calc, err_v] = energy_error_norm1(u, q, nx, edges, k1);

err_conv = err_true;
grid_size = dx;
end