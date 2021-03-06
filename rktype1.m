%%%%%%%%% Runge Kutta Type 1 estimate %%%%%%%%%%

%a = convergence rate

function [err1, err2] = rktype1(nx, a, k)
% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
k1 = @(x) 1;
k_const = true;


%RHS function
f = @(x) 1;
f_const = true;

%%% u_h1, h1 = dx %%%
%Determine number of cells
nx1 = nx;
dx = 1/(nx-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices1d(nx,f,k1);

% Solve linear system
u1 = A\b;

%Compute flux vector
X = zeros(nx+1,1);
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
q = flux1d(k1,u1,dx,nx, k_const, f_const);

grid1 = zeros(nx1+1,1);
grid1(2:(nx1)) = edges(:,1,1);
grid1(end) = 1;

[err_true,err_calc, err_v] = energy_error_norm1(u1, q, nx1, edges, k1);
err1 = err_true;

%%% u_h2, h2 = k*h1 %%%
nx2 = nx*k;
%Determine number of cells
dx2 = 1/(nx2-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices1d(nx2,f,k1);

% Solve linear system
u2 = A\b;

%Compute flux vector
X = zeros(nx2+1,1);
X(1) = 0;
X(end) = 1;
for i=2:nx2
    X(i) = 0.5*dx2 + dx2*(i-2);
end
q = flux1d(k1,u2,dx2,nx2, k_const, f_const);

%%%%%% Calculate energy norm of u_h1-u_h2 %%%%%

%define new grid of edges which includes 0 and 1
grid2 = zeros(nx2+1,1);
grid2(2:(nx2)) = edges(:,1,1);
grid2(end) = 1;

%define potential for u1 on grid points as averages of the adjacents cells
pot1 = zeros(size(grid1));
pot1(1) = 0; %Boundary condition
for i=2:size(grid1)-1
    pot1(i) = 0.5*(u1(i)+u1(i-1));
end
pot1(end) = 0; %Boundary condition

%define potential on grid points as averages of the adjacents cells
pot2 = zeros(size(grid2));
pot2(1) = 0; %Boundary condition
for i=2:size(grid2)-1
    pot2(i) = 0.5*(u2(i)+u2(i-1));
end
pot2(end) = 0; %Boundary condition

err2 = 0;

for i=1:nx
    a1_1 = (pot1(i+1)-pot1(i))/(grid1(i+1)-grid1(i));
    for j=1:k
        a1_2 = (pot2((i-1)*k + j + 1) - pot2((i-1)*k + j))/(grid2(i+1)-grid2(i));
        %ki = (k(grid2(i)))^(-1);
        ki = 1;
        err = (a1_1 + a1_2)^2*((grid2(i+1)-grid2(i)));
        err2 = err2 + ki*err;
    end
end

% %left boundary
% i=1;
% a1_1 = (pot1(i+1)-pot1(i))/(dx/2);
% for j=2:k
%     a1_2 = (pot2((i-1)*k + j + 1) - pot2((i-1)*k + j))/(dx2/2);
%     %ki = (k(grid2(i)))^(-1);
%     ki = 1;
%     err = (a1_1 + a1_2)^2*(dx2/2);
%     err2 = err2 + ki*err;
% end
% j = 1;
% a1_2 = (pot2((i-1)*k + j + 1) - pot2((i-1)*k + j))/(dx2/2);
% %ki = (k(grid2(i)))^(-1);
% ki = 1;
% err = (a1_1 + a1_2)^2*(dx2/2);
% err2 = err2 + ki*err;
% 
% %right boundary
% i=nx;
% a1_1 = (pot1(i+1)-pot1(i))/(dx/2);
% for j=1:k-1
%     a1_2 = (pot2((i-1)*k + j + 1) - pot2((i-1)*k + j))/(dx2/2);
%     %ki = (k(grid2(i)))^(-1);
%     ki = 1;
%     err = (a1_1 + a1_2)^2*(dx2/2);
%     err2 = err2 + ki*err;
% end
% j = k;
% a1_2 = (pot2((i-1)*k + j + 1) - pot2((i-1)*k + j))/(dx2/2);
% %ki = (k(grid2(i)))^(-1);
% ki = 1;
% err = (a1_1 + a1_2)^2*(dx2/2);
% err2 = err2 + ki*err;


err2 = sqrt(err2);

err2 = 1/(k^a - 1)*err2;
        
        
        