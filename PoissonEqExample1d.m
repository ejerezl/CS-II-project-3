% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
%k1 = @(x) 1+(x-1).*(1/4-1);
k1 = @(x) 1;

% Exact solution
u_fabricated = @(x) sin(2*pi*x);

%RHS function
f = @(x) 8*pi*pi*sin(2*pi*x).*k1(x)-(2*pi*cos(2*x*pi).*(1/4-1));

%Determine number of cells
nx = 17;
dx = 1/(nx-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices1d(nx,f,k1);

% Solve linear system
u = A\b;

% Reshape and plot
X = (0:dx:1);
figure(1)
subplot(1,2,1)
plot(X,u)

% Calculate fabricated solution vector and compute error in "L2-norm"
u_fabricated_vect = u_fabricated(cells(:,1));
error = norm(u-u_fabricated_vect,2)*sqrt(dx);
xlabel('x')
title('u')
colorbar()


% Compute the flux vector [F_1,F_2] at the cells and plot with quiver
[F_1] = flux1d(G,K,u,cells,nx);
subplot(1,2,2)
quiver(cells(:,1),F_1)
xlabel('x')
title('Flux')

% Plot the fabricated solution
U_fabricated = u_fabricated_vect;
figure(2)
subplot(1,1,1)
plot(X, U_fabricated)
xlabel('x')
title('u_{fabricated}')
colorbar()