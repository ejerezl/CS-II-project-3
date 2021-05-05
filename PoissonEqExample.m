% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
%k = @(x,y) 1+(x-1).*(y-1);
k = @(x,y) 1;
% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);

%RHS function
%f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
f = @(x,y) 1;

%Determine number of cells in each direction
nx = 17;
ny = 17;
dx = 1/(nx-1);
dy = 1/(ny-1);

% Build matrices
[A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);

% Solve linear system
u = A\b;

% Reshape and plot
[X,Y] = meshgrid(0:dx:1,0:dy:1);
% to have coordinates x and y separated
U = reshape(u,nx,ny);
figure(1)
%gcf returns the current figure it would be the same to do figuree=fig(1) and
%then in the set function instead of gcf write figure
%position params: [left, bottom, width, height]
set(gcf,'Position',[100 100 1200 500])
%subplot(m,n,p) divide la figura actual en una cuadrícula de m por n y crea ejes en la posición especificada por p
subplot(1,2,1);
surf(X,Y,U)
axis([0 1 0 1])

% Calculate fabricated solution vector and compute error in "L2-norm"
u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
error = norm(u-u_fabricated_vect,2)*sqrt(dx*dy)
xlabel('x')
ylabel('y')
title('u')
colorbar()

% Compute error without knowing real solution
[energy_error, conservation_integral] = energyError(U,-(K.*G)*u, 1, dx, dy, cells,nx, ny)
majorant = energy_error + (1/(sqrt(1)*pi))*conservation_integral
energy_norm = norm(G*(u-u_fabricated_vect),2)*sqrt(dx*dy)
relative_error = energy_error/norm(G*u,2)*sqrt(dx*dy)


% Compute the flux vector [F_1,F_2] at the cells and plot with quiver
[F_1, F_2] = flux(G,K,u,cells,nx,ny);
subplot(1,2,2)
quiver(cells(:,1),cells(:,2),F_1,F_2)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Flux')

% Plot the fabricated solution
U_fabricated = reshape(u_fabricated_vect,nx,ny);
figure(2)
subplot(1,1,1)
surf(X,Y,U_fabricated)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('u_{fabricated}')
colorbar()


%Convergence rate
convergence_error = zeros(4,1);
grid_refinery = zeros(4,1);
for i=1:4
    [conv_rate, grid_size] = convergence(nx);
    err1 = [conv_rate, grid_size];
    convergence_error(i) = err1(1);
    grid_refinery(i) = err1(2);
    nx = 2*nx;
end

loglog(grid_refinery, convergence_error)
title('Convergence Plot');
hold on
y=grid_refinery;
loglog(grid_refinery,y)
hold off
legend({'numerical scheme','10^1'}, 'Location','northwest')
xlabel('gridsize')
ylabel('error')
xlim([0 0.07])