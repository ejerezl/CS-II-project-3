
% Poisson equation: -Div(k(x,y)Grad(u)) = f

%Diffusivity
%k1 = @(x) 1+(x-1).*(1/4-1);
k1 = @(x) 1;

% Exact solution
%u_fabricated = @(x) sin(2*pi*x);
u_fabricated = @(x) -0.5*x.^2 + 0.5*x;

%RHS function
%f = @(x) 4*pi*pi*sin(2*pi*x).*k1(x); %-(2*pi*cos(2*x*pi).*(1/4-1));
f = @(x) 1;

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
error = norm(u-u_fabricated_vect,2)*sqrt(dx)
xlabel('x')
title('u')
<<<<<<< HEAD
=======
colorbar()

energyError = (energyError1d(u, -(K*G)*u, 1, dx, cells, nx))^(1/2)
>>>>>>> c94fd49664fa6790b9d39b8ff665b7936bd309c9

%Compute and plot flux vector
X = zeros(nx+1,1);
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
q = flux1d(k1,u,dx,nx);
subplot(1,2,2)
plot(X,q)
xlabel('x')
title('Flux')


% Plot the fabricated solution
X = (0:dx:1);
U_fabricated = u_fabricated_vect;
figure(2)
subplot(1,1,1)
plot(X, U_fabricated)
xlabel('x')
title('u_{fabricated}')


[err_true,err_calc, err_v] = energy_error_norm1(u, q, nx, edges, k1);
disp(err_true)
disp(err_calc)
disp(err_true/err_v)
disp(err_calc/err_v)

% Runge Kutta type 1 error estimate

%2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
[err1, err2] = rktype1(nx, 1.4, 8); 
disp(err1);
disp(err2);

%Convergence rate
convergence_error = zeros(4,1);
grid_refinery = zeros(4,1);
for i=1:4
    [conv_rate, grid_size] = convergence1(nx);
    err1 = [conv_rate, grid_size];
    convergence_error(i) = err1(1);
    grid_refinery(i) = err1(2);
    nx = 2*nx;
end

loglog(grid_refinery, convergence_error)
title('Convergence Plot');
hold on
y=grid_refinery.^(1.4);
loglog(grid_refinery,y)
hold off
legend({'numerical scheme','10^{1.4}'}, 'Location','northwest')
xlabel('gridsize')
ylabel('error')
xlim([0 0.07])



