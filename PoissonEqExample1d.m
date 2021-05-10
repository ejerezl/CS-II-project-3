
% Poisson equation: -Div(k(x,y)Grad(u)) = f

%%%%% Choose set of initial data %%%%%

k1 = @(x) 2 - x;
u_fabricated = @(x) sin(2*pi*x);
q_fabricated = @(x) -2*pi*cos(2*pi*x).*(2-x);
f = @(x) 4*pi*pi*sin(2*pi*x).*k1(x) + 2*pi*cos(2*pi*x);
k_const = false;
f_const = false;

% k1 = @(x) 1;
% u_fabricated = @(x) sin(2*pi*x);
% q_fabricated = @(x) -2*pi*cos(2*pi*x);
% f = @(x) 4*pi*pi*sin(2*pi*x).*k1(x);
% k_const = true;
% f_const = false;

% k1 = @(x) 1;
% u_fabricated = @(x) -0.5*x.^2 + 0.5*x;
% q_fabricated = @(x) x - 0.5;
% f = @(x) 1;
% k_const = true;
% f_const = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%Compute and plot flux vector
X = zeros(nx+1,1); %new grid
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
q = flux1d(k1,u,dx,nx, k_const, f_const);
subplot(1,2,2)
plot(X,q)
xlabel('x')
title('Flux')

% Plot the fabricated solution
X = (0:dx:1);
U_fabricated = u_fabricated_vect;
figure(2)
subplot(1,2,1)
plot(X, U_fabricated)
xlabel('x')
title('u_{fabricated}')

%plot fabricated flux
X = zeros(nx+1,1); %new grid
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
q_fabricated_vect = q_fabricated(X);
Q_fabricated = q_fabricated_vect;
figure(2)
subplot(1,2,2)
plot(X,Q_fabricated)
xlabel('x')
title('q_{fabricated}')

%%%%% Error calculation & Convergence %%%%%

if k_const && f_const
    %Energy error
    [err_true,err_calc, err_v, err_pot] = energy_error_norm1(u, q, nx, edges, k1);

    % Runge Kutta type 1 error estimate

    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1(nx, 1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);

    %Convergence rate
    convergence_error = zeros(4,1);
    grid_refinery = zeros(4,1);
    for i=1:4
        [conv_rate, grid_size, err_calc, err_pot, err_v] = convergence1(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('calculated error: %d.\n', err_calc);
        fprintf('error of potential: %d.\n', err_pot);
        fprintf('relative true error: %d.\n', conv_rate/err_v);
        fprintf('relative calculated error: %d.\n', err_calc/err_v);
        fprintf('efficiency index: %d.\n', (err_calc)/conv_rate)
    end
    
    %Plot convergence
    figure(3)
    loglog(grid_refinery, convergence_error)
    title('Convergence Plot');
    hold on
    y=grid_refinery.^(1);
    loglog(grid_refinery,y)
    loglog(grid_refinery, convergence_error, 'b.', 'Markersize', 20)
    hold off
    legend({'numerical scheme','10^{1}'}, 'Location','northwest')
    xlabel('gridsize')
    ylabel('error')
    xlim([0 0.07])

elseif k_const 
    c = 1/pi; %Poincare constant for k=1
    %Energy error
    [err_true,err_calc, err_cons, err_v, err_pot] = energy_error_norm2(u, q, nx, edges, k1);
    
    %Runge Kutta estimate type 1
    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1_2(nx, 1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);
    
    %Convergence rate
    convergence_error = zeros(4,1);
    grid_refinery = zeros(4,1);
    for i=1:4
        [conv_rate, grid_size, err_calc, err_cons, err_pot, err_v] = convergence2(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('calculated error: %d.\n', err_calc);
        fprintf('error of potential: %d.\n', err_pot);
        fprintf('conservation error %d.\n', c*err_cons);
        fprintf('relative true error: %d.\n', conv_rate/err_v);
        fprintf('relative calculated error: %d.\n', (err_calc+c*err_cons)/err_v);
        fprintf('efficiency index: %d.\n', (err_calc+c*err_cons)/conv_rate)
        fprintf('combined efficiency index: %d.\n', 3*(err_calc+c*err_cons)/(conv_rate + err_pot + c*err_cons))
    end
    
    %Plot convergence rate
    figure(3)
    loglog(grid_refinery, convergence_error)
    title('Convergence Plot');
    hold on 
    y=grid_refinery.^(1);
    loglog(grid_refinery,y)
    loglog(grid_refinery, convergence_error, 'b.', 'Markersize', 20)
    hold off
    legend({'numerical scheme' ,'10^{1}'}, 'Location','northwest')
    xlabel('gridsize')
    ylabel('error')
    xlim([0 0.07])
else
    c=1/sqrt(14.3249); %14.3249 is smallest eigenvalue of A which is not 1
    [err_true, err_calc, err_cons, err_v, err_pot] = energy_error_norm3(u, q, nx, edges, k1);
    
    %Runge Kutta estimate type 1
    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1_3(nx, 1.1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);
    
    %Convergence rate
    convergence_error = zeros(4,1);
    grid_refinery = zeros(4,1);
    for i=1:4
        [conv_rate, grid_size, err_calc, err_cons, err_pot, err_v] = convergence3(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('calculated error: %d.\n', err_calc);
        fprintf('error of potential: %d.\n', err_pot);
        fprintf('conservation error %d.\n', c*err_cons);
        fprintf('relative true error: %d.\n', conv_rate/err_v);
        fprintf('relative calculated error: %d.\n', (err_calc+c*err_cons)/err_v);
        fprintf('efficiency index: %d.\n', (err_calc+c*err_cons)/conv_rate)
        fprintf('combined efficiency index: %d.\n', 3*(err_calc+c*err_cons)/(conv_rate + err_pot + c*err_cons))
    end
    
    %Plot convergence
    figure(3)
    loglog(grid_refinery, convergence_error)
    title('Convergence Plot');
    hold on
    y=grid_refinery.^(1.1);
    loglog(grid_refinery,y)
    loglog(grid_refinery, convergence_error, 'b.', 'Markersize', 20)
    hold off
    legend({'numerical scheme','10^{1.1}'}, 'Location','northwest')
    xlabel('gridsize')
    ylabel('error')
    xlim([0 0.07])
end


