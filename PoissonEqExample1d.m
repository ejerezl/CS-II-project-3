
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

%Plot difference of solutions
figure(3)
X = (0:dx:1);
subplot(1,2,1)
plot(X,abs(u-U_fabricated))
xlabel('x', 'FontSize', 16)
ylabel('u_{true} - u_{analytical}','FontSize',16)
title({'absolute difference between analytical','and numerical potential'},'FontSize',14)

X = zeros(nx+1,1); %new grid
X(1) = 0;
X(end) = 1;
for i=2:nx
    X(i) = 0.5*dx + dx*(i-2);
end
subplot(1,2,2)
plot(X,abs(Q_fabricated-q))
xlabel('x', 'FontSize', 16)
ylabel('q_{true} - q_{analytical}','FontSize',16)
title({'absolute difference between analytical','and numerical flux'},'FontSize',14)

%%%%% Error calculation & Convergence %%%%%

if k_const && f_const
    %Energy error
    [err_true,err_calc, err_cons, err_v, err_flux] = energy_error_norm1(u, q, nx, edges, k1);

    % Runge Kutta type 1 error estimate

    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1(nx, 1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);

    %Convergence rate
    c = 1/pi; %Poincare constant for k=1
    convergence_error = zeros(6,1);
    grid_refinery = zeros(6,1);
    majorant = zeros(6,1);
    combined_error = zeros(6,1);
    efficiency = zeros(6,1);
    efficiency_star = zeros(6,1);
    for i=1:6
        [conv_rate, grid_size, err_calc, err_cons, err_v, err_flux] = convergence1(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        majorant(i) =  err_calc+c*err_cons;
        combined_error(i) = (conv_rate + err_flux + c*err_cons);
        efficiency(i) = (err_calc + c*err_cons)/conv_rate;
        efficiency_star(i) = 3*(err_calc + c*err_cons)/(conv_rate + err_flux + c*err_cons);        
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('energy error: %d.\n', err_calc);
        fprintf('conservation error: %d.\n', c*err_cons);
        fprintf('error of flux: %d.\n', err_flux);
        fprintf('Majorant: %d.\n', err_calc+c*err_cons);
        fprintf('3*Majorant: %d.\n', 3*(err_calc+c*err_cons));
        fprintf('efficiency index: %d.\n', (err_calc + c*err_cons)/conv_rate)
        fprintf('combined efficiency index: %d.\n', 3*(err_calc + c*err_cons)/(conv_rate + err_flux + c*err_cons))
    end
   
    
    %Plot convergence
    figure(4)
    yyaxis left
    loglog(grid_refinery, convergence_error,'-o')
    title('Convergence Plot');
    hold on
    loglog(grid_refinery,majorant, ':o')
    loglog(grid_refinery, 3*majorant, '-.o')
    loglog(grid_refinery, combined_error, '--o')
    ylabel('error', 'FontSize', 16)
    
    yyaxis right
    plot(grid_refinery, efficiency, '--o')
    plot(grid_refinery, efficiency_star, '-.o')
    yline(1)
    yline(3)
    hold off
    
    legend({'$\| k^(1/2)\nabla (u-v)\|$', '$\mathcal{M}(r,v,f)$', '$3\mathcal{M}(r,v,f)$', '$\| (e_v, e_r)\|_\ast$', '$I_v$', '$I_v^\ast$' }, 'Location','northwest', 'Interpreter', 'latex')
    xlabel('gridsize', 'FontSize', 16)
    ylabel('efficiency index', 'FontSize', 16)
    
    xlim([0 0.07])
    ylim([0.8 3.5])

elseif k_const 
    c = 1/pi; %Poincare constant for k=1
    %Energy error
    [err_true,err_calc, err_cons, err_v, err_flux] = energy_error_norm2(u, q, nx, edges, k1);
    
    %Runge Kutta estimate type 1
    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1_2(nx, 1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);
    
    %Convergence rate
    convergence_error = zeros(7,1);
    grid_refinery = zeros(7,1);
    majorant = zeros(7,1);
    combined_error = zeros(7,1);
    efficiency = zeros(7,1);
    efficiency_star = zeros(7,1);
    for i=1:7
        [conv_rate, grid_size, err_calc, err_cons, err_v, err_flux] = convergence2(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        majorant(i) = (err_calc+c*err_cons);
        combined_error(i) = (conv_rate + err_flux + c*err_cons);
        efficiency(i) = (err_calc + c*err_cons)/conv_rate;
        efficiency_star(i) = 3*(err_calc + c*err_cons)/(conv_rate + err_flux + c*err_cons); 
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('energy error: %d.\n', err_calc);
        fprintf('error of flux: %d.\n', err_flux);
        fprintf('conservation error: %d.\n', c*err_cons);
        fprintf('combined error norm: %d.\n', conv_rate + err_flux + c*err_cons);
        fprintf('Majorant: %d.\n', err_calc+c*err_cons);
        fprintf('3*Majorant: %d.\n', 3*(err_calc+c*err_cons));
        fprintf('efficiency index: %d.\n', (err_calc+c*err_cons)/conv_rate)
        fprintf('combined efficiency index: %d.\n', 3*(err_calc+c*err_cons)/(conv_rate + err_flux + c*err_cons))
    end
    
    
    %Plot convergence
    figure(4)
    yyaxis left
    loglog(grid_refinery, convergence_error, '-o')
    title({'Convergence Plot 1D'}, 'FontSize', 14);
    hold on
    loglog(grid_refinery, majorant, ':o')
    loglog(grid_refinery, 3*majorant, '-.o')
    loglog(grid_refinery, combined_error, '--o')
    ylabel({'error'}, 'FontSize', 16)

    yyaxis right
    plot(grid_refinery, efficiency, '--o')
    plot(grid_refinery, efficiency_star, '-.o')
    yline(1)
    yline(3)
    legend({'$\| k^(1/2)\nabla (u-v)\|$', '$\mathcal{M}(r,v,f)$', '$3\mathcal{M}(r,v,f)$', '$\| (e_v, e_r)\|_\ast$', '$I_v$','$I^\ast_v$'}, 'Location','northwest', 'Interpreter', 'latex')
    hold off
    
    xlabel({'gridsize'}, 'FontSize', 16)
    ylabel({'efficiency index'},  'FontSize', 16)
    xlim([0 0.07])
    ylim([0.8 3.5])
   
else
    c=1/sqrt(14.3249); %14.3249 is smallest eigenvalue of A which is not 1
    [err_true, err_calc, err_cons, err_v, err_flux] = energy_error_norm3(u, q, nx, edges, k1);
    
    %Runge Kutta estimate type 1
    %2nd entry: Order of convergence, 3rd entry: ratio between h1 and h2
    [err1, err2] = rktype1_3(nx, 1.1, 8); 
    fprintf('True error: %d.\n', err1);
    fprintf('Runge-Kutta estimate: %d.\n', err2);
    
    %Convergence rate
    convergence_error = zeros(6,1);
    grid_refinery = zeros(6,1);
    majorant = zeros(6,1);
    combined_error = zeros(6,1);
    efficiency = zeros(6,1);
    efficiency_star = zeros(6,1);
    for i=1:6
        [conv_rate, grid_size, err_calc, err_cons, err_v, err_flux] = convergence3(nx);
        err1 = [conv_rate, grid_size];
        convergence_error(i) = err1(1);
        grid_refinery(i) = err1(2);
        majorant(i) = (err_calc+c*err_cons);
        combined_error(i) = (conv_rate + err_flux + c*err_cons); 
        efficiency(i) = (err_calc + c*err_cons)/conv_rate;
        efficiency_star(i) = 3*(err_calc + c*err_cons)/(conv_rate + err_flux + c*err_cons); 
        nx = 2*nx;
        disp(i)
        fprintf('true error: %d.\n', conv_rate);
        fprintf('energy error: %d.\n', err_calc);
        fprintf('error of flux: %d.\n', err_flux);
        fprintf('conservation error: %d.\n', c*err_cons);
        fprintf('combined error norm: %d.\n', conv_rate + err_flux + c*err_cons);
        fprintf('Majorant: %d.\n', err_calc+c*err_cons);
        fprintf('3*Majorant: %d.\n', 3*(err_calc+c*err_cons));
        fprintf('efficiency index: %d.\n', (err_calc+c*err_cons)/conv_rate)
        fprintf('combined efficiency index: %d.\n', 3*(err_calc+c*err_cons)/(conv_rate + err_flux + c*err_cons))
    end
    %Plot convergence
    figure(4)
    yyaxis left
    loglog(grid_refinery, convergence_error, '-o')
    title({'Convergence Plot 1D'}, 'FontSize', 16);
    hold on
    loglog(grid_refinery, majorant,':o')
    loglog(grid_refinery, 3*majorant, '--o')
    loglog(grid_refinery, combined_error, '-.o')
    ylabel({'error'}, 'FontSize', 16)

    yyaxis right
    plot(grid_refinery, efficiency, '--o')
    plot(grid_refinery, efficiency_star, '-o')
    yline(1)
    yline(3)
    hold off
    
    legend({'$\| k^{1/2}\nabla (u-v)\|$', '$\mathcal{M}(r,v,f)$', '$3\mathcal{M}(r,v,f)$', '$\| (e_v, e_r)\|_\ast$', '$I_v$','$I^\ast_v$'}, 'Location','northwest', 'Interpreter', 'latex')
    
    xlabel({'gridsize'}, 'FontSize', 16)
    ylabel({'efficiency index'},  'FontSize', 16)
    xlim([0 0.07])
    ylim([0.8 3.5])
end



