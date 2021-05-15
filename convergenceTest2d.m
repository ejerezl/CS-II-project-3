%Diffusivity
k = @(x,y) 1+(x-1).*(y-1);
%k = @(x,y) 1;
% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);

%RHS function
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
%f = @(x,y) 1;

convergence_error = zeros(10,1);
grid_refinery = zeros(10,1);
n = zeros(10,1);
j = 1;
for i = 10:10:100
    fprintf('Loop with i = %d\n', i);
    %Determine number of cells in each direction
    nx = i;
    ny = i;
    dx = 1/(nx-1);
    dy = 1/(ny-1);
    
    
    grid_refinery(j) = dx;
    n(j) = nx;
    % Build matrices
    [A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);

    % Solve linear system
    u = A\b;
    U = reshape(u,nx,ny);
    % Calculate fabricated solution vector and compute error in "L2-norm"
    u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
    error = norm(u-u_fabricated_vect,2)*sqrt(dx*dy);
    error1 = norm(u-u_fabricated(cells(:,1),cells(:,2)),2)*dx;
    convergence_error(j) = error;
    j = j + 1;
end

y1 = (grid_refinery).^2;
loglog(grid_refinery, y1);
hold on
loglog(grid_refinery, convergence_error)
title('Convergence Plot');
legend({'x^2', 'Real error'}, 'Location','northwest')
xlabel('\Delta x')
ylabel('Error')
hold off
P = polyfit(log(grid_refinery), log(convergence_error),1);
slope = P(1);
fprintf("Convergence rate: %f\n", slope);
