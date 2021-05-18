%Diffusivity
k = @(x,y) 1+(x-1).*(y-1);
%k = @(x,y) 1;
% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);
q_fabricated = {@(x,y) -2*pi*(1+(x-1).*(y-1)).*cos(2*pi*x).*sin(2*pi*y), @(x,y) -2*pi*(1+(x-1).*(y-1)).*sin(2*pi*x).*cos(2*pi*y)};
gradu_1 = @(x,y) 2*pi*cos(2*pi*x).*sin(2*pi*y);
gradu_2 = @(x,y) 2*pi*sin(2*pi*x).*cos(2*pi*y);
%RHS function
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
%f = @(x,y) 1;

fileID = fopen('computations2D.txt', 'w');

majorant_values = zeros(5,1);
convergence_error = zeros(5,1);
grid_refinery = zeros(5,1);
star_norm_values = zeros(5,1);
energy_norm_values = zeros(5,1);
flag = 0;
j = 1;
for i = 20
    fprintf('Loop with i = %d\n', i);
    %Determine number of cells in each direction
    nx = i;
    ny = i;
    dx = 1/(nx-1);
    dy = 1/(ny-1);
    
    grid_refinery(j) = dx;
    [X,Y] = meshgrid(0:dy:1,0:dx:1);
    % Build matrices
    [A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);
    
    % Solve linear system
    v = A\b;
    V = reshape(v,nx,ny);
    
    % Calculate fabricated solution vector and compute error in "L2-norm"
    u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
    error = norm(u_fabricated_vect-v,2)*sqrt(dx*dy);
    convergence_error(j) = error;
    R = -(K.*G)*v;
    % Compute error without knowing real solution
    [energy_error, conservation_integral, energy_error_flux, energy_potential_norm, v_norm, energy_flux_norm] = energyError(V,R, k, dx, dy, cells,nx, ny, f, gradu_1, gradu_2);
    
    %Poincare constant
    eigenvalues = eigs(A, 10, "smallestabs");
    poincareConstant = 1/sqrt(eigenvalues(find(eigenvalues > 1.0001, 1)));
    
    
    majorant = energy_error + poincareConstant*conservation_integral;
    
    
    relative_error = energy_potential_norm/v_norm;
    relative_error_bound = majorant/v_norm;
    
    energy_norm_values(j) = energy_potential_norm;
    bound = majorant*poincareConstant;
    
    norm_star = energy_potential_norm + energy_error_flux + poincareConstant*conservation_integral;
    ef_index_combined = (3*majorant)/norm_star;
    ef_index_potential = majorant/energy_potential_norm;
    
    majorant_values(j) = majorant;
    star_norm_values(j) = norm_star;
    
    fprintf(fileID, 'Computations with values: nx = %d, ny = %d\n', nx, ny);
    fprintf (fileID, 'Majorant: %f \n Real_error_potential %f\nEnergy_potential_norm %f\n Norm_star %f \nRelative_error %f\nEfficiency index (potential) %f\n Efficiency index (combined) %f\nBound: %f\n\n', majorant, error, energy_potential_norm, norm_star, relative_error, ef_index_potential, ef_index_combined, bound);
    j = j + 1;
end
fclose(fileID);

if flag == 0
    figure(1);
    subplot(1,2,1);
    U = reshape(u_fabricated_vect,nx,ny);
    surf(X,Y, abs(U-V)./abs(U));
    axis([0.1 0.9 0.1 0.9 0 0.04])
    xlabel('x', 'FontSize', 18);
    ylabel('y', 'FontSize', 18);
    zlabel('Relative difference', 'FontSize', 18);
    title({'Relative difference between analytical','and numerical potential'},'FontSize',17)
    subplot(1,2,2);
    q_1 = q_fabricated{1}(cells(:,1),cells(:,2));
    q_2 = q_fabricated{2}(cells(:,1),cells(:,2));
    Q = [q_1, q_2];
    [r_1, r_2] = flux(G,K,v,cells,nx,ny);
    R = [r_1, r_2];
    norm_vector = zeros(400,1);
    norm_r = zeros(400,1);
    for i = 1:400
        norm_vector(i) = norm(q_1(i)-r_1(i),2);
        norm_r(i) = norm(q_1(i),2);
    end
    norma = reshape(norm_vector, nx, ny);
    surf(X,Y, norma);
    axis([0 1 0.25 0.75 0 0.1])
    xlabel('x', 'FontSize', 18);
    ylabel('y', 'FontSize', 18);
    zlabel('Relative difference', 'FontSize', 18);
    title({'Relative difference between analytical','and numerical flux'},'FontSize',17)

else
    figure(2)
    loglog(grid_refinery, 3*majorant_values, '-o');
    hold on
    loglog(grid_refinery, majorant_values, '-o');
    loglog(grid_refinery, star_norm_values, '-o');
    loglog(grid_refinery, convergence_error, '-o');
    loglog(grid_refinery, energy_norm_values, '-o');
    title('Error comparison', 'FontSize', 24);
    legend({'3Majorant', 'Majorant', 'Star norm', 'Real error', 'Energy norm'}, 'Location','northwest', 'FontSize', 16)
    xlabel('\Delta x', 'FontSize', 24)
    ylabel('Error', 'FontSize', 24)
    hold off
    P = polyfit(log(grid_refinery), log(convergence_error),1);
    slope = P(1);
    fprintf("Convergence rate: %f\n", slope);
end