%Diffusivity
k = @(x,y) 1+(x-1).*(y-1);
%k = @(x,y) 1;
% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);

gradu_1 = @(x,y) 2*pi*cos(2*pi*x).*sin(2*pi*y);
gradu_2 = @(x,y) 2*pi*sin(2*pi*x).*cos(2*pi*y);
%RHS function
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
%f = @(x,y) 1;

fileID = fopen('computations2D.txt', 'w');

majorant_values = zeros(5,1);
grid_refinery = zeros(5,1);
%j = 1;
for i = 4:7
    fprintf('Loop with i = %d\n', 2^i);
    %Determine number of cells in each direction
    nx = 2^i;
    ny = 2^i;
    dx = 1/(nx-1);
    dy = 1/(ny-1);
    
    %grid_refinery(j) = dx;
    
    % Build matrices
    [A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);
    
    % Solve linear system
    v = A\b;
    V = reshape(v,nx,ny);
    
    % Calculate fabricated solution vector and compute error in "L2-norm"
    u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
    error = norm(u_fabricated_vect-v,2)*sqrt(dx*dy);
    
    
    % Compute error without knowing real solution
    [energy_error, conservation_integral, energy_error_flux, energy_potential_norm, v_norm, energy_flux_norm] = energyError(V,-(K.*G)*v, k, dx, dy, cells,nx, ny, f, gradu_1, gradu_2);
    
    %Poincare constant
    eigenvalues = eigs(A, 10, "smallestabs");
    poincareConstant = 1/sqrt(eigenvalues(find(eigenvalues > 1.0001, 1)));
    
    
    majorant = energy_error + poincareConstant*conservation_integral;
    
    
    relative_error = energy_potential_norm/v_norm
    relative_error_bound = majorant/v_norm
    
    
    bound = majorant*poincareConstant;
    
    norm_star = energy_potential_norm + energy_error_flux + poincareConstant*conservation_integral;
    ef_index_combined = (3*majorant)/norm_star
    ef_index_potential = majorant/energy_potential_norm
    
    %majorant_values(j) = majorant;
    
    fprintf(fileID, 'Computations with values: nx = %d, ny = %d\n', nx, ny);
    fprintf (fileID, 'Majorant: %f \n Real_error_potential %f\nEnergy_potential_norm %f\n Norm_star %f \nRelative_error %f\nEfficiency index (potential) %f\n Efficiency index (combined) %f\nBound: %f\n\n', majorant, error, energy_potential_norm, norm_star, relative_error, ef_index_potential, ef_index_combined, bound);
    %j = j + 1;
end

fclose(fileID);
% loglog(grid_refinery, majorant_values)
% title('Majorant Plot');
% P = polyfit(log(grid_refinery), log(majorant_values),1);
% slope = P(1);
% fprintf("Convergence rate: %f\n", slope);