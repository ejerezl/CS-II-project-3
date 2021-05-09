%Diffusivity
%k = @(x,y) 1+(x-1).*(y-1);
k = @(x,y) 1;
% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);

%RHS function
%f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
f = @(x,y) 1;
fileID = fopen('computations2D.txt', 'w');
for i = 4:7
    fprintf('Loop with i = %d\n', i);
    %Determine number of cells in each direction
    nx = 2^i;
    ny = 2^i;
    dx = 1/(nx-1);
    dy = 1/(ny-1);

    % Build matrices
    [A,b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);

    % Solve linear system
    u = A\b;

    % Calculate fabricated solution vector and compute error in "L2-norm"
    u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
    error = norm(u-u_fabricated_vect,2)*sqrt(dx*dy);


    % Compute error without knowing real solution
    [energy_error, conservation_integral] = energyError(U,-(K.*G)*u, 1, dx, dy, cells,nx, ny, 1);
    majorant = energy_error + (1/(sqrt(2)*pi))*conservation_integral;
    energy_norm = norm((k(0)^(1/2))*G*(u-u_fabricated_vect),2)*sqrt(dx*dy);
    relative_error = energy_error/norm((k(0)^(1/2))*G*u_fabricated_vect,2)*sqrt(dx*dy);
    
    ef_index_potential = majorant/energy_norm;
    
    fprintf(fileID, 'Computations with values: nx = %d, ny = %d\n', nx, ny);
    fprintf (fileID, 'Majorant: %f \nReal_error %f\nEnergy_norm %f\nRelative_error %f\nEfficiency index (potential) %f\n\n', majorant, error, energy_norm, relative_error, ef_index_potential);
end
fclose(fileID);