%Calculate error norm, k=1, f = @(x) 4*pi*pi*sin(2*pi*x)
%err_true = true error of potential |k^1/2 * (u-v)|
%err_calc = error estimate, |k^^-1/2 * (r- dxv)| (dxv = derivative of v)
%err_v = error of approx. v, |k^1/2 * dxv|
%err_pot = true error of flux,|k^^-1/2 * (q-r)|

function [err_true, err_calc, err_cons, err_v, err_flux] = energy_error_norm2(u, q, num_cells,edges, k)
%dx = 1/(num_cells-1);

%define new grid of edges which includes 0 and 1
grid = zeros(num_cells+1,1);
grid(2:(num_cells)) = edges(:,1,1);
grid(end) = 1;

%define potential on grid points as averages of the adjacents cells
pot = zeros(size(grid));
pot(1) = 0; %Boundary condition
for i=2:size(grid)-1
    pot(i) = 0.5*(u(i)+u(i-1));
end
pot(end) = 0; %Boundary condition

%Calculate energy error norm with true solution
%u_fabricated = @(x) @(x) sin(2*pi*x);

err_true = 0;

%Calculate error for each cell and sum it up
for i=1:num_cells
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    err1 = (pi*cos(2*pi*grid(i+1)) - 2*a1)*sin(2*pi*grid(i+1)) - (pi*cos(2*pi*grid(i)) - 2*a1)*sin(2*pi*grid(i));
    err2 = (a1^2 + 2*pi^2)*(grid(i+1)-grid(i));
    err_true = err_true + err1 + err2;
end

err_true = sqrt(err_true);


%Calculate energy error norm without true solution
err_calc = 0;

%Calculate error for each cell and sum it up

%%%%%%from the integral calculator :)%%%%%
for i=1:num_cells
    b0 = q(i,1);
    b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    ki = (k(grid(i)));
    err1 = ((b0+a1*ki) + b1*(grid(i+1)-grid(i)))^3/(3*b1);
    err2 = -(b0+a1*ki)^3/(3*b1);
    err_calc = err_calc + ki^(-1)*(err1 + err2);
end

err_calc = sqrt(err_calc);

% Calculate conservation integral
err_cons = 0;
for i=1:num_cells
    b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
    err1 = 8*pi^4*(grid(i+1)-grid(i));
    err2 = -2*pi^3*(sin(4*pi*grid(i+1))-sin(4*pi*grid(i)));
    err3 = 4*pi*b1*(cos(2*pi*grid(i+1))-cos(2*pi*grid(i)));
    err4 = b1^2*(grid(i+1)-grid(i));
    err_cons = err_cons + err1 + err2 + err3 + err4;
end

err_cons = sqrt(err_cons);


%Calculatate true error of potential
err_flux = 0;

for i=1:num_cells
    b0 = q(i,1);
    b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
    err1 = (pi*cos(2*pi*grid(i+1))+2*(b1*(grid(i+1)-grid(i))+b0))*sin(2*pi*grid(i+1))+(b1*cos(2*pi*grid(i+1)))/pi+(b1^2*grid(i+1)^3)/3-b1*(b1*grid(i)-b0)*grid(i+1)^2+((b1*grid(i)-b0)^2+2*pi^2)*grid(i+1);
    err2 = -((pi*cos(2*pi*grid(i))+2*(b1*(grid(i)-grid(i))+b0))*sin(2*pi*grid(i))+(b1*cos(2*pi*grid(i)))/pi+(b1^2*grid(i)^3)/3-b1*(b1*grid(i)-b0)*grid(i)^2+((b1*grid(i)-b0)^2+2*pi^2)*grid(i));
    err_flux = err_flux + err1 +err2;
end
err_flux = sqrt(err_flux);


% Energy error integral of v
err_v = 0;
for i=1:num_cells
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    ki = (k(grid(i)));
    err1 = (grid(i+1)-grid(i))*(ki^(1/2)*a1)^2;
    err_v = err_v + err1;
end
err_v = sqrt(err_v);
