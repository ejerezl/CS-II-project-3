%Calculate error norm, k=1, f=1

function [err_true, err_calc, err_v] = energy_error_norm1(u, q, num_cells,edges, k)
dx = 1/(num_cells-1);

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
%u_fabricated = @(x) -0.5*x.^2 + 0.5*x;

err_true = 0;

%Calculate error for each cell and sum it up
for i=2:num_cells-1
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    ki = (k(grid(i)))^(-1);
    err1 = (0.5-a1)^2*dx;
    err2 = -(0.5-a1)*(grid(i+1)^2-grid(i)^2);
    err3 = 1/3*(grid(i+1)^3-grid(i)^3);
    err_true = err_true + ki*(err1 + err2+err3);
end

i=1;
a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
ki = (k(grid(i)))^(-1);
err1 = (0.5-a1)^2*dx/2;
err2 = -(0.5-a1)*(grid(i+1)^2-grid(i)^2);
err3 = 1/3*(grid(i+1)^3-grid(i)^3);
err_true = err_true + ki*(err1 + err2+err3);

i=num_cells;
a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
ki = (k(grid(i)))^(-1);
err1 = (0.5-a1)^2*dx/2;
err2 = -(0.5-a1)*(grid(i+1)^2-grid(i)^2);
err3 = 1/3*(grid(i+1)^3-grid(i)^3);
err_true = err_true + ki*(err1 + err2+err3);

err_true = sqrt(err_true*dx);

%Calculate energy error norm without true solution
err_calc = 0;

%Calculate error for each cell and sum it up

%%%%%%from the integral calculator :)%%%%%
for i=2:num_cells-1
    b0 = q(i,1);
    b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    ki = (k(grid(i)))^(-1);
    err1 = ((b0+a1*ki) + b1*dx)^3/(3*b1);
    err2 = -(b0+a1*ki)^3/(3*b1);
    err_calc = err_calc + ki*(err1 + err2);
end

%Left boundary (since here we only have dx/2)
i=1;
b0 = q(i,1);
b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
ki = (k(grid(i)))^(-1);
err1 = ((b0+a1*ki) + b1*(dx/2))^3/(3*b1);
err2 = -(b0+a1*ki)^3/(3*b1);
err_calc = err_calc + ki*(err1 + err2);

%Right boundary (since here we only have dx/2)
i=num_cells;
b0 = q(i,1);
b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
ki = (k(grid(i)))^(-1);
err1 = ((b0+a1*ki) + b1*(dx/2))^3/(3*b1);
err2 = -(b0+a1*ki)^3/(3*b1);
err_calc = err_calc + err1 + err2;

err_calc = sqrt(err_calc);


% Energy error integral of v
err_v = 0;
for i=1:num_cells
    a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
    ki = (k(grid(i)))^(-1);
    err1 = dx*(ki*a1)^2;
    err_v = err_v + err1;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some old calculations which we may want to
%%%%%%%%%%%%%%%%%%%%%%%%%%%% revisit at some point %%%%%%%%%%%%%%%%%%%%


%%%%%%My Calculations%%%%
% for i=2:num_cells-1
%     b0 = q(i,1);
%     b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
%     a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
%     ki = (k(grid(i)))^(-1);
%     err1 = (b0 + ki*a1)^2;
%     err2 = b1*dx^2/3;
%     err3 = dx*b1*(b0+ki*a1);
%     err_calc = err_calc + ki*(err1 + err2 + err3);
% end
% 
% i=1;
% b0 = q(i,1);
% b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
% a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
% ki = (k(grid(i)))^(-1);
% err1 = 0.5*(b0 + ki*a1)^2;
% err2 = b1*dx^2/(3*8);
% err3 = (dx/4)*b1*(b0+ki*a1);
% err_calc = err_calc + ki*(err1 + err2 + err3);
% 
% i=num_cells;
% b0 = q(i,1);
% b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
% a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
% ki = (k(grid(i)))^(-1);
% err1 = 0.5*(b0 + ki*a1)^2;
% err2 = b1*dx^2/(3*8);
% err3 = (dx/4)*b1*(b0+ki*a1);
% err_calc = err_calc + ki*(err1 + err2 + err3);
% 
% err_calc = sqrt(err_calc*dx);

%%%% Jans slides with correct sign%%%%%
% for i=1:num_cells
%     b0 = q(i,1);
%     b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
%     a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
%     ki = (k(grid(i)))^(-1);
%     err1 = (b0 + ki*a1)^2;
%     err2 = (b1^2)*dx/4;
%     err_calc = err_calc + ki*(err1 + err2);
% end
% err_calc = sqrt(err_calc*dx);





