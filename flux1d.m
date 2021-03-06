%Calculate flux at cells

function [q] = flux1d(k,u,dx,num_cells, k_const, f_const)
q = zeros(num_cells+1,1);
    
for i=2:num_cells
    q(i) = -k(0.5*dx+dx*(i-1))/dx * (u(i)-u(i-1));
end
     
%Calculation of boundary flux: depends on f and k

if k_const && f_const
    %left boundary
    q(1) = q(2) - dx*0.5;
    %right boundary
    q(end) = q(end-1) + dx*0.5; 
elseif k_const
    %left boundary, x_1/2 = dx/2
    q(1) = q(2) + k(dx/2)*pi*(cos(2*pi*dx/2)-1);
    %right boundary, x_(N-1/2) = 1-dx/2
    q(end) = q(end-1) - k(1-dx/2)*pi*(1-cos(2*pi*(1-dx/2)));
else
    %left boundary, x_1/2 = dx/2
    q(1) = q(2) + k(dx/2)*2*pi*cos(2*pi*dx/2) - 2*pi*k(0)*cos(0);
    %right boundary, x_(N-1/2) = 1-dx/2
    q(end) = q(end-1) - 2*pi*k(1)*cos(2*pi*1) + 2*pi*k(1-dx/2)*cos(2*pi*(1-dx/2));
end

    

     