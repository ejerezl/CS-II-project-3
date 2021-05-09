% Function that computes the error between the real solution and an
% approximation
% U = matrix with potentials
% F = flux
% k = from PDE (considering it a constant)
function [energy_error, conservation_integral, energy_error_flux] = energyError(U, F, k, dx, dy, cells, cells_x, cells_y, f)

energy_error = 0;
conservation_integral = 0;
energy_error_flux = 0;

%For computing the aproximation of the error we need to compute two
%integrals

%Computing the first term: ||k^(-1/2)(r+k*Deltav)||
%do a for iterating for each cell
%in each cell compute the cell potentials doing interpolations s.t. the
%vertices coincide with the lines of fluxes.
%Then compute a1,a2,a3,a4,b1,b2,b3,b4 doing interpolation for the general
%point (x,y) inside that cell.
%NOTE: It remains to compute the cells on the borders
for j=2:cells_y-1
    for i=2:cells_x-1
        
        %cell in which we are
        cell = cells_x*(j-1)+i;
        
        %Potentials in new cell(doing interpolation)
        A=(U(i,j)+U(i-1,j)+U(i,j-1)+U(i-1,j-1))/4;
        B=(U(i,j)+U(i+1,j)+U(i,j-1)+U(i+1,j-1))/4;
        C=(U(i,j)+U(i-1,j)+U(i,j+1)+U(i-1,j+1))/4;
        D=(U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1))/4;
        
        %Coordinates of new cell
        x1=cells(cell,1)-dx/2;
        x2=cells(cell,1)+dx/2;
        y1=cells(cell,2)-dy/2;
        y2=cells(cell,2)+dy/2;
        
        [a2, a3, a4, b1, b2, b3, b4] = interpolation(1,cells(cell,:),F,dx,dy,A,B,C,D,x1,y1);
        % this error correspond to the integral in that cell of r*r
        error1=(b1^2+b3^2)*dx*dy+b1*b2*(x2^2-x1^2)*dy+b3*b4*dx*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*dy+b4^2*dx*(y2^3-y1^3))/3;
        % This error correspond to the integral in that cell of delta v * delta v
        error2=(a2^2+a3^2)*dx*dy+a2*a4*dx*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*dy+a4^2*((x2^3-x1^3)*dy+dx*(y2^3-y1^3))/3;
        % This error correspond to the integral in that cell of r * delta v
        error3=(b1*a2+b3*a3)*dx*dy+(x2^2-x1^2)*dy*(b2*a2+b3*a4)/2+dx*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

        % Add error at each cell plus the previously computed
        energy_error=energy_error+error1/k+k*error2+2*error3;
        
        %Computing the conservation integral when f = 1
        conservation_integral = conservation_integral + (f-b2-b4)^2*(x2-x1)*(y2-y1);

        %Computing energy error of the flux er
        energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;
    end
end

% Bottom cells (Execpt corners)
for i = 2:cells_x-1
    j = 1;
    cell = cells_x*(j-1)+i;
    A = (U(i,j)+U(i-1, j))/2;
    B = (U(i,j) + U(i+1,j))/2;
    C = (U(i,j)+U(i-1,j)+U(i,j+1)+U(i-1,j+1))/4;
    D = (U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1))/4;
    
    %Coordinates of new cell
    x1=cells(cell,1)-dx/2;
    x2=cells(cell,1)+dx/2;
    y1=cells(cell,2);
    y2=cells(cell,2)+dy/2;
    
    [a2, a3, a4, ~,~,~,~] = interpolation(0, cells(cell,:),F,dx,dy/2,A,B,C,D,x1,y1);
    Fr = F(cells(cell, 4));
    Fl = F(cells(cell,6));
    Ft = F(cells(cell,5));
    Fb = Ft - f*(y2-y1);
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
    
    % this error correspond to the integral in that cell of r*r
    error1=(b1^2+b3^2)*dx*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*dx*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*dx*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of delta v * delta v
    error2=(a2^2+a3^2)*dx*(dy/2)+a2*a4*dx*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+dx*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of r * delta v
    error3=(b1*a2+b3*a3)*dx*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+dx*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

    % Add error at each cell plus the previously compute
    energy_error=energy_error+error1/k+k*error2+2*error3;
    %Computing the conservation integral when f = 1
    conservation_integral = conservation_integral + (f-b2-b4)^2*(x2-x1)*(y2-y1);
    
    %Computing energy error of the flux er
    energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;
end

% Right cells (Except corners)
for j = 2:cells_y - 1
    i = cells_x;
    cell = cells_x*(j-1)+i;
    A = (U(i,j)+U(i-1,j)+U(i,j-1)+U(i-1,j-1))/4;
    B = (U(i,j) + U(i,j-1))/2;
    C = (U(i,j)+U(i-1,j)+U(i,j+1)+U(i-1,j+1))/4;
    D = (U(i,j) + U(i, j+1))/2;
    
    %Coordinates of new cell
    x1=cells(cell,1)-dx/2;
    x2=cells(cell,1);
    y1=cells(cell,2)-dy/2;
    y2=cells(cell,2)+dy/2;
    
    [a2, a3, a4,~,~,~,~] = interpolation(0, cells(cell,:),F,dx/2,dy,A,B,C,D,x1,y1);
    
    Fl = F(cells(cell, 6));
    Fr = Fl + f*(x2-x1);
    Ft = F(cells(cell, 5));
    Fb = F(cells(cell,3));
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
    
    % this error correspond to the integral in that cell of r*r
    error1=(b1^2+b3^2)*(dx/2)*dy+b1*b2*(x2^2-x1^2)*dy+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*dy+b4^2*(dx/2)*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of delta v * delta v
    error2=(a2^2+a3^2)*(dx/2)*dy+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*dy+a4^2*((x2^3-x1^3)*dy+(dx/2)*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of r * delta v
    error3=(b1*a2+b3*a3)*(dx/2)*dy+(x2^2-x1^2)*dy*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

    % Add error at each cell plus the previously compute
    energy_error=energy_error+error1/k+k*error2+2*error3;
    
    %Computing the conservation integral when f = 1
    conservation_integral = conservation_integral + (f-b2-b4)^2*(x2-x1)*(y2-y1);
    
    %Computing energy error of the flux er
    energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;
end

%Top cells (Except corners)
for i = 2:cells_x-1
    j = cells_y;
    cell = cells_x*(j-1)+i;
    A= (U(i,j)+U(i-1,j)+U(i,j-1)+U(i-1,j-1))/4;
    B= (U(i,j)+U(i+1,j)+U(i,j-1)+U(i+1,j-1))/4;
    C = (U(i,j)+U(i-1,j))/2;
    D = (U(i,j)+U(i+1,j))/2;
    
    %Coordinates of new cell
    x1=cells(cell,1)-dx/2;
    x2=cells(cell,1)+dx/2;
    y1=cells(cell,2)-dy/2;
    y2=cells(cell,2);
    
    [a2, a3, a4, ~,~,~,~] = interpolation(0, cells(cell,:),F,dx,dy/2,A,B,C,D,x1,y1);
    
    Fr = F(cells(cell, 4));
    Fl = F(cells(cell,6));
    Fb = F(cells(cell,3));
    Ft = Fb + f*(y2-y1);
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
    
    % this error correspond to the integral in that cell of r*r
    error1=(b1^2+b3^2)*dx*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*dx*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*dx*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of delta v * delta v
    error2=(a2^2+a3^2)*dx*(dy/2)+a2*a4*dx*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+dx*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of r * delta v
    error3=(b1*a2+b3*a3)*dx*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+dx*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

    % Add error at each cell plus the previously compute
    energy_error=energy_error+error1/k+k*error2+2*error3;
    
    %Computing the conservation integral when f = 1
    conservation_integral = conservation_integral + (f-b2-b4)^2*(x2-x1)*(y2-y1);
    
    %Computing energy error of the flux er
    energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;
end

%Left cells (Except corners)
for j = 2:cells_y - 1
    i = 1;
    cell = cells_x*(j-1)+i;
    
    A = (U(i,j)+U(i, j-1))/2;
    B = (U(i,j)+U(i+1,j)+U(i,j-1)+U(i+1,j-1))/4;
    C = (U(i,j)+U(i,j+1))/2;
    D = (U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1))/4;
    
    %Coordinates of new cell
    x1=cells(cell,1);
    x2=cells(cell,1)+dx/2;
    y1=cells(cell,2)-dy/2;
    y2=cells(cell,2)+dy/2;
    
    [a2, a3, a4,~,~,~,~] = interpolation(0, cells(cell,:),F,dx/2,dy,A,B,C,D,x1,y1);
    
    Fr = F(cells(cell, 4));
    Fl = Fr - f*(x2-x1);
    Ft = F(cells(cell, 5));
    Fb = F(cells(cell,3));
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
    
    % this error correspond to the integral in that cell of r*r
    error1=(b1^2+b3^2)*(dx/2)*dy+b1*b2*(x2^2-x1^2)*dy+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*dy+b4^2*(dx/2)*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of delta v * delta v
    error2=(a2^2+a3^2)*(dx/2)*dy+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*dy+a4^2*((x2^3-x1^3)*dy+(dx/2)*(y2^3-y1^3))/3;
    % This error correspond to the integral in that cell of r * delta v
    error3=(b1*a2+b3*a3)*(dx/2)*dy+(x2^2-x1^2)*dy*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

    % Add error at each cell plus the previously compute
    energy_error=energy_error+error1/k+k*error2+2*error3;
    
    %Computing the conservation integral when f = 1
    conservation_integral = conservation_integral + (f-b2-b4)^2*(x2-x1)*(y2-y1);
    
    %Computing energy error of the flux er
    energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;
end

%Corners
%Left bottom corner
cell = 1;
i=1;
j=1;
A = U(i,j);
B = (U(i,j)+U(i+1,j))/2;
C = (U(i,j)+U(i,j+1))/2;
D = (U(i,j)+U(i+1,j)+U(i,j+1)+U(i+1,j+1))/4;
    
%Coordinates of new cell
x1=cells(cell,1);
x2=cells(cell,1)+dx/2;
y1=cells(cell,2);
y2=cells(cell,2)+dy/2;

[a2,a3,a4,~,~,~,~] = interpolation(0, cells(cell,:), F, dx/2,dy/2, A,B,C,D,x1,y1);

b4 = (x2*f)/(x2+y2);
b2 = f -b4;
b3 = -y2*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
error1=(b1^2+b3^2)*(dx/2)*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of delta v * delta v
error2=(a2^2+a3^2)*(dx/2)*(dy/2)+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of r * delta v
error3=(b1*a2+b3*a3)*(dx/2)*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1/k+k*error2+2*error3;
% conservation_integral = conservation_integral +
% (f-b2-b4)^2*(x2-x1)*(y2-y1); this will be zero cause one condition is f =
% b2+b4
%Computing energy error of the flux er
energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;

%right bottom corner
i = cells_x;
j = 1;
cell = cells_x*(j-1)+i;

A = (U(i-1,j)+U(i, j))/2;
B = U(i,j);
C = (U(i,j)+U(i-1,j)+U(i,j+1)+U(i-1,j+1))/4;
D = (U(i,j)+U(i,j+1))/2;
    
%Coordinates of new cell
x1=cells(cell,1)-dx/2;
x2=cells(cell,1);
y1=cells(cell,2);
y2=cells(cell,2)+dy/2;

[a2,a3,a4,~,~,~,~] = interpolation(0, cells(cell,:), F, dx/2,dy/2, A,B,C,D,x1,y1);

b4 = (x1*f)/(x1+y2);
b2 = f -b4;
b3 = -y2*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
error1=(b1^2+b3^2)*(dx/2)*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of delta v * delta v
error2=(a2^2+a3^2)*(dx/2)*(dy/2)+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of r * delta v
error3=(b1*a2+b3*a3)*(dx/2)*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1/k+k*error2+2*error3;
% conservation_integral = conservation_integral +
% (f-b2-b4)^2*(x2-x1)*(y2-y1); this will be zero cause one condition is f =
% b2+b4
%Computing energy error of the flux er
energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;

%Left upper corner
i = 1;
j = cells_y;
cell = cells_x*(j-1)+i;

A = (U(i,j)+U(i, j-1))/2;
B = (U(i,j)+U(i+1,j)+U(i,j-1)+U(i+1,j-1))/4;
C = U(i,j);
D = (U(i,j)+U(i+1,j))/2;
    
%Coordinates of new cell
x1=cells(cell,1);
x2=cells(cell,1)+dx/2;
y1=cells(cell,2)-dy/2;
y2=cells(cell,2);

[a2,a3,a4,~,~,~,~] = interpolation(0, cells(cell,:), F, dx/2,dy/2, A,B,C,D,x1,y1);

b4 = (f*x2)/(x2+y1);
b2 = f -b4;
b3 = -y1*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
error1=(b1^2+b3^2)*(dx/2)*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of delta v * delta v
error2=(a2^2+a3^2)*(dx/2)*(dy/2)+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of r * delta v
error3=(b1*a2+b3*a3)*(dx/2)*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1/k+k*error2+2*error3;
    
%Computing energy error of the flux er
energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;

%Computing the conservation integral when f = 1
% conservation_integral = conservation_integral +
% (f-b2-b4)^2*(x2-x1)*(y2-y1); this will be zero cause one condition is f =
% b2+b4


%Right upper corner
i = cells_x;
j = cells_y;
cell = cells_x*(j-1)+i;

A = (U(i,j)+U(i-1,j)+U(i,j-1)+U(i-1,j-1))/4;
B = (U(i,j)+U(i,j-1))/2;
C = (U(i,j)+U(i-1,j))/2;
D = U(i,j);
    
%Coordinates of new cell
x1=cells(cell,1)-dx/2;
x2=cells(cell,1);
y1=cells(cell,2)-dy/2;
y2=cells(cell,2);

[a2,a3,a4,~,~,~,~] = interpolation(0, cells(cell,:), F, dx/2,dy/2, A,B,C,D,x1,y1);

b4 = (f*x1)/(x1+y1);
b2 = f -b4;
b3 = -y1*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
error1=(b1^2+b3^2)*(dx/2)*(dy/2)+b1*b2*(x2^2-x1^2)*(dy/2)+b3*b4*(dx/2)*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*(dy/2)+b4^2*(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of delta v * delta v
error2=(a2^2+a3^2)*(dx/2)*(dy/2)+a2*a4*(dx/2)*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*(dy/2)+a4^2*((x2^3-x1^3)*(dy/2)+(dx/2)*(y2^3-y1^3))/3;
% This error correspond to the integral in that cell of r * delta v
error3=(b1*a2+b3*a3)*(dx/2)*(dy/2)+(x2^2-x1^2)*(dy/2)*(b2*a2+b3*a4)/2+(dx/2)*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1/k+k*error2+2*error3;
% conservation_integral = conservation_integral +
% (f-b2-b4)^2*(x2-x1)*(y2-y1); this will be zero cause one condition is f =
% b2+b4

%Computing energy error of the flux er
energy_error_flux = energy_error_flux + k*error2 + (1/k)*error1 + 2*error3;

conservation_integral = sqrt(conservation_integral);
energy_error = sqrt(energy_error);
energy_error_flux = sqrt(energy_error_flux);

end