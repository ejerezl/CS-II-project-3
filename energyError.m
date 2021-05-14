% Function that computes the error between the real solution and an
% approximation
% U = matrix with potentials
% F = flux
% k = from PDE (considering it a constant)
function [energy_error, conservation_integral, energy_error_flux, energy_norm, v_norm, energy_q_norm] = energyError(U, F, k, dx, dy, cells, cells_x, cells_y, f, u_1, u_2)

energy_error = 0;
conservation_integral = 0;
energy_error_flux = 0;
energy_norm = 0;
v_norm = 0;
energy_q_norm = 0;
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
        fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
        error1=integral2(fun,x1,x2,y1,y2);
        % This error correspond to the integral in that cell of delta v * delta v
        fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
        error2=integral2(fun,x1,x2,y1,y2);
        % This error correspond to the integral in that cell of r * delta v
        fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
        error3=integral2(fun,x1,x2,y1,y2);
        v_norm = v_norm + error3;
        
        % Add error at each cell plus the previously computed
        energy_error=energy_error+error1+error2+error3;
        
        %Computing the conservation integral
        fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
        conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
        q_norm = @(x,y) k(x,y).*fun(x,y);
        energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);
        
        %Computing energy error of the flux er
        fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
        error2_f = integral2(fun,x1,x2,y1,y2);
        fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
        error3_f = integral2(fun,x1,x2,y1,y2);
        energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;
        
        %Computing energy norm
        fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
        error_norm = integral2(fun,x1,x2,y1,y2);
        energy_norm = energy_norm + error2_f - error_norm + error3;
        
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
    Fb = Ft - f(x1,y1)*(y2-y1);
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
 
    % this error correspond to the integral in that cell of r*r
    fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
    error1=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of delta v * delta v
    fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
    error2=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of r * delta v
    fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
    error3=integral2(fun,x1,x2,y1,y2);
    v_norm = v_norm + error3;
    
    % Add error at each cell plus the previously computed
    energy_error=energy_error+error1+error2+error3;
        
    %Computing the conservation integral
    fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
    conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
    q_norm = @(x,y) k(x,y).*fun(x,y);
    energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);
    
    %Computing energy error of the flux er
    fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
    error2_f = integral2(fun,x1,x2,y1,y2);
    fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
    error3_f = integral2(fun,x1,x2,y1,y2);
    energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

    %Computing energy norm
    fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
    error_norm = integral2(fun,x1,x2,y1,y2);
    energy_norm = energy_norm + error2_f - error_norm + error3;
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
    Fr = Fl + f(x1,y1)*(x2-x1);
    Ft = F(cells(cell, 5));
    Fb = F(cells(cell,3));
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);
    
    % this error correspond to the integral in that cell of r*r
    fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
    error1=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of delta v * delta v
    fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
    error2=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of r * delta v
    fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
    error3=integral2(fun,x1,x2,y1,y2);
    v_norm = v_norm + error3;
    
    % Add error at each cell plus the previously computed
    energy_error=energy_error+error1+error2+error3;
        
    %Computing the conservation integral
    fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
    conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
    q_norm = @(x,y) k(x,y).*fun(x,y);
    energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);
    
    %Computing energy error of the flux er
    fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
    error2_f = integral2(fun,x1,x2,y1,y2);
    fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
    error3_f = integral2(fun,x1,x2,y1,y2);
    energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

    %Computing energy norm
    fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
    error_norm = integral2(fun,x1,x2,y1,y2);
    energy_norm = energy_norm + error2_f - error_norm + error3;
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
    Ft = Fb + f(x1,y1)*(y2-y1);
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);

    % this error correspond to the integral in that cell of r*r
    fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
    error1=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of delta v * delta v
    fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
    error2=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of r * delta v
    fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
    error3=integral2(fun,x1,x2,y1,y2);
    v_norm = v_norm + error3;
    
    % Add error at each cell plus the previously computed
    energy_error=energy_error+error1+error2+error3;
        
    %Computing the conservation integral
    fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
    conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
    q_norm = @(x,y) k(x,y).*fun(x,y);
    energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);
    
    %Computing energy error of the flux er
    fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
    error2_f = integral2(fun,x1,x2,y1,y2);
    fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
    error3_f = integral2(fun,x1,x2,y1,y2);
    energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

    %Computing energy norm
    fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
    error_norm = integral2(fun,x1,x2,y1,y2);
    energy_norm = energy_norm + error2_f - error_norm + error3;
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
    Fl = Fr - f(x1,y1)*(x2-x1);
    Ft = F(cells(cell, 5));
    Fb = F(cells(cell,3));
    b1 = Fl + x1*(Fl-Fr)/(x2-x1);
    b2 = (Fr-Fl)/(x2-x1);
    b3 = Fb + y1*(Fb-Ft)/(y2-y1);
    b4 = (Ft-Fb)/(y2-y1);

    % this error correspond to the integral in that cell of r*r
    fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
    error1=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of delta v * delta v
    fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
    error2=integral2(fun,x1,x2,y1,y2);
    % This error correspond to the integral in that cell of r * delta v
    fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
    error3=integral2(fun,x1,x2,y1,y2);
    v_norm = v_norm + error3;
    
    % Add error at each cell plus the previously computed
    energy_error=energy_error+error1+error2+error3;
        
    %Computing the conservation integral
    fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
    conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
    q_norm = @(x,y) k(x,y).*fun(x,y);
    energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);
    
    %Computing energy error of the flux er
    fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
    error2_f = integral2(fun,x1,x2,y1,y2);
    fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
    error3_f = integral2(fun,x1,x2,y1,y2);
    energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

    %Computing energy norm
    fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
    error_norm = integral2(fun,x1,x2,y1,y2);
    energy_norm = energy_norm + error2_f - error_norm + error3;
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

b4 = (x2*f(x1,y1))/(x2+y2);
b2 = f(x1,y1) -b4;
b3 = -y2*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
error1=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of delta v * delta v
fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
error2=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of r * delta v
fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
error3=integral2(fun,x1,x2,y1,y2);
v_norm = v_norm + error3;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1+error2+error3;
        
%Computing the conservation integral
fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
q_norm = @(x,y) k(x,y).*fun(x,y);
energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);

%Computing energy error of the flux er
fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
error2_f = integral2(fun,x1,x2,y1,y2);
fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
error3_f = integral2(fun,x1,x2,y1,y2);
energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

%Computing energy norm
fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
error_norm = integral2(fun,x1,x2,y1,y2);
energy_norm = energy_norm + error2_f - error_norm + error3;

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

b4 = (x1*f(x1,y1))/(x1+y2);
b2 = f(x1,y1) -b4;
b3 = -y2*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
error1=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of delta v * delta v
fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
error2=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of r * delta v
fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
error3=integral2(fun,x1,x2,y1,y2);
v_norm = v_norm + error3;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1+error2+error3;
        
%Computing the conservation integral
fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
q_norm = @(x,y) k(x,y).*fun(x,y);
energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);

%Computing energy error of the flux er
fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
error2_f = integral2(fun,x1,x2,y1,y2);
fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
error3_f = integral2(fun,x1,x2,y1,y2);
energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

%Computing energy norm
fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
error_norm = integral2(fun,x1,x2,y1,y2);
energy_norm = energy_norm + error2_f - error_norm + error3;

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

b4 = (f(x1,y1)*x2)/(x2+y1);
b2 = f(x1,y1) -b4;
b3 = -y1*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
error1=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of delta v * delta v
fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
error2=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of r * delta v
fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
error3=integral2(fun,x1,x2,y1,y2);
v_norm = v_norm + error3;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1+error2+error3;
        
%Computing the conservation integral
fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
q_norm = @(x,y) k(x,y).*fun(x,y);
energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);

%Computing energy error of the flux er
fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
error2_f = integral2(fun,x1,x2,y1,y2);
fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
error3_f = integral2(fun,x1,x2,y1,y2);
energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

%Computing energy norm
fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
error_norm = integral2(fun,x1,x2,y1,y2);
energy_norm = energy_norm + error2_f - error_norm + error3;


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

b4 = (f(x1,y1)*x1)/(x1+y1);
b2 = f(x1,y1) -b4;
b3 = -y1*b4;
b1 = b3;

% this error correspond to the integral in that cell of r*r
fun = @(x,y) ((b1+b2*x).^2+(b3+b4*y).^2)./k(x,y);
error1=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of delta v * delta v
fun = @(x,y) 2*((b1+b2*x).*(a2+a4*y)+(b3+b4*y).*(a3+a4*x));
error2=integral2(fun,x1,x2,y1,y2);
% This error correspond to the integral in that cell of r * delta v
fun = @(x,y) k(x,y).*((a2+a4*y).^2+(a3+a4*x).^2);
error3=integral2(fun,x1,x2,y1,y2);
v_norm = v_norm + error3;

% Add error at each cell plus the previously computed
energy_error=energy_error+error1+error2+error3;
        
%Computing the conservation integral
fun = @(x,y) (f(x,y)-b2-b4).^2 +x -x +y -y;
conservation_integral = conservation_integral + integral2(fun,x1,x2,y1,y2);
q_norm = @(x,y) k(x,y).*fun(x,y);
energy_q_norm = energy_q_norm + integral2(q_norm, x1,x2,y1,y2);

%Computing energy error of the flux er
fun = @(x,y) k(x,y).*((u_1(x,y)).^2 + (u_2(x,y)).^2);
error2_f = integral2(fun,x1,x2,y1,y2);
fun = @(x,y) 2*(u_1(x,y).*(b1+b2*x)+u_2(x,y).*(b3+b4*y));
error3_f = integral2(fun,x1,x2,y1,y2);
energy_error_flux = energy_error_flux + error1 + error2_f + error3_f;

%Computing energy norm
fun = @(x,y) 2*k(x,y).*(u_1(x,y).*(a2+a4*y)+u_2(x,y).*(a3+a4*x));
error_norm = integral2(fun,x1,x2,y1,y2);
energy_norm = energy_norm + error2_f - error_norm + error3;

%-------------- END OF COMPUTATIONS IN CELLS

%Final computations
conservation_integral = sqrt(conservation_integral);
energy_error = sqrt(energy_error);
energy_error_flux = sqrt(energy_error_flux);
energy_norm = sqrt(energy_norm);
v_norm = sqrt(v_norm);
energy_q_norm = sqrt(energy_q_norm);

end