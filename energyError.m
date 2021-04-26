% Function that computes the error between the real solution and an
% approximation
% U = matrix with potentials
% F = flux
% k = from PDE (considering it a constant)
function error = energyError(U, F, k, dx, dy, cells, cells_x, cells_y)

error = 0;

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
        
        [a2, a3, a4, b1, b2, b3, b4] = interpolation(cells(cell,:),F,dx,dy,A,B,C,D,x1,y1);
        % this error correspond to the integral in that cell of r*r
        error1=(b1^2+b3^2)*dx*dy+b1*b2*(x2^2-x1^2)*dy+b3*b4*dx*(y2^2-y1^2)+(b2^2*(x2^3-x1^3)*dy+b4^2*dx*(y2^3-y1^3))/3;
        % This error correspond to the integral in that cell of delta v * delta v
        error2=(a2^2+a3^2)*dx*dy+a2*a4*dx*(y2^2-y1^2)+a3*a4*(x2^2-x1^2)*dy+a4^2*((x2^3-x1^3)*dy+dx*(y2^3-y1^3))/3;
        % This error correspond to the integral in that cell of r * delta v
        error3=(b1*a2+b3*a3)*dx*dy+(x2^2-x1^2)*dy*(b2*a2+b3*a4)/2+dx*(y2^2-y1^2)*(b1*a4+b4*a3)/2+(x2^2-x1^2)*(y2^2-y1^2)*(a4*b2+a4*b4)/4;

        % Add error at each cell plus the previously computed
        error=error+error1/k+k*error2+2*error3;
    end
end
end