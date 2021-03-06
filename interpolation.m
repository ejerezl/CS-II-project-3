% U = matrix of potentials
% F = matrix of fluxes
% x1,x2, y1, y2 = vertices of the cell in which we are
function [a2,a3,a4,b1,b2,b3,b4] = interpolation(flag, cell, F, dx, dy, A, B, C, D,x,y)
    if flag == 0
        %In this case only compute the interpolation of the potential
        %a1 = x*y*(A-B-C+D)/(dx*dy) - y*(-A+C)/dy - x*(B-A)/dx + A;
        a2 = (B-A+y*(B-A-D+C)/dy)/dx;
        a3 = (C-A+x*(C-D-A+B)/dx)/dy;
        a4 = (D-C-B+A)/(dx*dy);
        b1 = 0;
        b2 = 0;
        b3 = 0;
        b4 = 0;
        return
    end
    %a1 = x*y*(A-B-C+D)/(dx*dy) - y*(-A+C)/dy - x*(B-A)/dx + A;
    a2 = (B-A+y*(B-A-D+C)/dy)/dx;
    a3 = (C-A+x*(C-D-A+B)/dx)/dy;
    a4 = (D-C-B+A)/(dx*dy);
    b1 = F(cell(6))+x*(F(cell(6))-F(cell(4)))/dx;
    b2 = (F(cell(4))-F(cell(6)))/dx;
    b3 = F(cell(3))+y*(F(cell(3))-F(cell(5)))/dy;
    b4 = (F(cell(5))-F(cell(3)))/dy;
end