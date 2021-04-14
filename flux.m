%Calculate flux at cells
function [F_1, F_2] = flux(G,K,u,cells, num_cells_x,num_cells_y) 
F=-(K.*G)*u;

F_1 = zeros(length(cells(:,1)),1);
F_2 = F_1;
for j = 2:num_cells_y-1
    for i = 2:num_cells_x-1
    F_1(num_cells_x*(j-1)+i) = 0.5*(F(cells(num_cells_x*(j-1)+i,4))+F(cells(num_cells_x*(j-1)+i,6)));
    F_2(num_cells_x*(j-1)+i) = 0.5*(F(cells(num_cells_x*(j-1)+i,3))+F(cells(num_cells_x*(j-1)+i,5)));
    end
end
for i = 2:num_cells_x-1
    F_2(i) = 0.5*F(cells(i,5)); 
    F_2(i+(num_cells_y-1)*num_cells_x) = 0.5*F(cells(i+(num_cells_y-1)*num_cells_x,3));
end
for i = 2:num_cells_y-1
    F_1((i-1)*num_cells_x+1) = 0.5*F(cells((i-1)*num_cells_x+1,4)); 
    F_1(i*num_cells_x) = 0.5*F(cells(i*num_cells_x,6));
end