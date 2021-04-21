%Calculate flux at cells
function [F_1] = flux1d(G,K,u,cells, num_cells) 
F=-(K*G)*u;

F_1 = zeros(length(cells(:,1)),1);


for i = 2:num_cells-1
F_1(i) = 0.5*(F(cells(i,2))+F(cells(i,3)));
end
