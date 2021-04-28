%Interpolating potential
function [v] = inter_pot(x, pot, grid, num_cells)

% Interpolate
for i=1:num_cells-1
    if (grid(i)<= x) && (x <= grid(i+1))
        a0 = pot(i)+pot(i+1);
        a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
        v = a0+a1*(x-grid(i));
    end
end

% function [v] = inter_pot(x,pot, num_cells)
% 
% % Determine where x is
% for i=1:num_cells-1
%      if (grid(i)<= x) && (x <= grid(i+1))
%          a0 = pot(i)+pot(i+1);
%          a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
%          gridl = grid(i);
%          gridr = grid(i+1);
%      end
% end

% function v = inter_pot(x)%,pot, grid, num_cells)
% nx = 17;
% dx = 1/(nx-1);
% num_edges = nx-1;
% %edges
% edges = zeros(num_edges, 3);
% for i = 1:num_edges
%     edges(i,:) = [0.5*dx + (i-1)*dx, i, i+1];
% end
% %define new grid of edges which includes 0 and 1
% grid = zeros(nx+1,1);
% grid(2:(nx)) = edges(:,1,1);
% grid(end) = 1;
% %define potential on grid points as averages of the adjacents cells
% pot = zeros(size(grid));
% pot(1) = 0; %Boundary condition
% for i=2:grid-1
%     pot(i) = 0.5*(u(i)+u(i+1));
% end
% pot(end) = 0; %Boundary condition
% 
% % Interpolate
% for i=1:nx-1
%     if (grid(i)<= x) <= grid(i+1)
%         a0 = pot(i)+pot(i+1);
%         a1 = (pot(i+1)-pot(i))/(grid(i+1)-grid(i));
%         v = a0+a1*(x-grid(i));
%     end
% end
% end

