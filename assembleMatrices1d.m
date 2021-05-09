%assembling of matrices for 1D
%cell = intervall in R, first coord: position of midpoint, 2nd&3rd:
%adjacent edges
%edges = points between cells, first coord: position, 2nd&3rd: adjacent
%cells

function [A, b, G, D, K, cells, edges] = assembleMatrices1d(num_cells, f, k1)

%Define mesh
dx = 1/(num_cells - 1);
cells = zeros(num_cells, 3);

%relevant edges
num_edges = num_cells - 1;
edges = zeros(num_edges, 3);

%set node values
for i = 1:num_cells
    cells(i,1) = dx*(i-1);
end

%set edges connected to interior cells; 1st entry: left edge, 2nd
%entry: right edge
for i = 2:num_cells-1
    cells(i,2:3) = [i-1, i];
end

%set edges on boundary
cells(1,2:3) = [0,1];
cells(num_cells,2:3) = [num_edges,0];

%set edge values and connectivity to cells, 1st entry: position, 2nd entry:
%left cell, 3rd entry: right cell
for i = 1:num_edges
    edges(i,:) = [0.5*dx + (i-1)*dx, i, i+1];
end

% Grad operator
G = zeros(num_edges,num_cells);
for i=1:num_edges
    G(i,i)=1/dx;
end
for i=1:(num_edges-1)
    G(i,i+1) = -1/dx;
end
G = sparse(G);


%K operator
K = zeros(num_edges,num_edges);
for i=2:num_edges-1
    K(i,i) = 2/(1/k1(cells(edges(i,2),1))+1/k1(cells(edges(i,3),1)));
end
K(1,1) = k1(cells(edges(1,3),1));
K(num_edges,num_edges) = k1(cells(edges(num_edges,2),1));

D=-G';
A = -D*(K*G);

A(:,1) = 0;
A(1,:) = 0;
A(1,1) =1;
A(end,:) = 0;
A(:,end) = 0;
A(end,end) = 1;


%Create RHS
b = ones(num_cells,1);
for i=1:num_cells
    b(i) = f(cells(i));
end
b(1) = 0;
b(end) = 0;

%Hack boundary cells
% boundary_cells = zeros(1,2);
% boundary_cells(1,1) = 1;
% boundary_cells(1,2) = num_cells;
% for i=boundary_cells
%     A(i,i) = 1;
%     b(i) = 0;
% end




