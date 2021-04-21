num_cells = 5;
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
cells(1,2:3) = [1,0];
cells(num_cells,2:3) = [0,num_edges];

%set edge values and connectivity to cells, 1st entry: position, 2nd entry:
%left cell, 3rd entry: right cell
for i = 1:num_edges
    edges(i,:) = [0.5*dx + (i-1)*dx, i, i+1];
end

%Grad operator
grad_occupation_i = zeros(1,2*num_edges);
grad_occupation_j = grad_occupation_i;
grad_values = grad_occupation_i;
diffusion_values = grad_occupation_i;

for i = 1:(num_edges)
    grad_occupation_j((2*i):2*i+1) = i;
    grad_occupation_i((2*i):2*i+1) = edges(i,2:3); %different than erlands code!!
    grad_values((2*i):2*i+1) = [-1,1]/dx; 
    diffusion_values((2*i):2*i+1) = 2*[1,1]/(1/k1(cells(edges(i,2),1))+1/k1(cells(edges(i,3),1)));
end
disp(diffusion_values)

%Remove impact from cells on the boundary LEFT/RIGHT
grad_occupation_i(1) = 0;
grad_occupation_j(1) = 0;
grad_values(1) = 0;
%only right cell influences K
diffusion_values(1:2) = [0,1]*k1(cells(edges(1,3),1));

grad_occupation_i(end) = 0;
grad_occupation_j(end) = 0;
grad_values(end) = 0;
%only left cell influences K
diffusion_values(end-1:end) = [1,0]*k1(cells(edges(num_edges,2),1));

%Build G, K
grad_occupation_i(grad_occupation_i ==0)= [];
grad_occupation_j(grad_occupation_j ==0)= [];
grad_values(grad_values ==0)= [];
diffusion_values(diffusion_values ==0)= [];
disp(grad_occupation_i)
disp(grad_occupation_j)
disp(grad_values)
disp(diffusion_values)
G = sparse(grad_occupation_i,grad_occupation_j,grad_values,num_edges,num_cells);
disp(G)
K = sparse(grad_occupation_i,grad_occupation_j,diffusion_values,num_edges,num_cells);
disp(K)


%Div operator
div_occupation_i = zeros(1,2*num_cells);
div_occupation_j = div_occupation_i;
div_values = div_occupation_i;

%interior cells
    for i = 1:num_cells-1
        div_occupation_i(2*i-1:2*i) = i;
        div_occupation_j(2*i-1:2*i) = cells(i,2:3);
        div_values(2*i-1:2*i) = [1/dx, -1/dx];
    end


div_occupation_i(div_occupation_i ==0)= [];
div_occupation_j(div_occupation_j ==0)= [];
div_values(div_values ==0)= [];

D = sparse(div_occupation_i, div_occupation_j,div_values,num_cells, num_edges);
disp(D)