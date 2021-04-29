%Interpolating flux
function[r] = inter_flux(x, q, grid, num_edges) 

for i=1:num_edges+1
    if (grid(i)<= x) && (x <= grid(i+1))
        b0 = q(i,1);
        b1 = (q(i+1,1)-q(i,1))/(grid(i+1)-grid(i));
        r = b0+b1*(x-edges_new(i));
    end
end

