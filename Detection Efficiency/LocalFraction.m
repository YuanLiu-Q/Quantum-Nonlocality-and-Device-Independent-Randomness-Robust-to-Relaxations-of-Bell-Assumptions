function lf = LocalFraction(p_cell)
% LocalFraction:  Compute the local fraction of the behavior p_cell.
%   - p_cell{x,y} is a |A|-by-|B| matrix of probabilities p(a,b|x,y)
% Dependencies: D_polytope.m and CVX

m = size(p_cell, 1);
k = size(p_cell{1,1}, 1);
D = D_polytope(m, k);

s = k * m;
p_mat = [];
for i = 1:m
    p_row = [];
    for j = 1:m
        p_row = [p_row, p_cell{i,j}];
    end
    p_mat = [p_mat; p_row];
end

cvx_begin quiet
    variables m(s, s)
    minimize sum(sum(m .* p_mat))
    subject to
        for i = 1:length(D)
            sum(sum(m .* D{i})) >= 1;
        end
        for i = 1:s
            for j = 1:s
                m(i,j) >= 0;
            end
        end
cvx_end

lf = cvx_optval;
end

