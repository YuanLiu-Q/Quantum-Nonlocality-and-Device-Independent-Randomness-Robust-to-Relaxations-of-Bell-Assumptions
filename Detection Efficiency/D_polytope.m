function D_mat = D_polytope(m, k)
% D_polytope:  Generate all LHV vertices for a bipartite Bell scenario with
%       - m inputs for Alice and Bob (x,y = 1,...,m),
%       - k outputs for each party (a,b = 1,...,k).
%
%   Each cell of the output D_mat{n} is an (m*k)-by-(m*k) matrix 
%   representing one deterministic behavior p(a,b | x,y)

A = combos_lex_lastfast(m, k);
B = combos_lex_lastfast(m, k);

D_mat = {};

len = size(A, 1);
for s_a = 1:len % output strategy of A
    for s_b = 1:len % output strategy of B
        oA = A(s_a, :);
        oB = B(s_b, :);
        curr = [];
        for i = 1:m
            curr_row = [];
            for j = 1:m
                curr_xy = zeros(k, k);
                curr_xy(oA(i), oB(j)) = 1;
                curr_row = [curr_row, curr_xy];
            end
            curr = [curr; curr_row];
        end
        D_mat{end+1} = curr;
    end
end

function A = combos_lex_lastfast(m, k)
    C = cell(1, m);
    [C{:}] = ndgrid(1:k);
    A = reshape(cat(m+1, C{:}), [], m);
end

end
