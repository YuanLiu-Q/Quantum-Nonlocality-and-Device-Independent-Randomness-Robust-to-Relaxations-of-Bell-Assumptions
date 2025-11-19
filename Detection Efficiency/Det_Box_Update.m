function New_Box = Det_Box_Update(eta, p_cell)
% Det_Box_Update:  Update a behavior under finite detection efficiency.
%       - eta is the detection efficiency parameter , 0 <= eta <= 1
%       - p_cell{x,y} is a |A|-by-|B| matrix of probabilities p(a,b|x,y)
%
%   output New_Box{x,y} is a |A|-by-|B| matrix of probabilities p(a,b|x,y),
%   which is the effective behavior when taking into account the detection efficiency eta
%   and all non-detected events are assigned to the last outcome 

m = size(p_cell, 1);
New_Box = cell(m, m);

for i = 1:m
    for j = 1:m
        curr_p = p_cell{i,j};
        % update the probabilities for Alice (row-wise)
        curr_p = Update_x(eta, curr_p); 
        % update the probabilities for Bob (column-wise, via transpose)
        curr_p = curr_p';
        curr_p = Update_x(eta, curr_p);
        curr_p = curr_p';
        
        New_Box{i,j} = curr_p;
    end
end

    function output = Update_x(eta, p)
        % always assign non-click events to the last outcome
        k = size(p, 1);
        output = zeros(k, k);
        
        % clicked events for outcomes 1,...,k-1
        for ii = 1:k-1
            for jj = 1:k
                output(ii,jj) = p(ii,jj) * eta;
            end
        end
        
        % last outcome collects its clicked part plus all non-clicks
        p_marginal = sum(p, 1);
        for jj = 1:k
            output(k,jj) = p(k,jj) * eta + p_marginal(1,jj) * (1 - eta);
        end
    end

end
