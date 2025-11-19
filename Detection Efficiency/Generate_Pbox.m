function p_cell = Generate_Pbox(ep)
% Generate_Pbox:  Generate the two-qubit quantum behavior for a given epsilon.
%       - ep is the scalar parameter epsilon, 0 <= ep < 1
%
%   Output p_cell{x,y} is a 2-by-2 matrix of probabilities p(a,b|x,y)


sigma_x = [0,1;1,0];
sigma_z = [1,0;0,-1];
theta = asin(3 - sqrt(5 + 4*ep));
state = cos(theta/2)*kron([1;0],[1;0]) - sin(theta/2)*kron([0;1],[0;1]);
rho = state*state';

A{1} = -(2+sin(theta))*sqrt(1-sin(theta))/(2-sin(theta))/sqrt(1+sin(theta))*sigma_z + ...
       -sqrt(2)*sin(theta)*sqrt(sin(theta))/(2-sin(theta))/sqrt(1+sin(theta))*sigma_x;
B{1} = A{1};

A{2} = -sqrt(1-sin(theta))/sqrt(1+sin(theta))*sigma_z + ...
        sqrt(2)*sqrt(sin(theta))/sqrt(1+sin(theta))*sigma_x;
B{2} = A{2};

I = eye(2);

p_cell = cell(2,2);
for i = 1:2
    for j = 1:2
        box_curr = zeros(2,2);
        box_curr(1,1) = 1/4*(1 + trace(kron(A{i},I)*rho) + trace(kron(I,B{j})*rho) + trace(kron(A{i},B{j})*rho));
        box_curr(1,2) = 1/4*(1 + trace(kron(A{i},I)*rho) - trace(kron(I,B{j})*rho) - trace(kron(A{i},B{j})*rho));
        box_curr(2,1) = 1/4*(1 - trace(kron(A{i},I)*rho) + trace(kron(I,B{j})*rho) - trace(kron(A{i},B{j})*rho));
        box_curr(2,2) = 1/4*(1 - trace(kron(A{i},I)*rho) - trace(kron(I,B{j})*rho) + trace(kron(A{i},B{j})*rho));
        p_cell{i,j} = box_curr;
    end
end

end
