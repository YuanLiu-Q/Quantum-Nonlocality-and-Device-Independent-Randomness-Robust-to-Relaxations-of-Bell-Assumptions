clear all
clc

% This script computes the quantum adversary Eve's guessing probability for
% Alice's first measurement outcomes, in the scenario that both Alice's and
% Bob's input information are partially leaked quantified by kappaA and 
% kappaB respectively, conditioned on an observed CHSH value s.
%
% The guessing probability is obtained by solving an SDP using the NPA
% hierarchy with localizing matrices that enforce the input-leakage
% constraints. 
%
% Dependencies: CVX and Moment (https://github.com/ajpgarner/moment)

% set the parameters kappaA and kappaB
kappaA = 0;
kappaB = 0.1;

% set the observed value for the CHSH Bell test
s = 2*sqrt(2);

% set the level for NPA hierarchy and local matrix level
level    = 3;
lm_level = 3;

% three-party scenario
setting = LocalityScenario(3);
Alice   = setting.Parties(1);
Bob     = setting.Parties(2);
Eve = setting.Parties(3);

% input and output settings
for i = 1:4
    Alice.AddMeasurement(2);
end
for i = 1:4
    Bob.AddMeasurement(2);
end
Eve.AddMeasurement(2);

% Make moment matrix
matrix = setting.MomentMatrix(level);

I = setting.id();

% differences between projectors Pi_{x,y=0}^a and Pi_{x,y=1}^a (Alice side)
c1 = setting.get(1) - setting.get(2) + kappaB * I;
c2 = -setting.get(1) + setting.get(2) + kappaB * I;
c3 = setting.get(3) - setting.get(4) + kappaB * I;
c4 = -setting.get(3) + setting.get(4) + kappaB * I;

% differences between projectors Pi_{y,x=0}^b and Pi_{y,x=1}^b (Bob side)
c5 = setting.get(5) - setting.get(6) + kappaA * I;
c6 = -setting.get(5) + setting.get(6) + kappaA * I;
c7 = setting.get(7) - setting.get(8) + kappaA * I;
c8 = -setting.get(7) + setting.get(8) + kappaA * I;

% localizing matrices corresponding to the above differences
l_matrix1 = c1.LocalizingMatrix(lm_level);
l_matrix2 = c2.LocalizingMatrix(lm_level);
l_matrix3 = c3.LocalizingMatrix(lm_level);
l_matrix4 = c4.LocalizingMatrix(lm_level);

l_matrix5 = c5.LocalizingMatrix(lm_level);
l_matrix6 = c6.LocalizingMatrix(lm_level);
l_matrix7 = c7.LocalizingMatrix(lm_level);
l_matrix8 = c8.LocalizingMatrix(lm_level);

% Make Alice–Bob behavior as a cell
BoxAB = cell(4,4);
for inputA = 1:4
    for inputB = 1:4
        box_curr = [];
        for outputA = 1:2
            box_row = [];
            for outputB = 1:2
                box_row = [box_row, ...
                    setting.getPMO([1,inputA,outputA]) * ...
                    setting.getPMO([2,inputB,outputB])];
            end
            box_curr = [box_curr; box_row];
        end
        BoxAB{inputA,inputB} = box_curr;
    end
end

% coefficients corresponding to the CHSH Bell test under information leakage
mat  = [1,-1; -1,1];
coee = cell(4,4);
for i = 1:4
    for j = 1:4
        coee{i,j} = zeros(2,2);
    end
end
coee{1,1} = mat;
coee{2,3} = mat;
coee{3,2} = mat;
coee{4,4} = -mat;

% guessing probability operators for Alice's first measurement outcome
p_guess = [ ...
    setting.getPMO([1,1,1]) * setting.getPMO([3,1,1]), ...
    setting.getPMO([1,1,2]) * setting.getPMO([3,1,2]) ...
];

% Define and solve SDP
cvx_begin sdp
    setting.cvxVars('a');
    M = matrix.Apply(a);
    M >= 0;
    l_matrix1.Apply(a) >= 0;
    l_matrix2.Apply(a) >= 0;
    l_matrix3.Apply(a) >= 0;
    l_matrix4.Apply(a) >= 0;
    l_matrix5.Apply(a) >= 0;
    l_matrix6.Apply(a) >= 0;
    l_matrix7.Apply(a) >= 0;
    l_matrix8.Apply(a) >= 0;
    a(1) == 1;

    chsh_val = 0;
    for inputA = 1:4
        for inputB = 1:4
            curr = BoxAB{inputA,inputB};
            for outputA = 1:2
                for outputB = 1:2
                    chsh_val = chsh_val + ...
                        coee{inputA,inputB}(outputA,outputB) * ...
                        curr(outputA,outputB).Apply(a);
                end
            end
        end
    end
    chsh_val == s;

    pg = p_guess(1).Apply(a) + p_guess(2).Apply(a);
    maximize(pg)
cvx_end
