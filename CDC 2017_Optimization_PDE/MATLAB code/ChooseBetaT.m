% Tunable parameter: betaT
% n: number of unstable eigenmode
% Choose betaT to find the desirable number of unstable eigenmodes
% betaT = 80   for n = 1


clear;clc;

betaU = 2;
gamma = 4;
betaT = 80;

test = zeros(1,5);
n = 0;
for i = 1:5
    test(1,i) = - i^2 - betaU + betaT*gamma*exp(-gamma);
    if test(1,i) > 0
        n = n + 1;
    end
end
n
% choice: betaT = 200