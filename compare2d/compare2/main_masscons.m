clear
load('../results/gaussian2d_compall');
format shortEng

mass = (sum(rho0s,'all')+sum(rho1s,'all'))/2;
%% FISTA
fprintf('fista\n');
disp(norm( sum(rho_fista,[2,3])/mass - 1 )*sqrt(dt))
disp(max( sum(rho_fista,[2,3])/mass - 1 ))

%% ALG
fprintf('\nAugmented Lagrangian\n');
disp(norm( sum(rho_alg,[2,3])/mass-1 )*sqrt(dt))
disp(max( sum(rho_alg,[2,3])/mass - 1 ))

%% G-prox
fprintf('\nG-prox\n');
disp(norm( sum(rho_gprox,[2,3])/mass - 1 )*sqrt(dt))
disp(max( sum(rho_gprox,[2,3])/mass - 1 ))

%% ML
fprintf('\nmultilevel\n');
disp(norm( sum(rho_ml,[2,3])/mass - 1 )*sqrt(dt))
disp(max( sum(rho_ml,[2,3])/mass - 1 ))

%% MG5
fprintf('\nmultigrid5\n');
disp(norm( sum(rho_mg5,[2,3])/mass - 1 )*sqrt(dt))
disp(max( sum(rho_mg5,[2,3])/mass - 1 ))

%% MG10
fprintf('\nmultigrid10\n');
disp(norm( sum(rho_mg10,[2,3])/mass - 1 )*sqrt(dt))
disp(max( sum(rho_mg10,[2,3])/mass - 1 ))

