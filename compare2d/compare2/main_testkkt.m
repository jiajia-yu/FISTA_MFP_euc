clear
load('../results/gaussian2d_compall');
format shortEng
%% FISTA
fprintf('fista\n');
[res1,res2] = kkt_ot2d(rho_fista,mx_fista,my_fista);
disp(res1');disp(res2')

%% ALG
fprintf('\nAugmented Lagrangian\n');
[res1,res2] = kkt_ot2d_alg(rho_alg,mx_alg,my_alg,outs_alg.phi,outs_alg.a,outs_alg.b);
disp(res1');disp(res2')

%% G-prox
fprintf('\nG-prox\n');
[res1,res2] = kkt_ot2d_gprox(rho_gprox,mx_gprox,my_gprox,outs_gprox.phi);
disp(res1');disp(res2')

%% ML
fprintf('\nmultilevel\n');
[res1,res2] = kkt_ot2d(rho_ml,mx_ml,my_ml);
disp(res1');disp(res2')

%% MG5
fprintf('\nmultigrid5\n');
[res1,res2] = kkt_ot2d(rho_mg5,mx_mg5,my_mg5);
disp(res1');disp(res2')

%% MG10
fprintf('\nmultigrid10\n');
[res1,res2] = kkt_ot2d(rho_mg10,mx_mg10,my_mg10);
disp(res1');disp(res2')


