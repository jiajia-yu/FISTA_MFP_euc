clear
load('../results/gaussian1d_compall');
format shortEng
%% FISTA
fprintf('fista\n');
[res1,res2] = kkt_ot1d(rho_fista,mx_fista);
disp(res1');disp(res2')

%% ALG
fprintf('\nAugmented Lagrangian\n');
[res1,res2] = kkt_ot1d_alg(rho_alg,mx_alg,outs_alg.phi,outs_alg.a,outs_alg.b);
disp(res1');disp(res2')

%% G-prox
fprintf('\nG-prox\n');
[res1,res2] = kkt_ot1d_gprox(rho_gprox,mx_gprox,outs_gprox.phi);
disp(res1');disp(res2')

%% ML
fprintf('\nmultilevel\n');
[res1,res2] = kkt_ot1d(rho_ml,mx_ml);
disp(res1');disp(res2')

%% MG5
fprintf('\nmultigrid5\n');
[res1,res2] = kkt_ot1d(rho_mg5,mx_mg5);
disp(res1');disp(res2')

%% MG10
fprintf('\nmultigrid10\n');
[res1,res2] = kkt_ot1d(rho_mg10,mx_mg10);
disp(res1');disp(res2')


