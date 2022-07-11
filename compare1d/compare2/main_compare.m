% test 
clear
close all
% clc

%% problem setting
nt = 64; ntp = nt+1; ntm = nt-1;
dt = 1/nt;
t = linspace(0,1,ntp)';

nx = 256; nxp = nx+1; nxm = nx-1;
dx = 1/nx;
x = linspace(0,1,nxp);

% gaussian example
eg = 'gaussian1d';
meanvalue0 = 1/3; sigma0 = 0.1;
meanvalue1 = 2/3; sigma1 = 0.1;
Normal = @(x,meanvalue,sigma) 1/(sqrt(2*pi)*sigma)*exp(-0.5*((x-meanvalue)/sigma).^2);
rho0 = Normal(x,meanvalue0,sigma0)+0.1;
rho1 = Normal(x,meanvalue1,sigma1)+0.1;

rho0s = (rho0(1:end-1)+rho0(2:end))/2;
rho1s = (rho1(1:end-1)+rho1(2:end))/2;
% fprintf('theoretical W^2=%f\n',(meanvalue0-meanvalue1)^2+(sigma0-sigma1)^2);

%% options
opts = [];
opts.nt = nt;
opts.L0 = 8;
opts.eta = 1.2;

opts.tol = 1e-4;

%% FISTA
fprintf('fista\n');
opts.maxit = 1e3;

tic
[rho_fista,mx_fista,outs_fista] = ot1d_fista(rho0,rho1,opts);
t_fista = toc;
rho_fista = [ rho0s; rho_fista; rho1s];
mx_fista = [zeros(nt,1),mx_fista,zeros(nt,1)];

fprintf('%d iterations in %.3f sec\n',length(outs_fista.objs),t_fista);
fprintf('numerical W^2=%f\n',outs_fista.objs(end));

% figure;contour(rho_fista);title('FISTA');

%% ALG
fprintf('\nAugmented Lagrangian\n');
opts.maxit = 1e3;

tic
[rho_alg,mx_alg,outs_alg] = ot1d_alg(rho0,rho1,opts);
t_alg = toc;

fprintf('%d iterations in %.3f sec\n',length(outs_alg.objs),t_alg);
fprintf('numerical W^2=%f\n',outs_alg.objs(end));

% figure;contour(rho_alg);title('ALG');

%% G-prox
fprintf('\nG-prox\n');
opts.maxit = 1e3;

tic
[rho_gprox,mx_gprox,outs_gprox] = ot1d_gprox(rho0,rho1,opts);
t_gprox = toc;

fprintf('%d iterations in %.3f sec\n',length(outs_gprox.objs),t_gprox);
fprintf('numerical W^2=%f\n',outs_gprox.objs(end));

% figure;contour(rho_gprox);title('G-prox');

%% ML
fprintf('\nmultilevel\n');
opts.maxit = 500;
t_start = tic;
[rho_ml,mx_ml,outs_ml] = ot1d_fista_ml(rho0,rho1,4,opts);
t_gprox = toc(t_start);

rho_ml = [ rho0s; rho_ml; rho1s];
mx_ml = [zeros(nt,1),mx_ml,zeros(nt,1)];

fprintf('total: %d iterations in %.3f sec\n',outs_ml.nit_total,outs_ml.t_total);
fprintf('total time: %.3f sec\n',t_gprox);
fprintf('numerical W^2=%f\n',outs_ml.objs(end));

% figure;contour(rho_ml);title('multilevel')

%% MG5
fprintf('\nmultigrid5\n');
opts.maxit = 1;
opts.maxit_des = 5;
opts.maxit_aec = 500;
t_start = tic;
[rho_mg5,mx_mg5,outs_mg5] = ot1d_fista_mg(rho0,rho1,4,opts);
t_gprox = toc(t_start);
rho_mg5 = [ rho0s; rho_mg5; rho1s];
mx_mg5 = [zeros(nt,1),mx_mg5,zeros(nt,1)];

fprintf('total: %d iterations in %.3f sec\n',outs_mg5.nit_total,outs_mg5.t_total);
fprintf('total time: %.3f sec\n',t_gprox);
fprintf('numerical W^2=%f\n',outs_mg5.objs(end));

% figure;contour(rho_mg5);title('multigrid5')

%% MG10
fprintf('\nmultigrid10\n');
opts.maxit = 1;
opts.maxit_des = 10;
opts.maxit_aec = 500;
t_start = tic;
[rho_mg10,mx_mg10,outs_mg10] = ot1d_fista_mg(rho0,rho1,4,opts);
t_gprox = toc(t_start);
rho_mg10 = [ rho0s; rho_mg10; rho1s];
mx_mg10 = [zeros(nt,1),mx_mg10,zeros(nt,1)];

fprintf('total: %d iterations in %.3f sec\n',outs_mg10.nit_total,outs_mg10.t_total);
fprintf('total time: %.3f sec\n',t_gprox);
fprintf('numerical W^2=%f\n',outs_mg10.objs(end));
% figure;contour(rho_mg10);title('multigrid10')

%%
eg = ['../results/',eg];
save([eg,'_compall']);
% show_evolution(rho_fista,[eg,'_fista']);
% show_evolution(rho_alg,[eg,'_alg']);
% show_evolution(rho_gprox,[eg,'_gprox']);
% show_evolution(rho_ml,[eg,'_ml']);
% show_evolution(rho_mg5,[eg,'_mg5']);
% show_evolution(rho_mg10,[eg,'_mg10']);

