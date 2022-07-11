% test 
clear
close all
% clc

%% problem setting
nt = 64; ntp = nt+1; ntm = nt-1;
dt = 1/nt;

nx = 256; nxp = nx+1; nxm = nx-1;
dx = 1/nx;

ny = 256; nyp = ny+1; nym = ny-1;
dy = 1/ny;

x = reshape( linspace(0,1,nxp),nxp,1);
x = repmat(x, 1,nyp);
y = reshape( linspace(0,1,nyp),1,nyp);
y = repmat(y, nxp,1);

% gaussian example
eg = 'gaussian2d';
meanvalue0 = [1/3,2/3]; sigma0 = [0.05,0;0,0.05];
meanvalue1 = [2/3,1/3]; sigma1 = [0.05,0;0,0.05];
Normal = @(x,y,meanvalue,sigmainv) sqrt(det(sigmainv))/(2*pi) * ...
    exp(-0.5* (sigmainv(1,1)*(x-meanvalue(1)).^2 + ...
               sigmainv(1,2)*(x-meanvalue(1)).*(y-meanvalue(2)) + ...
               sigmainv(2,2)*(y-meanvalue(2)).^2));

rho0 = Normal(x,y,meanvalue0,inv(sigma0))+0.1;
rho1 = Normal(x,y,meanvalue1,inv(sigma1))+0.1;

rho0s = (rho0(1:end-1,:)+rho0(2:end,:))/2;
rho0s = (rho0s(:,1:end-1)+rho0s(:,2:end))/2;
rho0s = reshape(rho0s,1,nx,ny);
rho1s = (rho1(1:end-1,:)+rho1(2:end,:))/2;
rho1s = (rho1s(:,1:end-1)+rho1s(:,2:end))/2;
rho1s = reshape(rho1s,1,nx,ny);
% fprintf('theoretical W^2=%f\n',(meanvalue0-meanvalue1)^2+(sigma0-sigma1)^2);

%% options
opts = [];
opts.nt = nt;
opts.L0 = 8;
opts.eta = 1.2;

opts.tol = 1e-3;

%% FISTA
fprintf('fista\n');
opts.maxit = 5e2;

tic
[rho_fista,mx_fista,my_fista,outs_fista] = ot2d_fista(rho0,rho1,opts);
t_fista = toc;
rho_fista = cat(1,rho0s,rho_fista,rho1s);
mx_fista = cat(2,zeros(nt,1,ny),mx_fista,zeros(nt,1,ny));
my_fista = cat(3,zeros(nt,nx,1),my_fista,zeros(nt,nx,1));

fprintf('%d iterations in %.3f sec\n',length(outs_fista.objs),t_fista);
fprintf('numerical W^2=%f\n',outs_fista.objs(end));

%% ALG
fprintf('\nAugmented Lagrangian\n');
opts.maxit = 1e3;

tic
[rho_alg,mx_alg,my_alg,outs_alg] = ot2d_alg(rho0,rho1,opts);
t_alg = toc;

fprintf('%d iterations in %.3f sec\n',length(outs_alg.objs),t_alg);
fprintf('numerical W^2=%f\n',outs_alg.objs(end));

%% G-prox
fprintf('\nG-prox\n');
opts.maxit = 1e3;

tic
[rho_gprox,mx_gprox,my_gprox,outs_gprox] = ot2d_gprox(rho0,rho1,opts);
t_gprox = toc;

fprintf('%d iterations in %.3f sec\n',length(outs_gprox.objs),t_gprox);
fprintf('numerical W^2=%f\n',outs_gprox.objs(end));

%% ML
fprintf('\nmultilevel\n')
opts.maxit = 5e2;
t_start = tic;
[rho_ml,mx_ml,my_ml,outs_ml] = ot2d_fista_ml(rho0,rho1,4,opts);
t_total = toc(t_start);
rho_ml = cat(1,rho0s,rho_ml,rho1s);
mx_ml = cat(2,zeros(nt,1,ny),mx_ml,zeros(nt,1,ny));
my_ml = cat(3,zeros(nt,nx,1),my_ml,zeros(nt,nx,1));

fprintf('total: %d iterations in %.3f sec\n',outs_ml.nit_total,outs_ml.t_total);
fprintf('total time %.3f sec\n',t_total);
fprintf('numerical W^2=%f\n',outs_ml.objs(end));

% show_evolution2d(rho_ml,'contour');
% show_evolution2d(rho_ml,'mesh');

%% MG5
fprintf('\nmultigrid5\n');
opts.maxit = 1;
opts.maxit_des = 5;
opts.maxit_aec = 500;
t_start = tic;
[rho_mg5,mx_mg5,my_mg5,outs_mg5] = ot2d_fista_mg(rho0,rho1,4,opts);
t_total = toc(t_start);
rho_mg5 = cat(1,rho0s,rho_mg5,rho1s);
mx_mg5 = cat(2,zeros(nt,1,ny),mx_mg5,zeros(nt,1,ny));
my_mg5 = cat(3,zeros(nt,nx,1),my_mg5,zeros(nt,nx,1));

fprintf('total: %d iterations in %.3f sec\n',outs_mg5.nit_total,outs_mg5.t_total);
fprintf('total time %.3f sec\n',t_total);
fprintf('numerical W^2=%f\n',outs_mg5.objs(end));

% show_evolution2d(rho_mg5,'contour');
% show_evolution2d(rho_mg,'mesh');

%% MG10
fprintf('\nmultigrid10\n');
opts.maxit = 1;
opts.maxit_des = 10;
opts.maxit_aec = 500;
t_start = tic;
[rho_mg10,mx_mg10,my_mg10,outs_mg10] = ot2d_fista_mg(rho0,rho1,4,opts);
t_total = toc(t_start);
rho_mg10 = cat(1,rho0s,rho_mg10,rho1s);
mx_mg10 = cat(2,zeros(nt,1,ny),mx_mg10,zeros(nt,1,ny));
my_mg10 = cat(3,zeros(nt,nx,1),my_mg10,zeros(nt,nx,1));

fprintf('total: %d iterations in %.3f sec\n',outs_mg10.nit_total,outs_mg10.t_total);
fprintf('total time %.3f sec\n',t_total);
fprintf('numerical W^2=%f\n',outs_mg10.objs(end));

% show_evolution2d(rho_mg10,'contour');
% show_evolution2d(rho_mg,'mesh');

%%
eg = ['../results/',eg];
save([eg,'_compall'])

% showfunc = 'imshow';
% showfunc = 'contour';
% showfunc = 'mesh';
% 
% show_evolution(rho_fista,showfunc,[eg,'_fista']);
% show_evolution(rho_alg,showfunc,[eg,'_alg']);
% show_evolution(rho_gprox,showfunc,[eg,'_gprox']);
% show_evolution(rho_ml,showfunc,[eg,'_ml']);
% show_evolution(rho_mg5,showfunc,[eg,'_mg5']);
% show_evolution(rho_mg10,showfunc,[eg,'_mg10']);
