%% problem setting
clear
close all
clc

eg = 'gaussian1d';
meanvalue0 = 1/3; sigma0 = 0.1;
meanvalue1 = 2/3; sigma1 = 0.1;
Normal = @(x,meanvalue,sigma) 1/(sqrt(2*pi)*sigma)*exp(-0.5*((x-meanvalue)/sigma).^2);

nt = 64; ntp = nt+1; ntm = nt-1;
dt = 1/nt;
t = linspace(0,1,ntp)';

%% main iteration

nx = 128;
p_max = 5;

fprintf('nx   FISTA                              ALG                                G-prox\n')
fprintf('     NumIt Time(s) Time(s)/Iter  Res    NumIt Time(s) Time(s)/Iter  Res    NumIt Time(s) Time(s)/Iter  Res  \n')  
for p = 1:p_max
    
    nx = nx*2; nxp = nx+1; nxm = nx-1;
    dx = 1/nx;
    x = linspace(0,1,nxp);
    fprintf('%4d ',nx);

    rho0 = Normal(x,meanvalue0,sigma0)+0.1;
    rho1 = Normal(x,meanvalue1,sigma1)+0.1;
    % fprintf('theoretical W^2=%f\n',(meanvalue0-meanvalue1)^2+(sigma0-sigma1)^2);

    rho0 = (rho0(1:end-1)+rho0(2:end))/2;
    rho1 = (rho1(1:end-1)+rho1(2:end))/2;

    %% options
    opts = [];
    opts.nt = nt;

    opts.tol = 1e-4;

    %% FISTA
    opts.L0 = 8;
    opts.eta = 1.2;
    opts.sub_maxit = 5;
    opts.maxit = 1e3;
    tic
    [rho_fista,mx_fista,outs_fista] = ot1d_fista(rho0,rho1,opts);
    t_fista = toc;
    nit_fista = length(outs_fista.objs);
    [res1,res2] = kkt_ot1d(rho_fista,mx_fista);
    res_fista = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e ',nit_fista,t_fista,t_fista/nit_fista,res_fista);
    
    
    %% augmented lagrangian
    opts.maxit = 1e3;
    tic
    [rho_alg,mx_alg,outs_alg] = ot1d_alg(rho0,rho1,opts);
    t_alg = toc;
    nit_alg = length(outs_alg.objs);
    [res1,res2] = kkt_ot1d_alg(rho_alg,mx_alg,outs_alg.phi,outs_alg.a,outs_alg.b);
    res_alg = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e ',nit_alg,t_alg,t_alg/nit_alg,res_alg);
    

    %% G-prox
    opts.maxit = 1e3;
    tic
    [rho_gprox,mx_gprox,outs_gprox] = ot1d_gprox(rho0,rho1,opts);
    t_gprox = toc;
    nit_gprox = length(outs_gprox.objs);
    [res1,res2] = kkt_ot1d_gprox(rho_gprox,mx_gprox,outs_gprox.phi);
    res_gprox = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e \n',nit_gprox,t_gprox,t_gprox/nit_gprox,res_gprox);
    
    %%
    save(['../results/',eg,'_',num2str(nt),'_',num2str(nx)])
    % show_evolution(rho_fista);
    % show_evolution(rho_alg);
    % show_evolution(rho_gprox);

end

