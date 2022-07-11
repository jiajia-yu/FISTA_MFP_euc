%% problem setting
clear
% clc

nt = 32; ntp = nt+1; ntm = nt-1;
% nt = 64; ntp = nt+1; ntm = nt-1;
dt = 1/nt;

eg = 'gaussian2d';
meanvalue0 = [1/3,2/3]; sigma0 = [1e-2,0;0,1e-2];
meanvalue1 = [2/3,1/3]; sigma1 = [1e-2,0;0,1e-2];
Gaussian = @(x,y, mu_x,mu_y,sigma) exp(-((x-mu_x).^2+(y-mu_y).^2)/(2*sigma^2))/(2*pi*sigma);
% meanvalue0 = [1/3,2/3]; sigma0 = [5e-2,0;0,5e-2];
% meanvalue1 = [2/3,1/3]; sigma1 = [5e-2,0;0,5e-2];
% Normal = @(x,y,meanvalue,sigmainv) sqrt(det(sigmainv))/(2*pi) * ...
%     exp(-0.5* (sigmainv(1,1)*(x-meanvalue(1)).^2 + ...
%                sigmainv(1,2)*(x-meanvalue(1)).*(y-meanvalue(2)) + ...
%                sigmainv(2,2)*(y-meanvalue(2)).^2));

%% main iteration
nx = 24; ny = 24;
% nx = 64; ny = 64;
p_max = 3;

fprintf(' nx  ny  FISTA                              ALG                                G-prox\n')
fprintf('         NumIt Time(s) Time(s)/Iter  Res    NumIt Time(s) Time(s)/Iter  Res    NumIt Time(s) Time(s)/Iter  Res  \n')  
for p = 1:p_max
    
    nx = nx*2; nxp = nx+1; nxm = nx-1;
    dx = 1/nx;
    fprintf('%3d ',nx);

    ny = ny*2; nyp = ny+1; nym = ny-1;
    dy = 1/ny;
    fprintf('%3d ',ny);

    x = reshape( linspace(0,1,nxp),nxp,1);
    x = repmat(x, 1,nyp);
    y = reshape( linspace(0,1,nyp),1,nyp);
    y = repmat(y, nxp,1);

%     rho0 = Normal(x,y,meanvalue0,inv(sigma0))+0.1;
%     rho1 = Normal(x,y,meanvalue1,inv(sigma1))+0.1;
    rho0 = Gaussian(x,y,1/3,2/3,0.1)+0.1;
    rho1 = Gaussian(x,y,2/3,1/3,0.1)+0.1;
    
    rho0 = (rho0(1:end-1,:)+rho0(2:end,:))/2;
    rho0 = (rho0(:,1:end-1)+rho0(:,2:end))/2;
    rho1 = (rho1(1:end-1,:)+rho1(2:end,:))/2;
    rho1 = (rho1(:,1:end-1)+rho1(:,2:end))/2;

    %% options
    opts = [];
    opts.nt = nt;

    opts.maxit = 1e3;
    opts.tol = 1e-3;

    %% FISTA
    opts.L0 = 8;
    opts.eta = 1.5;
    opts.sub_maxit = 5;
    opts.maxit = 1e3;
    tic
    [rho_fista,mx_fista,my_fista,outs_fista] = ot2d_fista(rho0,rho1,opts);
    t_fista = toc;
    nit_fista = length(outs_fista.objs);
    [res1,res2] = kkt_ot2d(rho_fista,mx_fista,my_fista);
    res_fista = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e ',nit_fista,t_fista,t_fista/nit_fista,res_fista);
    

    %%
    opts.maxit = 2e2;
    tic
    [rho_alg,mx_alg,my_alg,outs_alg] = ot2d_alg(rho0,rho1,opts);
    t_alg = toc;
    nit_alg = length(outs_alg.objs);
    [res1,res2] = kkt_ot2d_alg(rho_alg,mx_alg,my_alg,outs_alg.phi,outs_alg.a,outs_alg.b);
    res_alg = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e ',nit_alg,t_alg,t_alg/nit_alg,res_alg);
    
    
    %% G-prox
    opts.maxit = 2e2;
    tic
    [rho_gprox,mx_gprox,my_gprox,outs_gprox] = ot2d_gprox(rho0,rho1,opts);
    t_gprox = toc;
    nit_gprox = length(outs_gprox.objs);
    [res1,res2] = kkt_ot2d_gprox(rho_gprox,mx_gprox,my_gprox,outs_gprox.phi);
    res_gprox = max(abs([res1(2,:),res2(2)]));
    fprintf('%3d    %3.2f    %3.2e  %3.2e \n',nit_gprox,t_gprox,t_gprox/nit_gprox,res_gprox);
    
    %%
%     save(['../results/',eg,'_',num2str(nt),'_',num2str(nx)])
    
end
