%% Comparison of OT, MFP, MFG
clear
clc

% problem setting
eg = 'image_20';
nt = 20; ntp = nt+1; ntm = nt-1;
dt = 1/nt;

% image preprocessing
rho0 = imread('centaur1.bmp');
rho1 = imread('man1.bmp');
rho0 = imresize(rho0,[256,256]);
rho1 = imresize(rho1,[256,256]);
Qx = double(rho0>0&rho1>0);

rho0 = double(255-rho0)/255 +0.1;
rho1 = double(255-rho1)/255 +0.1;
Gx = -rho1;

t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile

imshow(rho0,[]);
colorbar
exportgraphics(t,'results\rho0.pdf','BackgroundColor','none')

imshow(rho1,[]);
colorbar
exportgraphics(t,'results\rho1.pdf','BackgroundColor','none')

imshow(Qx,[]);
colorbar
exportgraphics(t,'results\Qx.pdf','BackgroundColor','none')

imshow(Gx,[]);
colorbar
exportgraphics(t,'results\Rx.pdf','BackgroundColor','none')

% figure(1);imshow(rho0,[])
% figure(2);imshow(rho1,[])
% figure(3);imshow(Qx,[])
% figure(4);imshow(Gx,[])


%% MFG
fprintf('\nmean field game\n')
% options
opts_mfg = [];
opts_mfg.nt = nt;

opts_mfg.maxit = 5e2;
opts_mfg.tol = 1e-8;
opts_mfg.L0 = 10;
opts_mfg.eta = 1.5;
opts_mfg.sub_maxit = 2;

opts_mfg.lambda_L = 1;
opts_mfg.lambda_E = 0.01;
opts_mfg.lambda_P = 0.1;
opts_mfg.lambda_G = 1;

opts_mfg.Qx = Qx;
opts_mfg.Gx = Gx;

% optimization
tic
[rho_mfg,mx_mfg,my_mfg,outs_mfg] = func_mfg(rho0,rho1,opts_mfg);
t_mfg = toc;
nit_mfg = length(outs_mfg.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_mfg,t_mfg,t_mfg/nit_mfg);

% [res1,res2] = kkt_mfg(rho_mfg,mx_mfg,my_mfg,opts_mfg);
% disp(res1);
% disp(res2');

%% Normalize rho0, rho1
rho0 = imread('centaur1.bmp');
rho1 = imread('man1.bmp');
rho0 = imresize(rho0,[256,256]);
rho1 = imresize(rho1,[256,256]);
Qx = double(rho0>0&rho1>0);

rho0 = double(255-rho0)/255;
rho1 = double(255-rho1)/255;
rho0 = rho0/sum(rho0(:))*4e4 + 0.1;
rho1 = rho1/sum(rho1(:))*4e4 + 0.1;

%% OT
fprintf('\noptimal transport\n')
% options
opts_ot = [];
opts_ot.nt = nt;

opts_ot.maxit = 5e2;
opts_ot.tol = 1e-8;
opts_ot.L0 = 10;
opts_ot.eta = 1.5;
opts_ot.sub_maxit = 2;

% optimization
tic
[rho_ot,mx_ot,my_ot,outs_ot] = func_ot(rho0,rho1,opts_ot);
t_ot = toc;
nit_ot = length(outs_ot.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_ot,t_ot,t_ot/nit_ot);

% [res1,res2] = kkt_ot(rho_ot,mx_ot,my_ot);
% disp(res1);
% disp(res2');

%% MFP1 interaction \rho*log(\rho)
fprintf('\nmean field planning 1\n')
% options
opts_mfp1 = [];
opts_mfp1.nt = nt;

opts_mfp1.maxit = 2e2;
opts_mfp1.tol = 1e-8;
opts_mfp1.L0 = 100;
opts_mfp1.eta = 1.5;
opts_mfp1.sub_maxit = 10;

opts_mfp1.lambda_L = 1;
opts_mfp1.lambda_E = 0.01;%0.01;
opts_mfp1.lambda_P = 0.1;

opts_mfp1.Qx = Qx;

% optimization
tic
[rho_mfp1,mx_mfp1,my_mfp1,outs_mfp1] = func_mfp1(rho0,rho1,opts_mfp1);
t_mfp1 = toc;
nit_mfp1 = length(outs_mfp1.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_mfp1,t_mfp1,t_mfp1/nit_mfp1);

% [res1,res2] = kkt_mfp1(rho_mfp1,mx_mfp1,my_mfp1,opts_mfp1);
% disp(res1);
% disp(res2');

%% MFP2 interaction \rho*\rho/2
fprintf('\nmean field planning 2\n')
% options
opts_mfp2 = [];
opts_mfp2.nt = nt;

opts_mfp2.maxit = 5e2;
opts_mfp2.tol = 1e-8;
opts_mfp2.L0 = 10;
opts_mfp2.eta = 1.5;
opts_mfp2.sub_maxit = 5;

opts_mfp2.lambda_L = 1;
opts_mfp2.lambda_E = 0.01;%0.01;
opts_mfp2.lambda_P = 0.1;

opts_mfp2.Qx = Qx;

% optimization
tic
[rho_mfp2,mx_mfp2,my_mfp2,outs_mfp2] = func_mfp2(rho0,rho1,opts_mfp2);
t_mfp2 = toc;
nit_mfp2 = length(outs_mfp2.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_mfp2,t_mfp2,t_mfp2/nit_mfp2);

% [res1,res2] = kkt_mfp2(rho_mfp2,mx_mfp2,my_mfp2,opts_mfp2);
% disp(res1);
% disp(res2');

%% MFP3 interaction 1/\rho
fprintf('\nmean field planning3\n')
% options
opts_mfp3 = [];
opts_mfp3.nt = nt;

opts_mfp3.maxit = 5e2;
opts_mfp3.tol = 1e-8;
opts_mfp3.L0 = 10;
opts_mfp3.eta = 1.5;
opts_mfp3.sub_maxit = 5;

opts_mfp3.lambda_L = 1;
opts_mfp3.lambda_E = 0.01;%0.01;
opts_mfp3.lambda_P = 0.1;

opts_mfp3.Qx = Qx;

% optimization
tic
[rho_mfp3,mx_mfp3,my_mfp3,outs_mfp3] = func_mfp3(rho0,rho1,opts_mfp3);
t_mfp3 = toc;
nit_mfp3 = length(outs_mfp3.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_mfp3,t_mfp3,t_mfp3/nit_mfp3);

% [res1,res2] = kkt_mfp3(rho_mfp3,mx_mfp3,my_mfp3,opts_mfp3);
% disp(res1);
% disp(res2');

%%
save(['results\',eg])
% show_evolution(rho_mfg,['results\',eg,'_mfg']);
% show_evolution(rho_ot,['results\',eg,'_ot']);
% show_evolution(rho_mfp1,['results\',eg,'_mfp1']);
% show_evolution(rho_mfp2,['results\',eg,'_mfp2']);
% show_evolution(rho_mfp3,['results\',eg,'_mfp3']);
