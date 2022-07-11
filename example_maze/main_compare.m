%% Comparison MFP with different parameters
clear
% close all
% clc

% problem setting
nt = 20; ntp = nt+1; ntm = nt-1;
dt = 1/nt;

nx = 256; nxp = nx+1; nxm = nx-1;
dx = 1/nx;

ny = 256; nyp = ny+1; nym = ny-1;
dy = 1/ny;

x = reshape( linspace(-1/2,1/2,nxp),nxp,1);
x = repmat(x, 1,nyp);
x = (x(1:end-1,:)+x(2:end,:))/2;
x = (x(:,1:end-1)+x(:,2:end))/2;
y = reshape( linspace(-1/2,1/2,nyp),1,nyp);
y = repmat(y, nxp,1);
y = (y(1:end-1,:)+y(2:end,:))/2;
y = (y(:,1:end-1)+y(:,2:end))/2;

% example
eg = 'maze_20';
rho0 = 10*exp(-( (x+0.3).^2 + (y+0.3).^2)*100 ) + 0.1;
rho1 = 10*exp(-( (x-0.3).^2 + (y-0.3).^2)*100 ) + 0.1;
% Qx = double((x>-0.1) & (x<-0.05) & (y<0)) + double((x>0.05) & (x<0.1) & (y>0));
% Qx = double((x>-0.02) & (x<0.02) & (y<-0.2));
Qx = double( x.^2+y.^2 < 0.02);
rho0(Qx>0) = 0;
rho1(Qx>0) = 0;
% Gx = -rho1;

%%
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile

imshow(rho0,[]);
colormap default
colorbar
exportgraphics(t,'results\rho30.pdf','BackgroundColor','none')

imshow(rho1,[]);
colormap default
colorbar
exportgraphics(t,'results\rho31.pdf','BackgroundColor','none')

imshow(Qx,[]);
colormap default
colorbar
exportgraphics(t,'results\Qx3.pdf','BackgroundColor','none')

eg_illustration(max(rho0,rho1),Qx);
exportgraphics(t,'results\eg3.pdf','BackgroundColor','none')

% imshow(Gx,[]);
% colormap default
% colorbar
% exportgraphics(t,'results\Rx.pdf','BackgroundColor','none')

%%
fprintf('\nmean field planning lambda_P=20000\n')
% options
opts = [];
opts.nt = nt;

opts.maxit = 1e3;
opts.tol = 1e-8;
opts.L0 = 1e3;
opts.eta = 2;
opts.sub_maxit = 2;

opts.lambda_L = 1;
opts.lambda_E = 0;%0.01;
opts.lambda_P = 1;

opts.Qx = 80000*Qx;

% optimization
tic
[rho_mfp,mx_mfp,my_mfp,outs_mfp] = func_mfp(rho0,rho1,opts);
t_mfp = toc;
nit_mfp= length(outs_mfp.objs);
fprintf('%d iterations in %.2f s -> %.2e s/iter\n',nit_mfp,t_mfp,t_mfp/nit_mfp);

% show_evolution(rho_mfp3,'imshow',['results\',eg,'_mfp3']);
figure; 
subplot(121); plot(outs_mfp.objs); 
subplot(122); semilogy(outs_mfp.ress);

[res1,res2] = kkt_mfp(rho_mfp,mx_mfp,my_mfp,opts);
disp(res1);
disp(res2');

%%
save('results\maze3');
% save(eg)

show_movement(rho_mfp,mx_mfp,my_mfp,Qx,[],'results\maze3')
% show_movement(rho_mfp,mx_mfp,my_mfp,[])
