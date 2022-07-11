% test 
clear
% close all
clc

load('linear8-32-5_50000','rhos','mxs');
rhos_old = rhos;
mxs_old = mxs;
%% eg
eg = 'linear8-32-5_50000';
nt0 = 8;
nx0 = 32;
p_max = 5;
k = 1; 

%% options
opts = [];

opts.maxit = 5e4;
opts.tol = 1e-11;
opts.L0 = 10;
opts.eta = 1;
opts.sub_maxit = 1;

%%
rhos = cell(p_max,1);
mxs = cell(p_max,1);
outs = cell(p_max,1);
for p = 1:p_max
%% problem setting
nt = nt0*2^p; ntp = nt+1; ntm = nt-1;
dt = 1/nt;
t = linspace(0,1,ntp)';
opts.nt = nt;

nx = nx0*2^p; nxp = nx+1; nxm = nx-1;
dx = 1/nx;
x = linspace(0,1,nxp);

fprintf('nt=%d, nx=%d, ',nt,nx);

rho0 = k*(x-1/2)+1;
rho1 = ones(size(x));
rho0 = (rho0(1:end-1)+rho0(2:end))/2;
rho1 = (rho1(1:end-1)+rho1(2:end))/2;

opts.rho = rhos_old{p}(2:end-1,:);
opts.mx = mxs_old{p}(:,2:end-1);
%%
tic
[rhos{p},mxs{p},outs{p}] = ot1d_fista(rho0,rho1,opts);
t_fista = toc;
nit_fista = length(outs{p}.objs);
fprintf('%d iterations in %.2f s with residue %.2e\n',nit_fista,t_fista,outs{p}.ress(end));

end

save([eg,'_',num2str(opts.maxit)])

%% check convergence rate
fprintf('\n');
rho_errors = cell(p_max,1);
mx_errors = cell(p_max,1);
for p = 1:p_max
    % --- mesh size
    nt = nt0*2^p; ntp = nt+1; 
    dt = 1/nt;
    nx = nx0*2^p; nxp = nx+1; 
    dx = 1/nx;
    
    % --- rho_exact
    t = linspace(0,1,ntp)';
    x = linspace(dx/2,1-dx/2,nx);
    T = repmat(t,1,nx);
    X = repmat(x,ntp,1);

    sqrtterm = sqrt(2*T.*X+(T/2-1).^2);
    rho_exact = (sqrtterm + T - 1)./(T.*sqrtterm);
    rho_exact(1,:) = x + 0.5;

    % --- mx_exact
    t = linspace(dt/2,1-dt/2,nt)';
    x = linspace(0,1,nxp);
    T = repmat(t,1,nxp);
    X = repmat(x,nt,1);

    sqrtterm = sqrt(2*T.*X+(T/2-1).^2);
    Tsq = T.^2;
    Tcub = T.^3;
    mx_exact = X./Tsq + (T-3)./(2*Tcub).*sqrtterm ...
             - (T-1).*(Tsq-4)./(8*Tcub)./sqrtterm - (3*T-4)./(2*Tcub);
    
    [res_stat,res_feas] = kkt_ot1d(rho_exact,mx_exact);
    fprintf('nt=%d, nx=%d, residue of exact solution: ', nt,nx);
    disp([res_stat(:)',res_feas(:)']);
    rho_errors{p} = rhos{p} - rho_exact;
    mx_errors{p} = mxs{p} - mx_exact;
    
%     figure(1); mesh(rho_errors{p});title(['rho error, nt=',num2str(nt),...
%                                                     ' nx=',num2str(nx)]);
%     figure(2); mesh(mx_errors{p}); title(['mx  error, nt=',num2str(nt),...
%                                                     ' nx=',num2str(nx)]);
    
%     figure(1); 
%     for k = 1:nx
%         plot(linspace(0,1,ntp),rhos{p}(:,k)-rho_exact(:,k));
%         title(['rho error at x=',num2str((k-0.5)*dx)]);
%         pause(0.1);
%     end
%     
%     figure(2);
%     for k = 1:nt
%         plot(linspace(0,1,nxp),mxs{p}(k,:)-mx_exact(k,:));
%         title(['mx error at t=',num2str((k-0.5)*dt)]);
%         pause(0.1);
%     end
end

%---
W2_exact = 1/120;
maxerrors = zeros(1,p_max);
L2errors = zeros(1,p_max);
W2errors = zeros(1,p_max);
for p = 1:p_max
    nt = nt0*2^p; dt = 1/nt;
    nx = nx0*2^p; dx = 1/nx;
    
    rho_error = rho_errors{p};
    mx_error = mx_errors{p};
    maxerrors(p) = max( max(abs(rho_error(:))),max(abs(mx_error(:))) );
    L2errors(p) = ( norm(rho_error(:))+norm(mx_error(:)) )*sqrt(dt*dx);    
    W2errors(p) = abs( outs{p}.objs(end) - W2_exact );
    
    fprintf('nt=%d, nx=%d, max error=%f, L2 error=%e, W2 error=%e\n',...
        nt,nx,maxerrors(p),L2errors(p),W2errors(p));

end

% -diff(log(maxerrors))/log(2)
% -diff(log(L2errors))/log(2)
fprintf('\n');
fprintf('error max norm order: ');
disp(log(maxerrors(1:end-1)./maxerrors(2:end))/log(2));
fprintf('error L2 norm order: ');
disp(log(L2errors(1:end-1) ./L2errors(2:end) )/log(2));
fprintf('W2^2 error order: ');
disp(log(W2errors(1:end-1) ./W2errors(2:end) )/log(2));




