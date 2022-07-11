function [rho,mx,outs] = ot1d_fista(rho0,rho1,opts)

%% parameters
% time and space parameters
nx = length(rho0);
nxm = nx-1;
dx = 1/nx;
if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
ntm = nt-1;
dt = 1/nt;

% inital value of rho and mx
if isfield(opts,'rho') xrho_km1 = opts.rho; else xrho_km1 = ones(ntm,nx); end
if isfield(opts,'mx')  xmx_km1  = opts.mx;  else xmx_km1  = ones(nt,nxm); end

% stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 5e3;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end

% parameter for backtracking
if isfield(opts,'L0')  L0  = opts.L0;  else L0  = 10; end
if isfield(opts,'eta') eta = opts.eta; else eta = 1.2; end 
if isfield(opts,'sub_maxit') sub_maxit = opts.sub_maxit; else sub_maxit = 5; end

% for solving Poisson equation
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,  nx);
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt, 1);
lap = repmat(lap_t,1,nx) + repmat(lap_x,nt,1);

%% initialization
objs = zeros(maxit,1);
ress = zeros(maxit,1);
projerrs = zeros(maxit,1);
stepsizes = zeros(maxit,1);

yrho = xrho_km1;
ymx = xmx_km1;
% obj_km1 = compObj(yrho,ymx);

t_km1 = 1;
%% main iteration
for nit = 1:maxit
    % backtracking
    [grad_rho,grad_mx] = compGrad(yrho,ymx);
    
    sub_nit = 0;
    L = L0;
    while sub_nit < sub_maxit
        [xrho_k,xmx_k,projerr] = compProj(yrho - 1/L*grad_rho, ymx - 1/L*grad_mx);
        
        obj_k = compObj(xrho_k,xmx_k);
        diff_rho = xrho_k - yrho;
        diff_mx  = xmx_k - ymx;
        G = compObj(yrho,ymx) + sum(grad_rho.*diff_rho,'all') ...
                    + sum(grad_mx .*diff_mx, 'all') ...
                + L/2*sum(diff_rho.*diff_rho,'all') ...
                + L/2*sum(diff_mx .*diff_mx, 'all');
        
%         G = obj_km1 + sum(grad_rho.*diff_rho,'all') ...
%                     + sum(grad_mx .*diff_mx, 'all') ...
%                 + L/2*sum(diff_rho.*diff_rho,'all') ...
%                 + L/2*sum(diff_mx .*diff_mx, 'all');
        
        if obj_k <= G
            break
        end
        sub_nit = sub_nit + 1;
        L = L*eta;
    end
        
    % update variables
    t_k = (1 + sqrt(1+4*t_km1^2))/2;
    w_k = (t_km1-1)/t_k;
    drho = xrho_k-xrho_km1;
    dmx = xmx_k -xmx_km1;
    yrho = max( xrho_k + w_k*drho,0.1);
%     yrho = xrho_k + w_k*drho;
    ymx  = xmx_k  + w_k*dmx;
    
    t_km1 = t_k;
%     obj_km1 = obj_k;
%     xrho_km1 = max( xrho_k,0.1);
    xrho_km1 = xrho_k;
    xmx_km1 = xmx_k;
    
    objs(nit) = obj_k*dt*dx;
    projerrs(nit) = projerr;
    stepsizes(nit) = 1/L;
%     L = max(min(L, 1e5), 1e-5);

    % stop criteria
    if mod(nit,10000)==0
        [res_stat,res_feas] = kkt_ot1d([rho0;xrho_k;rho1],...
                                       [zeros(nt,1),xmx_k,zeros(nt,1)]);
        res = max( max(res_stat(:)), max(res_feas(:)) );
        ress(nit) = res;
        if res < tol
            break 
        end
    end

    if isnan(obj_k) || isinf(obj_k)
        fprintf('blow up at iteration %d\n', nit);
        rho = [rho0;xrho_k;rho1];
        mx = [zeros(nt,1),xmx_k,zeros(nt,1)];
        outs.objs = objs(1:nit);
        outs.ress = ress(1:nit);
        outs.projerrs = projerrs(1:nit);
        outs.stepsizes = stepsizes(1:nit);
        return
    end
    
end

%% copy results
rho = [rho0;xrho_k;rho1];
mx = [zeros(nt,1),xmx_k,zeros(nt,1)];
outs.objs = objs(1:nit);
outs.ress = ress(1:nit);
outs.projerrs = projerrs(1:nit);
outs.stepsizes = stepsizes(1:nit);

%% functions

    function obj = compObj(rho,mx)
        rho = [rho0;rho;rho1];
        mx = [zeros(nt,1),mx,zeros(nt,1)];
        rho = It(rho);
        mx = Ix(mx);
        ind = rho > 1e-8;
        
        obj = sum( mx(ind).^2./rho(ind) );
    end

    function [grad_rho,grad_mx] = compGrad(rho,mx)
        rho = [rho0;rho;rho1];
        mx = [zeros(nt,1),mx,zeros(nt,1)];
        rho = It(rho);
        mx = Ix(mx);
        ind = rho > 1e-8;
        
        grad_rho = zeros(size(rho));
        grad_mx = zeros(size(mx));
        grad_rho(ind) = - mx(ind).^2./rho(ind).^2/2;
        grad_mx(ind) = mx(ind)./rho(ind);
        
        grad_rho = It(grad_rho);
        grad_mx = Ix(grad_mx);
    end

    function [rho_proj,mx_proj,projerr] = compProj(rho,mx)
        rho_proj = rho;
        mx_proj  = mx;
        rho = [rho0;rho;rho1];
        mx = [zeros(nt,1),mx,zeros(nt,1)];
        
        phi = mirt_dctn( Dt(rho)+Dx(mx) )./lap;
        phi(1,1) = 0;
        phi = mirt_idctn(phi);
        
        rho_proj = rho_proj + Dt(phi);
        mx_proj  = mx_proj  + Dx(phi);
        
        projerr = Dt([rho0;rho_proj;rho1])+Dx([zeros(nt,1),mx_proj,zeros(nt,1)]);
        projerr = max(abs(projerr),[],'all');
    end

    function DtA = Dt(A)
        DtA = (A(2:end,:) - A(1:end-1,:))/dt;
    end

    function DxA = Dx(A)
        DxA = (A(:,2:end) - A(:,1:end-1))/dx;
    end

%     function DtadjA = Dtadj(A)
%         DtadjA = [-A(1,:); ...
%                    A(1:end-1,:)-A(2:end,:); ...
%                    A(end,:)] /dt;
%     end
% 
%     function DxadjA = Dxadj(A)
%         DxadjA = cat(2, -A(:,1), ...
%                          A(:,1:end-1)-A(:,2:end), ...
%                          A(:,end)) /dx;
%     end

    function ItA = It(A)
        ItA = (A(1:end-1,:) + A(2:end,:))/2;
    end

    function IxA = Ix(A)
        IxA = (A(:,1:end-1) + A(:,2:end))/2;
    end

end