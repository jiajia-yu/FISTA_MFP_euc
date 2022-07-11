function [rho, mx, outs] = ot1d_alg(rho0,rho1,opts)
rho0 = Ix(rho0);
rho1 = Ix(rho1);
%% parameters
% time and space parameters
nx = length(rho0);
nxp = nx+1;
dx = 1/nx;
if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
ntp = nt+1;
dt = 1/nt;

% initial value of iteration
if isfield(opts,'rho') rho = opts.rho; else rho = ones(ntp,nx); end
if isfield(opts,'mx')  mx = opts.mx;   else mx = ones(nt,nxp);  end
rho(1,:,:) = rho0; rho(end,:,:) = rho1;
mx(:,[1,end],:) = 0;

if isfield(opts,'a') a = opts.a; else a = ones(ntp,nx); end
if isfield(opts,'b') b = opts.b; else b = ones(nt,nxp); end
a([1,end],:) = 0;
b(:,[1,end]) = 0;

% number of iteration, penalty parameters and stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 500; end
if isfield(opts,'tol')   tol = opts.tol;     else tol = 1e-6; end
if isfield(opts,'r')     r = opts.r;         else r = 1; end

%% to solve the possion equation
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,  nx);
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt, 1);
lap = repmat(lap_t,1,nx) + repmat(lap_x,nt,1);

%% main iteration
nit = 0;
objs = zeros(nit,1);
conss = zeros(nit,1);
while nit<maxit 
    %% step A: update phi
    phi = mirt_dctn( Dt(rho/r-a) + Dx(mx/r-b) )./lap;
    phi(1,1) = 0;
    phi = mirt_idctn(phi);
    
    %% step B: update q
    % calculate a,b
    grad_t_phi = a;
    grad_x_phi = b;
        
    grad_t_phi(2:end-1,:) = Dt(phi);
    grad_x_phi(:,2:end-1) = Dx(phi);
    
    a = grad_t_phi + rho/r;
    b = grad_x_phi + mx /r;
    
    % projection        
    [a,b] = proj_K(a,b);
    a([1,end],:) = 0;
    b(:,[1,end]) = 0;
    
    %% step C: update mu
    drho = r*(grad_t_phi - a);
    dmx = r*(grad_x_phi - b);
    rho = rho + drho;
    mx = mx + dmx;
    
    %% update
    nit = nit+1;
    
    objs(nit) = objective(It(rho),Ix(mx));   
    conss(nit) = max(abs( Dt(rho) + Dx(mx) ),[],'all');
    
    res = sqrt( dt*dx*( norm(drho(:))^2+norm(dmx(:))^2 ) );
    if res < tol
        break
    end
end

outs = [];
outs.objs = objs(1:nit);
outs.conss = conss(1:nit);
outs.phi = phi;
outs.a = a;
outs.b = b;

%% ---------------------------------
    function A_t = Dt(A)
        A_t = ( A(2:end,:) - A(1:end-1,:) )/dt;
    end

    function A_x = Dx(A)
        A_x = ( A(:,2:end) - A(:,1:end-1) )/dx;
    end

    function A_it = It(A)
        A_it = ( A(1:end-1,:) + A(2:end,:) )/2;
    end

    function A_ix = Ix(A)
        A_ix = ( A(:,1:end-1) + A(:,2:end) )/2;
    end

    function obj = objective(rho,mx)
        ind = rho > 1e-8;
        obj = sum(mx(ind).^2./rho(ind))*dt*dx;
    end

end