function [rho, mx, my, outs] = ot2d_alg(rho0,rho1,opts)
%% parameters
% time and space parameters
if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
dt = 1/nt; ntp = nt+1;
if length(size(rho0))==2
    [nx,ny] = size(rho0);
    rho0 = reshape(rho0,1,nx,ny);
    rho1 = reshape(rho1,1,nx,ny);
else
    [~,nx,ny] = size(rho0);
end
dx = 1/nx; nxp = nx+1;
dy = 1/ny; nyp = ny+1;

% initial value of iteration
if isfield(opts,'rho') rho = opts.rho; else rho = ones(ntp,nx,ny); end
if isfield(opts,'mx')  mx = opts.mx;   else mx = ones(nt,nxp,ny);  end
if isfield(opts,'my')  my = opts.my;   else my = ones(nt,nx,nyp);  end
rho(1,:,:) = rho0; rho(end,:,:) = rho1;
mx(:,[1,end],:) = 0;
my(:,:,[1,end]) = 0;

if isfield(opts,'a')  a  = opts.a;  else a = ones(ntp,nx,ny); end
if isfield(opts,'bx') bx = opts.bx; else bx = ones(nt,nxp,ny); end
if isfield(opts,'by') by = opts.by; else by = ones(nt,nx,nyp); end
a([1,end],:,:) = 0;
bx(:,[1,end],:) = 0;
by(:,:,[1,end]) = 0;

% number of iteration, penalty parameters and stop criteria
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 500; end
if isfield(opts,'tol')   tol = opts.tol;     else tol = 1e-6; end
if isfield(opts,'r')     r = opts.r;         else r = 1; end

%% to solve the possion equation
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt,1,1);
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,nx,1);
lap_y = reshape((2 - 2*cos(pi*(0:ny-1)/ny))/dy/dy,1,1,ny);
lap = repmat(lap_t,1,nx,ny) + repmat(lap_x,nt,1,ny) + repmat(lap_y,nt,nx,1);

%% main iteration
nit = 0;
objs = zeros(nit,1);
conss = zeros(nit,1);
while nit<maxit 
    %% step A: update phi
    phi = mirt_dctn( Dt(rho/r-a) + Dx(mx/r-bx) + Dy(my/r-by) )./lap;
    phi(1,1,1) = 0;
    phi = mirt_idctn(phi);
    
    %% step B: update q
    % calculate a,b
    grad_t_phi = a;
    grad_x_phi = bx;
    grad_y_phi = by;
        
    grad_t_phi(2:end-1,:,:) = Dt(phi);
    grad_x_phi(:,2:end-1,:) = Dx(phi);
    grad_y_phi(:,:,2:end-1) = Dy(phi);
    
    a = grad_t_phi + rho/r;
    bx= grad_x_phi + mx /r;
    by= grad_y_phi + my /r;
    
    % projection        
    [a,bx,by] = proj_K(a,bx,by);
    a([1,end],:,:)  = 0;
    bx(:,[1,end],:) = 0;
    by(:,:,[1,end]) = 0;
    
    %% step C: update mu
    drho = r*(grad_t_phi - a);
    dmx = r*(grad_x_phi - bx);
    dmy = r*(grad_y_phi - by);
    rho = rho + drho;
    mx = mx + dmx;
    my = my + dmy;
    
    %% update
    nit = nit+1;
    objs(nit) = objective(It(rho),Ix(mx),Iy(my));   
    conss(nit) = max(abs( Dt(rho) + Dx(mx) + Dy(my)),[],'all');
    
    res = sqrt(dt*dx*dy*(norm(drho(:))^2+norm(dmx(:))^2+norm(dmy(:))^2));
    if res < tol
        break
    end
end

outs = [];
outs.objs = objs(1:nit);
outs.conss = conss(1:nit);
outs.phi = phi;
outs.a = a;
outs.b = {bx,by};

%% ---------------------------------
    function DtA = Dt(A)
        DtA = (A(2:end,:,:) - A(1:end-1,:,:))/dt;
    end

    function DxA = Dx(A)
        DxA = (A(:,2:end,:) - A(:,1:end-1,:))/dx;
    end

    function DyA = Dy(A)
        DyA = (A(:,:,2:end) - A(:,:,1:end-1))/dy;
    end

    function ItA = It(A)
        ItA = (A(1:end-1,:,:) + A(2:end,:,:))/2;
    end

    function IxA = Ix(A)
        IxA = (A(:,1:end-1,:) + A(:,2:end,:))/2;
    end

    function IyA = Iy(A)
        IyA = (A(:,:,1:end-1) + A(:,:,2:end))/2;
    end

    function obj = objective(rho,mx,my)
        ind = rho > 1e-8;
        obj = sum( (mx(ind).^2+my(ind).^2)./rho(ind) )*dt*dx*dy;
    end

end