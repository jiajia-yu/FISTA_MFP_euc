function [rho_k, mx_k, my_k, outs] = ot2d_gprox(rho0,rho1,opts)
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

% inital value of rho and mx
if isfield(opts,'rho') rho_k = opts.rho; else rho_k = ones(ntp,nx,ny); end
if isfield(opts,'mx')  mx_k = opts.mx;   else mx_k = ones(nt,nxp,ny);  end
if isfield(opts,'my')  my_k = opts.my;   else my_k = ones(nt,nx,nyp);  end
if isfield(opts,'phi') phi_k = opts.phi;   else phi_k = ones(nt,nx,ny);    end

rho_k(1,:,:) = rho0; rho_k(end,:,:) = rho1;
mx_k(:,[1,end],:) = 0;
my_k(:,:,[1,end]) = 0;

% number of iteration and stepsize
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 1e5;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end
if isfield(opts,'alpha') alpha = opts.alpha; else alpha = 1e-2; end

% penalty parameter
if isfield(opts,'tau1')  tau1 = opts.tau1;   else tau1 = 1;  end
if isfield(opts,'tau2')  tau2 = opts.tau2;   else tau2 = 1;  end

%% to solve the possion equation
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt,1,1);
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,nx,1);
lap_y = reshape((2 - 2*cos(pi*(0:ny-1)/ny))/dy/dy,1,1,ny);
lap = repmat(lap_t,1,nx,ny) + repmat(lap_x,nt,1,ny) + repmat(lap_y,nt,nx,1);

%% main iteration
nit = 0;
objs = zeros(maxit,1);
conss = zeros(maxit,1);
while nit < maxit 
    % update rho
    rho_kp1 = rho_k;
    polyb = -(tau1*Dt(phi_k) + rho_k(2:end-1,:,:));
    polyd = -tau1/2* It( (Ix(mx_k)).^2 + (Iy(my_k)).^2 );
    rho_kp1(2:end-1,:,:) = solve_cubic(1,polyb,0,polyd);
    drho = rho_kp1 - rho_k;
    
    % update mx,my
    mx_kp1 = mx_k;
    my_kp1 = my_k;
    rho_it = It(rho_kp1);
    rho_temp = Ix(rho_it);
    mx_kp1(:,2:end-1,:) = rho_temp ./ (rho_temp + tau1) .* ...
                        ( mx_k(:,2:end-1,:) + tau1*Dx(phi_k) );
    rho_temp = Iy(rho_it);
    my_kp1(:,:,2:end-1) = rho_temp ./ (rho_temp + tau1) .* ...
                        ( my_k(:,:,2:end-1) + tau1*Dy(phi_k) );
    dmx = mx_kp1 - mx_k;
    dmy = my_kp1 - my_k;
    
    % update rho_tilde, mx_tilde, my_tilde
    rho_tilde = drho + rho_kp1;
    mx_tilde = dmx + mx_kp1;
    my_tilde = dmy + my_kp1;
    
    % update phi
    phi_kp1 = tau2* (Dt(rho_tilde) + Dx(mx_tilde) + Dy(my_tilde));
    phi_kp1 = mirt_dctn(phi_kp1)./lap;
    phi_kp1(1,1) = 0;
    phi_kp1 = mirt_idctn(phi_kp1);
    phi_kp1 = phi_kp1 + phi_k;
    
    % update step
    nit = nit + 1;
    rho_k = rho_kp1;
    mx_k = mx_kp1;
    my_k = my_kp1;
    phi_k = phi_kp1;
    objs(nit) = objective( It(rho_k), Ix(mx_k), Iy(my_k) );
    conss(nit) = max(abs(Dt(rho_k) + Dx(mx_k) + Dy(my_k)), [], 'all');
    
    res = sqrt( dt*dx*dy*( norm(drho(:))^2+norm(dmx(:))^2+norm(dmy(:))^2 ) );
    if res < tol
        break
    end
end

outs.phi = phi_k;
outs.objs = objs(1:nit);
outs.conss = conss(1:nit);

%% 
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
