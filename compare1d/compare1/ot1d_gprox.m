function [rho_k, mx_k, outs] = ot1d_gprox(rho0,rho1,opts)
rho0 = Ix(rho0);
rho1 = Ix(rho1);
%% parameters
% time and space parameters
nx = length(rho0);
nxp = nx + 1;
dx = 1/nx;
if isfield(opts,'nt') nt = opts.nt; else nt = 20; end
ntp = nt+1;
dt = 1/nt;

% inital value of rho and mx
if isfield(opts,'rho') rho_k = opts.rho; else rho_k = ones(ntp,nx); end
if isfield(opts,'mx')  mx_k = opts.mx;   else mx_k = ones(nt,nxp);  end
if isfield(opts,'phi') phi_k = opts.phi;   else phi_k = ones(nt,nx); end

rho_k(1,:,:) = rho0; rho_k(end,:,:) = rho1;
mx_k(:,[1,end],:) = 0;

% number of iteration and stepsize
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 1e5;  end
if isfield(opts,'tol')   tol   = opts.tol;   else tol   = 1e-6; end
if isfield(opts,'alpha') alpha = opts.alpha; else alpha = 1e-2; end

% penalty parameter
if isfield(opts,'tau1')  tau1 = opts.tau1;   else tau1 = 1;end
if isfield(opts,'tau2')  tau2 = opts.tau2;   else tau2 = 1;        end

%% to solve the possion equation
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,  nx);
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt, 1);
lap = repmat(lap_t,1,nx) + repmat(lap_x,nt,1);

%% main iteration
nit = 0;
objs = zeros(maxit,1);
conss = zeros(maxit,1);
while nit < maxit 
    % update rho
    rho_kp1 = rho_k;
    rho_kp1(2:end-1,:) = solve_cubic(1, -(tau1*Dt(phi_k) + rho_k(2:end-1,:)),...
                                     0, -tau1/2 * It(Ix(mx_k)).^2);
    drho = rho_kp1 - rho_k;
                                 
    % update mx
    mx_kp1 = mx_k;
    rho_ixit = Ix(It(rho_kp1));
    mx_kp1(:,2:end-1) = rho_ixit ./ (rho_ixit + tau1) .* ...
                        ( mx_k(:,2:end-1) + tau1*Dx(phi_k) );
    dmx = mx_kp1 - mx_k;
                    
    % update rho_tilde, mx_tilde
    rho_tilde = drho + rho_kp1;
    mx_tilde = dmx + mx_kp1;
    
    % update phi
    phi_kp1 = tau2* (Dt(rho_tilde)+Dx(mx_tilde));
    phi_kp1 = mirt_dctn(phi_kp1)./lap;
    phi_kp1(1,1) = 0;
    phi_kp1 = mirt_idctn(phi_kp1);
    phi_kp1 = phi_kp1 + phi_k;
    
    % update step
    nit = nit + 1;
    rho_k = rho_kp1;
    mx_k = mx_kp1;
    phi_k = phi_kp1;
    objs(nit) = objective( It(rho_k), Ix(mx_k) );
    conss(nit) = max(abs(Dt(rho_k) + Dx(mx_k)), [], 'all');
    
    res = sqrt( dt*dx*( norm(drho(:))^2+norm(dmx(:))^2 ) );
    if res < tol
        break
    end
end
outs.phi = phi_kp1;
outs.objs = objs(1:nit);
outs.conss = conss(1:nit);

%% 
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
        obj = sum( mx(ind).^2./rho(ind) , 'all') * dt*dx;
    end

end
