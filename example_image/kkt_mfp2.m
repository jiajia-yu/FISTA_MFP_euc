function [res_stat, res_feas] = kkt_mfp2(rho,mx,my,opts)
res_feas = zeros(2,1);

[ntp,nx,ny] = size(rho);
nt = ntp - 1;
% ntm = nt - 1; nxm = nx - 1; 
dt = 1/nt;    dx = 1/nx;    dy = 1/ny;  

lambda_L = opts.lambda_L ;
lambda_E = opts.lambda_E;
lambda_P = opts.lambda_P;
Qx = reshape(opts.Qx,1,nx,ny);
Qx = repmat(Qx,nt,1,1);

lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt,1,1);
lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,nx,1);
lap_y = reshape((2 - 2*cos(pi*(0:ny-1)/ny))/dy/dy,1,1,ny);
lap = repmat(lap_t,1,nx,ny) + repmat(lap_x,nt,1,ny) + repmat(lap_y,nt,nx,1);

%% stationarity
[grad_rho,grad_mx,grad_my] = compGrad(rho,mx,my);

phi = mirt_dctn( Dtadj(grad_rho)+Dxadj(grad_mx)+Dyadj(grad_my) )./lap;
phi(1,1,1) = 0;
phi = mirt_idctn(phi);

err_stat = {grad_rho - Dt(phi), grad_mx - Dx(phi), grad_my - Dy(phi)};
res_stat = cellfun(@(A) [max(abs(A),[],'all');norm(A(:))*sqrt(dt*dx*dy)], ...
                        err_stat, 'UniformOutput', false); 
res_stat = cell2mat(res_stat);


%% primal feasibility
err_feas = Dt(rho) + Dx(mx) +Dy(my);

res_feas(1) = max(abs(err_feas),[],'all');
res_feas(2) = norm(err_feas(:))*sqrt(dt*dx*dy);

%%

    function [grad_rho,grad_mx,grad_my] = compGrad(rho,mx,my)
        rho = It(rho);
        mx = Ix(mx);
        my = Iy(my);
%         ind = rho > 1e-8;
        ind = rho > 0;
%         ind = rho ~=0;
        
        grad_rho = zeros(size(rho));
        grad_mx = zeros(size(mx));
        grad_my = zeros(size(my));
        
        grad_rho(ind) = - lambda_L*(mx(ind).^2+my(ind).^2)./rho(ind).^2/2 ...
            + lambda_E*((rho(ind))) + lambda_P*Qx(ind);
        grad_mx(ind) = lambda_L*mx(ind)./rho(ind);
        grad_my(ind) = lambda_L*my(ind)./rho(ind);
        
        grad_rho = It(grad_rho);
        grad_mx = Ix(grad_mx);
        grad_my = Iy(grad_my);
    end

    function DtadjA = Dtadj(A)
        DtadjA = cat(1, -A(1,:,:), ...
                         A(1:end-1,:,:)-A(2:end,:,:), ...
                         A(end,:,:)) /dt;
    end

    function DxadjA = Dxadj(A)
        DxadjA = cat(2, -A(:,1,:), ...
                         A(:,1:end-1,:)-A(:,2:end,:), ...
                         A(:,end,:)) /dx;
    end

    function DyadjA = Dyadj(A)
        DyadjA = cat(3, -A(:,:,1), ...
                         A(:,:,1:end-1)-A(:,:,2:end), ...
                         A(:,:,end)) /dy;
    end

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

end