function [res_stat, res_feas] = kkt_ot1d(rho,mx)
res_feas = zeros(2,1);

[ntp,nx] = size(rho);
nt = ntp - 1;
% ntm = nt - 1; nxm = nx - 1; 
dt = 1/nt;    dx = 1/nx;    

lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,  nx);
lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt, 1);
lap = repmat(lap_t,1,nx) + repmat(lap_x,nt,1);

%% stationarity
[grad_rho,grad_mx] = compGrad(rho,mx);

phi = mirt_dctn( Dtadj(grad_rho)+Dxadj(grad_mx) )./lap;
phi(1,1,1) = 0;
phi = mirt_idctn(phi);

err_stat = {grad_rho - Dt(phi), grad_mx - Dx(phi)};
res_stat = cellfun(@(A) [max(abs(A),[],'all');norm(A(:))*sqrt(dt*dx)], ...
                        err_stat, 'UniformOutput', false); 
res_stat = cell2mat(res_stat);


%% primal feasibility
err_feas = Dt(rho) + Dx(mx);

res_feas(1) = max(abs(err_feas),[],'all');
res_feas(2) = norm(err_feas(:))*sqrt(dt*dx);

%%

    function [grad_rho,grad_mx] = compGrad(rho,mx)
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

    function DtadjA = Dtadj(A)
        DtadjA = [-A(1,:); ...
                   A(1:end-1,:)-A(2:end,:); ...
                   A(end,:)] /dt;
    end

    function DxadjA = Dxadj(A)
        DxadjA = cat(2, -A(:,1), ...
                         A(:,1:end-1)-A(:,2:end), ...
                         A(:,end)) /dx;
    end

    function DtA = Dt(A)
        DtA = (A(2:end,:) - A(1:end-1,:))/dt;
    end

    function DxA = Dx(A)
        DxA = (A(:,2:end) - A(:,1:end-1))/dx;
    end

    function ItA = It(A)
        ItA = (A(1:end-1,:) + A(2:end,:))/2;
    end

    function IxA = Ix(A)
        IxA = (A(:,1:end-1) + A(:,2:end))/2;
    end


end