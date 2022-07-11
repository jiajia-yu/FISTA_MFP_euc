function [res_stat, res_feas] = kkt_ot2d_alg(rho,mx,my,phi,a,b)
res_feas = zeros(2,1);

[ntp,nx,ny] = size(rho);
nt = ntp - 1;
% ntm = nt - 1; nxm = nx - 1; 
dt = 1/nt;    dx = 1/nx;    dy = 1/ny;  

%% stationarity
err_stat = {a(2:end-1,:,:) - Dt(phi), b{1}(:,2:end-1,:) - Dx(phi), b{2}(:,:,2:end-1) - Dy(phi)};
res_stat = cellfun(@(A) [max(abs(A),[],'all');norm(A(:))*sqrt(dt*dx*dy)], ...
                        err_stat, 'UniformOutput', false); 
res_stat = cell2mat(res_stat);


%% primal feasibility
err_feas = Dt(rho) + Dx(mx) +Dy(my);

res_feas(1) = max(abs(err_feas),[],'all');
res_feas(2) = norm(err_feas(:))*sqrt(dt*dx*dy);

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

end

