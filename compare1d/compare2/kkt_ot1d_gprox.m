function [res_stat, res_feas] = kkt_ot1d_gprox(rho,mx,phi)
res_feas = zeros(2,1);

[ntp,nx] = size(rho);
nt = ntp - 1;
% ntm = nt - 1; nxm = nx - 1; 
dt = 1/nt;    dx = 1/nx;    

%% stationarity
[grad_rho,grad_mx] = compGrad(rho,mx);

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