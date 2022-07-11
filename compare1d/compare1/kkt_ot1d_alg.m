function [res_stat, res_feas] = kkt_ot1d_alg(rho,mx,phi,a,b)
res_feas = zeros(2,1);

[ntp,nx] = size(rho);
nt = ntp - 1;
% ntm = nt - 1; nxm = nx - 1; 
dt = 1/nt;    dx = 1/nx;    

%% stationarity

err_stat = {a(2:end-1,:) - Dt(phi), b(:,2:end-1) - Dx(phi)};
res_stat = cellfun(@(A) [max(abs(A),[],'all');norm(A(:))*sqrt(dt*dx)], ...
                        err_stat, 'UniformOutput', false); 
res_stat = cell2mat(res_stat);


%% primal feasibility
err_feas = Dt(rho) + Dx(mx);

res_feas(1) = max(abs(err_feas),[],'all');
res_feas(2) = norm(err_feas(:))*sqrt(dt*dx);

%%

    function DtA = Dt(A)
        DtA = (A(2:end,:) - A(1:end-1,:))/dt;
    end

    function DxA = Dx(A)
        DxA = (A(:,2:end) - A(:,1:end-1))/dx;
    end


end