function rho_new = inter1d_rho(rho0,rho1,rho)

% nt = size(rho,1)+1;
% nx = size(rho,2);
% 
% nt_new = nt;
% nx_new = nx*2;
% rho_new = ones(nt_new-1,nx_new);
% 
% rho_new(:,1:2:end) = rho;
% rho_new(:,2:2:end) = rho;

%% -----------------------------
rho0 = 0.5*(rho0(1:end-1)+rho0(2:end));
rho1 = 0.5*(rho1(1:end-1)+rho1(2:end));
nt = size(rho,1)+1;
nx = size(rho,2);

nt_new = nt*2;
nx_new = nx*2;
rho_new = zeros(nt_new-1,nx_new);

rho_new(2:2:end,1:2:end) = rho;
rho_new(2:2:end,2:2:end) = rho;

rho = [rho0;rho;rho1];
rho = 0.5*(rho(1:end-1,:)+rho(2:end,:));
rho_new(1:2:end,1:2:end) = rho;
rho_new(1:2:end,2:2:end) = rho;

end