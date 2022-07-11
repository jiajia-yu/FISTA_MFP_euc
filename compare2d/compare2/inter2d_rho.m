function rho_new = inter2d_rho(rho0,rho1,rho)

% nt = size(rho,1)+1;
% nx = size(rho,2);
% ny = size(rho,3);
% 
% nt_new = nt;
% nx_new = nx*2;
% ny_new = ny*2;
% rho_new = ones(nt_new-1,nx_new,ny_new);
% 
% rho_new(:,1:2:end,1:2:end) = rho;
% rho_new(:,1:2:end,2:2:end) = rho;
% rho_new(:,2:2:end,1:2:end) = rho;
% rho_new(:,2:2:end,2:2:end) = rho;

%% ----------------------------------
if length(size(rho0)) == 2
    [nxp,nyp] = size(rho0);
    rho0 = reshape(rho0,1,nxp,nyp);
    rho1 = reshape(rho1,1,nxp,nyp);    
end

rho0 = 0.25*(rho0(:,1:end-1,1:end-1) + rho0(:,1:end-1,2:end) ...
           + rho0(:,2:end,1:end-1) + rho0(:,2:end,2:end));
rho1 = 0.25*(rho1(:,1:end-1,1:end-1) + rho1(:,1:end-1,2:end) ...
           + rho1(:,2:end,1:end-1) + rho1(:,2:end,2:end));

nt = size(rho,1)+1;
nx = size(rho,2);
ny = size(rho,3);

nt_new = nt*2;
nx_new = nx*2;
ny_new = ny*2;
rho_new = ones(nt_new-1,nx_new,ny_new);

rho_new(2:2:end,1:2:end,1:2:end) = rho;
rho_new(2:2:end,1:2:end,2:2:end) = rho;
rho_new(2:2:end,2:2:end,1:2:end) = rho;
rho_new(2:2:end,2:2:end,2:2:end) = rho;

rho = cat(1,rho0,rho,rho1);
rho = 0.5*(rho(1:end-1,:,:)+rho(2:end,:,:));
rho_new(1:2:end,1:2:end,1:2:end) = rho;
rho_new(1:2:end,1:2:end,2:2:end) = rho;
rho_new(1:2:end,2:2:end,1:2:end) = rho;
rho_new(1:2:end,2:2:end,2:2:end) = rho;

end