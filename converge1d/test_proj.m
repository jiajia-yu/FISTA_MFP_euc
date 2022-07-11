% test poisson and projection
function test_proj
clear
clc
%% exact solution
phi_exact  = @(t,x) cos(pi*t).*cos(2*pi*x);
phit_exact = @(t,x) -  pi*sin(pi*t).*cos(2*pi*x);
phix_exact = @(t,x) -2*pi*cos(pi*t).*sin(2*pi*x);
phi_lap_exact = @(t,x) -5*pi^2 * phi_exact(t,x);

rho_tilde_exact = @(t,x) 3*pi*sin(pi*t).*cos(2*pi*x) + 3*pi*exp(t-x);
mx_tilde_exact  = @(t,x)   pi*cos(pi*t).*sin(2*pi*x) + 3*pi*exp(t-x);
rhot_tilde_exact = @(t,x) 3*pi^2*cos(pi*t).*cos(2*pi*x) + 3*pi*exp(t-x);
mxx_tilde_exact  = @(t,x) 2*pi^2*cos(pi*t).*cos(2*pi*x) - 3*pi*exp(t-x);

rho_exact       = @(t,x) 2*pi*sin(pi*t).*cos(2*pi*x) + 3*pi*exp(t-x);
mx_exact        = @(t,x) - pi*cos(pi*t).*sin(2*pi*x) + 3*pi*exp(t-x);
rhot_exact      = @(t,x)  2*pi^2*cos(pi*t).*cos(2*pi*x) + 3*pi*exp(t-x);
mxx_exact       = @(t,x) -2*pi^2*cos(pi*t).*cos(2*pi*x) - 3*pi*exp(t-x);

%% mesh size
nt0 = 8;
nx0 = 8;
p_max = 8;
error_phi_maxnorms = zeros(p_max,3);
error_phi_l2norms = zeros(p_max,3);
error_rhomx_maxnorms = zeros(p_max,3);
error_rhomx_l2norms = zeros(p_max,3);
%% grid
for p = 1:p_max
    nt = nt0*2^p; ntp = nt+1; dt = 1/nt;
    nx = nx0*2^p; nxp = nx+1; dx = 1/nx;
    fprintf('\nnt=%d,nx=%d\n',nt,nx);

    lap_x = reshape((2 - 2*cos(pi*(0:nx-1)/nx))/dx/dx,1,  nx);
    lap_t = reshape((2 - 2*cos(pi*(0:nt-1)/nt))/dt/dt,nt, 1);
    lap = repmat(lap_t,1,nx) + repmat(lap_x,nt,1);

    t_r = linspace(0,1,ntp)';
    t_s = (t_r(1:end-1)+t_r(2:end))/2;
    x_r = linspace(0,1,nxp);
    x_s = (x_r(1:end-1)+x_r(2:end))/2;

    Trs = repmat(t_r,1,nx);
    Xrs = repmat(x_s,ntp,1);

    Tsr = repmat(t_s,1,nxp);
    Xsr = repmat(x_r,nt,1);

    Tss = repmat(t_s,1,nx);
    Xss = repmat(x_s,nt,1);

    % %----------------------------------------------------------------------
    % %% test exact solution
    % error = Dt(phi_exact(Trs,Xrs)) - phit_exact(Tss,Xss);
    % disp(max(abs(error(:))));
    % error = Dx(phi_exact(Tsr,Xsr)) - phix_exact(Tss,Xss);
    % disp(max(abs(error(:))));
    % error = Dt(phit_exact(Trs,Xrs)) + Dx(phix_exact(Tsr,Xsr)) - phi_lap_exact(Tss,Xss);
    % disp(max(abs(error(:))))
    % 
    % error = Dt(rho_tilde_exact(Trs,Xrs)) - rhot_tilde_exact(Tss,Xss);
    % disp(max(abs(error(:))));
    % error = Dx(mx_tilde_exact(Tsr,Xsr)) - mxx_tilde_exact(Tss,Xss);
    % disp(max(abs(error(:))));
    % 
    % error = -phi_lap_exact(Tss,Xss) - ( rhot_tilde_exact(Tss,Xss)+mxx_tilde_exact(Tss,Xss) );
    % disp(max(abs(error(:))));
    % 
    % error = rho_exact(Tss,Xss) - ( rho_tilde_exact(Tss,Xss)+phit_exact(Tss,Xss) );
    % disp(max(abs(error(:))));
    % error = mx_exact(Tss,Xss)  - ( mx_tilde_exact(Tss,Xss) +phix_exact(Tss,Xss) );
    % disp(max(abs(error(:))));
    % error = rho_exact(Tss,Xss) - ( rho_tilde_exact(Tss,Xss)+Dt(phi_exact(Trs,Xrs)) );
    % disp(max(abs(error(:))));
    % error = mx_exact(Tss,Xss)  - ( mx_tilde_exact(Tss,Xss) +Dx(phi_exact(Tsr,Xsr)) );
    % disp(max(abs(error(:))));
    % 
    % error = rhot_exact(Tss,Xss) + mxx_exact(Tss,Xss);
    % disp(max(abs(error(:))));
    % error = Dt(rho_exact(Trs,Xrs)) + Dx(mx_exact(Tsr,Xsr));
    % disp(max(abs(error(:))));
    % %-----------------------------------------------------------------------
    %% test projection

    % ----rho_t+m_x exact-----------------------------------------
    RHS = rhot_tilde_exact(Tss,Xss) + mxx_tilde_exact(Tss,Xss);
    phi = mirt_dctn( RHS )./lap;
    phi(1,1) = 0;
    phi = mirt_idctn(phi);

    error_phi = phi - phi_exact(Tss,Xss);
    error_phi_maxnorms(p,1) = max(abs(error_phi(:)));
    error_phi_l2norms(p,1) = sqrt(dt*dx)*norm(error_phi(:));
    fprintf('rho_t+m_x exact, solve phi\n');
    fprintf('phi error: max norm = %e, L2_norm = %e\n',...
                        error_phi_maxnorms(p,1), error_phi_l2norms(p,1) );


    % ----
    rho = rho_tilde_exact(Trs(2:end-1,:),Xrs(2:end-1,:)) + Dt(phi);
    mx  = mx_tilde_exact(Tsr(:,2:end-1),Xsr(:,2:end-1))  + Dx(phi);

    error_rho = rho - rho_exact(Trs(2:end-1,:),Xrs(2:end-1,:));
    error_mx  = mx  - mx_exact(Tsr(:,2:end-1),Xsr(:,2:end-1));
    fprintf('rho error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_rho(:))), sqrt(dt*dx)*norm(error_rho(:)) );
    fprintf('mx  error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_mx(:))), sqrt(dt*dx)*norm(error_mx(:)) );
    error_rhomx = [error_rho(:);error_mx(:)];
    error_rhomx_maxnorms(p,1) = max(abs(error_rhomx(:)));
    error_rhomx_l2norms(p,1) = sqrt(dt*dx)*norm(error_rhomx(:));
    fprintf('{rho,mx} error: max norm = %e, L2_norm = %e\n',...
                        error_rhomx_maxnorms(p,1), error_rhomx_l2norms(p,1) );


    % ----rho_t+m_x finite difference---------------
    RHS = Dt(rho_tilde_exact(Trs,Xrs)) + Dx(mx_tilde_exact(Tsr,Xsr));
    phi = mirt_dctn( RHS )./lap;
    phi(1,1) = 0;
    phi = mirt_idctn(phi);

    error_phi = phi - phi_exact(Tss,Xss);
    error_phi_maxnorms(p,2) = max(abs(error_phi(:)));
    error_phi_l2norms(p,2) = sqrt(dt*dx)*norm(error_phi(:),2);
    fprintf('rho_t+m_x finite difference, solve phi\n');
    fprintf('phi error: max norm = %e, L2_norm = %e\n',...
                        error_phi_maxnorms(p,2), error_phi_l2norms(p,2) );

    % ----
    rho = rho_tilde_exact(Trs(2:end-1,:),Xrs(2:end-1,:)) + Dt(phi);
    mx  = mx_tilde_exact(Tsr(:,2:end-1),Xsr(:,2:end-1))  + Dx(phi);

    error_rho = rho - rho_exact(Trs(2:end-1,:),Xrs(2:end-1,:));
    error_mx  = mx  - mx_exact(Tsr(:,2:end-1),Xsr(:,2:end-1));
    fprintf('rho error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_rho(:))), sqrt(dt*dx)*norm(error_rho(:),2) );
    fprintf('mx  error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_mx(:))), sqrt(dt*dx)*norm(error_mx(:),2) );
    error_rhomx = [error_rho(:);error_mx(:)];
    error_rhomx_maxnorms(p,2) = max(abs(error_rhomx(:)));
    error_rhomx_l2norms(p,2) = sqrt(dt*dx)*norm(error_rhomx(:),2);
    fprintf('{rho,mx} error: max norm = %e, L2_norm = %e\n',...
                        error_rhomx_maxnorms(p,2), error_rhomx_l2norms(p,2) );

    % ----rho, m, perturbed ---------------
    pert_rho_tilde = rho_tilde_exact(Trs,Xrs) + (randi(1,ntp,nx)-0.5)*sqrt(dt*dx);
    pert_mx_tilde  = mx_tilde_exact(Tsr,Xsr)  + (randi(1,nt,nxp)-0.5)*sqrt(dt*dx);
    RHS = Dt(pert_rho_tilde) + Dx(pert_mx_tilde) ;
    phi = mirt_dctn( RHS )./lap;
    phi(1,1) = 0;
    phi = mirt_idctn(phi);

    error_phi = phi - phi_exact(Tss,Xss);
    error_phi_maxnorms(p,3) = max(abs(error_phi(:)));
    error_phi_l2norms(p,3) = sqrt(dt*dx)*norm(error_phi(:),2);
    fprintf('rho_t+m_x finite difference, rho_tilde,mx_tilde perturbed, solve phi\n');
    fprintf('phi error: max norm = %e, L2_norm = %e\n',...
                        error_phi_maxnorms(p,2), error_phi_l2norms(p,2) );

    % ----
    rho = pert_rho_tilde(2:end-1,:) + Dt(phi);
    mx  = pert_mx_tilde(:,2:end-1)  + Dx(phi);

    error_rho = rho - rho_exact(Trs(2:end-1,:),Xrs(2:end-1,:));
    error_mx  = mx  - mx_exact(Tsr(:,2:end-1),Xsr(:,2:end-1));
    fprintf('rho error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_rho(:))), sqrt(dt*dx)*norm(error_rho(:),2) );
    fprintf('mx  error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_mx(:))), sqrt(dt*dx)*norm(error_mx(:),2) );
    error_rhomx = [error_rho(:);error_mx(:)];
    error_rhomx_maxnorms(p,3) = max(abs(error_rhomx(:)));
    error_rhomx_l2norms(p,3) = sqrt(dt*dx)*norm(error_rhomx(:),2);
    fprintf('{rho,mx} error: max norm = %e, L2_norm = %e\n',...
                        error_rhomx_maxnorms(p,2), error_rhomx_l2norms(p,2) );
                    
                    
                    
end
order_phi_maxnorms = log(error_phi_maxnorms(1:end-1,:)./error_phi_maxnorms(2:end,:))/log(2);
order_phi_l2norms = log(error_phi_l2norms(1:end-1,:)./error_phi_l2norms(2:end,:))/log(2);
order_rhomx_maxnorms = log(error_rhomx_maxnorms(1:end-1,:)./error_rhomx_maxnorms(2:end,:))/log(2);
order_rhomx_l2norms = log(error_rhomx_l2norms(1:end-1,:)./error_rhomx_l2norms(2:end,:))/log(2);

fprintf('rhot_tilde+mxx_tilde exact\n')
fprintf('phi error max norm order');
disp(order_phi_maxnorms(:,1)');
fprintf('phi error L2 norm order');
disp(order_phi_l2norms(:,1)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,1)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,1)');

fprintf('rhot_tilde+mxx_tilde finite difference\n')
fprintf('phi error max norm order');
disp(order_phi_maxnorms(:,2)');
fprintf('phi error L2 norm order');
disp(order_phi_l2norms(:,2)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,2)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,2)');

fprintf('rho_tilde, mx_tilde randomly perturbed\n')
fprintf('phi error max norm order');
disp(order_phi_maxnorms(:,3)');
fprintf('phi error L2 norm order');
disp(order_phi_l2norms(:,3)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,3)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,3)');

%%
function DtA = Dt(A)
    DtA = (A(2:end,:) - A(1:end-1,:))/dt;
end

function DxA = Dx(A)
    DxA = (A(:,2:end) - A(:,1:end-1))/dx;
end

% function ItA = It(A)
%     ItA = (A(1:end-1,:) + A(2:end,:))/2;
% end
% 
% function IxA = Ix(A)
%     IxA = (A(:,1:end-1) + A(:,2:end))/2;
% end

end
