% test poisson and projection
function test_gd
clear
clc
%% exact solution
rho_tilde_exact = @(t,x) (1-t).*(x-0.5) + 1;
mx_tilde_exact = @(t,x) x.*(1-x);
eta = 0.1;

% L_exact = @(rho,mx) sum( (mx.^2)./rho,'all')/2;
grad_rho_exact = @(rho,mx) -(mx.^2)./(rho.^2)/2;
grad_mx_exact  = @(rho,mx)  mx./rho;
rho_exact = @(t,x) rho_tilde_exact(t,x) ...
            - eta*grad_rho_exact( rho_tilde_exact(t,x),mx_tilde_exact(t,x) );
mx_exact  = @(t,x) mx_tilde_exact(t,x) ...
            - eta*grad_mx_exact( rho_tilde_exact(t,x),mx_tilde_exact(t,x) );

%% mesh size
nt0 = 8;
nx0 = 64;
p_max = 8;
error_grad_maxnorms = zeros(p_max,3);
error_grad_l2norms = zeros(p_max,3);
error_rhomx_maxnorms = zeros(p_max,3);
error_rhomx_l2norms = zeros(p_max,3);
%% grid
for p = 1:p_max
    nt = nt0*2^p; ntp = nt+1; dt = 1/nt;
    nx = nx0*2^p; nxp = nx+1; dx = 1/nx;
    fprintf('\nnt=%d,nx=%d\n',nt,nx);

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

    % ----rho(t,x), m(t,x) exact on all grid points---------------------------
    rho_tilde = rho_tilde_exact(Tss,Xss);
    mx_tilde = mx_tilde_exact(Tss,Xss);
    grad_rho = grad_rho_exact(rho_tilde,mx_tilde);
    grad_mx  = grad_mx_exact(rho_tilde,mx_tilde);
    
    error_gradrho = grad_rho ...
      - grad_rho_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_gradmx  = grad_mx  ...
      -  grad_mx_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_grad = [error_gradrho(:);error_gradmx(:)];
    error_grad_maxnorms(p,1) = max(abs(error_grad(:)));
    error_grad_l2norms(p,1) = sqrt(dt*dx)*norm(error_grad(:));
    fprintf('rho(t,x), m(t,x) exact on all grid points\n');
    fprintf('grad error: max norm = %e, L2_norm = %e\n',...
                        error_grad_maxnorms(p,1), error_grad_l2norms(p,1) );


    % ----
    rho = rho_tilde_exact(Trs(2:end-1,:),Xrs(2:end-1,:)) - eta*It( grad_rho );
    mx  = mx_tilde_exact(Tsr(:,2:end-1),Xsr(:,2:end-1))  - eta*Ix( grad_mx );

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


    % ----rho(t,x), m(t,x) exact on rs-grid, sr-grid---------------
    rho_tilde = rho_tilde_exact(Trs,Xrs) ;
    mx_tilde =  mx_tilde_exact(Tsr,Xsr) ;
    grad_rho = grad_rho_exact( It(rho_tilde),Ix(mx_tilde) );
    grad_mx  = grad_mx_exact( It(rho_tilde),Ix(mx_tilde) );
    
    error_gradrho = grad_rho ...
      - grad_rho_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_gradmx  = grad_mx  ...
      -  grad_mx_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_grad = [error_gradrho(:);error_gradmx(:)];
    error_grad_maxnorms(p,2) = max(abs(error_grad(:)));
    error_grad_l2norms(p,2) = sqrt(dt*dx)*norm(error_grad(:));
    fprintf('rho(t,x), m(t,x) exact on rs-grid, sr-grid\n');
    fprintf('grad error: max norm = %e, L2_norm = %e\n',...
                        error_grad_maxnorms(p,1), error_grad_l2norms(p,1) );


    % ----
    rho = rho_tilde(2:end-1,:) - eta*It( grad_rho );
    mx  = mx_tilde(:,2:end-1)  - eta*Ix( grad_mx );

    error_rho = rho - rho_exact(Trs(2:end-1,:),Xrs(2:end-1,:));
    error_mx  = mx  - mx_exact(Tsr(:,2:end-1),Xsr(:,2:end-1));
    fprintf('rho error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_rho(:))), sqrt(dt*dx)*norm(error_rho(:)) );
    fprintf('mx  error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_mx(:))), sqrt(dt*dx)*norm(error_mx(:)) );
    error_rhomx = [error_rho(:);error_mx(:)];
    error_rhomx_maxnorms(p,2) = max(abs(error_rhomx(:)));
    error_rhomx_l2norms(p,2) = sqrt(dt*dx)*norm(error_rhomx(:));
    fprintf('{rho,mx} error: max norm = %e, L2_norm = %e\n',...
                        error_rhomx_maxnorms(p,1), error_rhomx_l2norms(p,1) );

    % ----rho(t,x), m(t,x) perturbed on rs-grid, sr-grid---------------
    rho_tilde = rho_tilde_exact(Trs,Xrs)+(rand(ntp,nx)-0.5)*(dt+dx);
    mx_tilde =  mx_tilde_exact(Tsr,Xsr)+(rand(nt,nxp)-0.5)*(dt+dx);
    grad_rho = grad_rho_exact( It(rho_tilde),Ix(mx_tilde) );
    grad_mx  = grad_mx_exact ( It(rho_tilde),Ix(mx_tilde) );
    
    error_gradrho = grad_rho ...
      - grad_rho_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_gradmx  = grad_mx  ...
      -  grad_mx_exact( rho_tilde_exact(Tss,Xss),mx_tilde_exact(Tss,Xss) );
    error_grad = [error_gradrho(:);error_gradmx(:)];
    error_grad_maxnorms(p,3) = max(abs(error_grad(:)));
    error_grad_l2norms(p,3) = sqrt(dt*dx)*norm(error_grad(:));
    fprintf('rho(t,x), m(t,x) perturbed on rs-grid, sr-grid\n');
    fprintf('grad error: max norm = %e, L2_norm = %e\n',...
                        error_grad_maxnorms(p,1), error_grad_l2norms(p,1) );

    
    % ----
    rho = rho_tilde(2:end-1,:) - eta*It( grad_rho );
    mx  = mx_tilde(:,2:end-1)  - eta*Ix( grad_mx );

    error_rho = rho - rho_exact(Trs(2:end-1,:),Xrs(2:end-1,:));
    error_mx  = mx  - mx_exact(Tsr(:,2:end-1),Xsr(:,2:end-1));
    fprintf('rho error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_rho(:))), sqrt(dt*dx)*norm(error_rho(:)) );
    fprintf('mx  error: max norm = %e, L2_norm = %e\n',...
                        max(abs(error_mx(:))), sqrt(dt*dx)*norm(error_mx(:)) );
    error_rhomx = [error_rho(:);error_mx(:)];
    error_rhomx_maxnorms(p,3) = max(abs(error_rhomx(:)));
    error_rhomx_l2norms(p,3) = sqrt(dt*dx)*norm(error_rhomx(:));
    fprintf('{rho,mx} error: max norm = %e, L2_norm = %e\n',...
                        error_rhomx_maxnorms(p,1), error_rhomx_l2norms(p,1) );    
                    
                    
end
order_grad_maxnorms = log(error_grad_maxnorms(1:end-1,:)./error_grad_maxnorms(2:end,:))/log(2);
order_grad_l2norms = log(error_grad_l2norms(1:end-1,:)./error_grad_l2norms(2:end,:))/log(2);
order_rhomx_maxnorms = log(error_rhomx_maxnorms(1:end-1,:)./error_rhomx_maxnorms(2:end,:))/log(2);
order_rhomx_l2norms = log(error_rhomx_l2norms(1:end-1,:)./error_rhomx_l2norms(2:end,:))/log(2);

fprintf('rho(t,x), m(t,x) exact on all grid points\n')
fprintf('grad error max norm');
disp(error_grad_maxnorms(:,1)');
fprintf('grad error L2 norm');
disp(error_grad_l2norms(:,1)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,1)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,1)');

fprintf('rho(t,x), m(t,x) exact on rs-grid, sr-grid\n')
fprintf('grad error max norm order');
disp(order_grad_maxnorms(:,2)');
fprintf('grad error L2 norm order');
disp(order_grad_l2norms(:,2)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,2)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,2)');

fprintf('rho(t,x), m(t,x) perturbed on rs-grid, sr-grid\n')
fprintf('grad error max norm order');
disp(order_grad_maxnorms(:,3)');
fprintf('grad error L2 norm order');
disp(order_grad_l2norms(:,3)');
fprintf('rho,mx error max norm order');
disp(order_rhomx_maxnorms(:,3)');
fprintf('rho,mx error L2 norm order');
disp(order_rhomx_l2norms(:,3)');

%%
% function DtA = Dt(A)
%     DtA = (A(2:end,:) - A(1:end-1,:))/dt;
% end
% 
% function DxA = Dx(A)
%     DxA = (A(:,2:end) - A(:,1:end-1))/dx;
% end

function ItA = It(A)
    ItA = (A(1:end-1,:) + A(2:end,:))/2;
end

function IxA = Ix(A)
    IxA = (A(:,1:end-1) + A(:,2:end))/2;
end

end
