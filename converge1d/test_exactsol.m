% test exact solution
p_max = 3;
ress = zeros(6,p_max);
W2sq = zeros(1,p_max);
W2sqerror = zeros(1,p_max);
for p = 1:p_max

nt = 16*2^p; ntp = nt+1; 
dt = 1/nt;

nx = 128*2^p; nxp = nx+1; 
dx = 1/nx;

%% rho_exact
t = linspace(0,1,ntp)';
x = linspace(dx/2,1-dx/2,nx);
T = repmat(t,1,nx);
X = repmat(x,ntp,1);

sqrtterm = sqrt(2*T.*X+(T/2-1).^2);
rho_exact = (sqrtterm + T - 1)./(T.*sqrtterm);
rho_exact(1,:) = x + 0.5;
figure(1);mesh(X,T,rho_exact);title('\rho');xlabel('x');ylabel('t')

%% mx_exact
t = linspace(dt/2,1-dt/2,nt)';
x = linspace(0,1,nxp);
T = repmat(t,1,nxp);
X = repmat(x,nt,1);

sqrtterm = sqrt(2*T.*X+(T/2-1).^2);
Tsq = T.^2;
Tcub = T.^3;
mx_exact = X./Tsq + (T-3)./(2*Tcub).*sqrtterm ...
         - (T-1).*(Tsq-4)./(8*Tcub)./sqrtterm - (3*T-4)./(2*Tcub);
figure(2);mesh(X,T,mx_exact);title('m');xlabel('x');ylabel('t')

%% kkt
[res_stat,res_feas] = kkt_ot1d(rho_exact,mx_exact);
ress(:,p) = [res_stat(:);res_feas(:)];

%% energy
rho_exact = ( rho_exact(1:end-1,:)+rho_exact(2:end,:) )/2;
mx_exact = ( mx_exact(:,1:end-1)+mx_exact(:,2:end) )/2;
W2sq(p) = sum( mx_exact.^2./rho_exact, 'all' )*dt*dx;

end
W2sqerror = abs(W2sq-1/120);

disp( ress' );
disp( (log(ress(:,1:end-1)./ress(:,2:end))./log(2))' );
disp( W2sqerror );
disp( log(W2sqerror(1:end-1)./W2sqerror(2:end))./log(2) );

