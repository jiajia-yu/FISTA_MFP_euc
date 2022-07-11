function rho_new = restr1d_rho(rho)

rho_new = 0.5*( rho(:,1:2:end)+rho(:,2:2:end) );
rho_new = 0.25*( rho_new(1:2:end-2,:)+2*rho_new(2:2:end-1,:)+rho_new(3:2:end,:) );

end