function rho_new = restr2d_rho(rho)

rho = 0.5*( rho(:,1:2:end,:)+rho(:,2:2:end,:) );
rho = 0.5*( rho(:,:,1:2:end)+rho(:,:,2:2:end) );
rho_new = 0.25*( rho(1:2:end-2,:,:)+2*rho(2:2:end-1,:,:)+rho(3:2:end,:,:) );

end