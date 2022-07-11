function mx_new = restr2d_mx(mx)

mx = 0.5*( mx(1:2:end,:,:)+mx(2:2:end,:,:) );
mx = 0.5*( mx(:,:,1:2:end)+mx(:,:,2:2:end) );
mx_new = 0.25*( mx(:,1:2:end-2,:)+2*mx(:,2:2:end-1,:)+mx(:,3:2:end,:) );

end