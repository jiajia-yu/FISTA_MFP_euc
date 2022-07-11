function mx_new = restr1d_mx(mx)

mx_new = 0.5*( mx(1:2:end,:)+mx(2:2:end,:) );
mx_new = 0.25*( mx_new(:,1:2:end-2)+2*mx_new(:,2:2:end-1)+mx_new(:,3:2:end) );

end


