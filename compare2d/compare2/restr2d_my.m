function my_new = restr2d_my(my)

my = 0.5*( my(1:2:end,:,:)+my(2:2:end,:,:) );
my = 0.5*( my(:,1:2:end,:)+my(:,2:2:end,:) );
my_new = 0.25*( my(:,:,1:2:end-2)+2*my(:,:,2:2:end-1)+my(:,:,3:2:end) );

end