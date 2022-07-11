function my_new = inter2d_my(my)

% nt = size(my,1);
% nx = size(my,2);
% ny = size(my,3)+1;
% 
% nt_new = nt;
% nx_new = nx*2;
% ny_new = ny*2;
% my_new = zeros(nt_new,nx_new,ny_new-1);
% 
% my_new(:,1:2:end,2:2:end) = my;
% my_new(:,2:2:end,2:2:end) = my;
% 
% my = cat(3,zeros(nt,nx,1),my,zeros(nt,nx,1));
% my = 0.5*(my(:,:,1:end-1)+my(:,:,2:end));
% my_new(:,1:2:end,1:2:end) = my;
% my_new(:,2:2:end,1:2:end) = my;

%% -----------------------------------

nt = size(my,1);
nx = size(my,2);
ny = size(my,3)+1;

nt_new = nt*2;
nx_new = nx*2;
ny_new = ny*2;
my_new = zeros(nt_new,nx_new,ny_new-1);

my_new(1:2:end,1:2:end,2:2:end) = my;
my_new(1:2:end,2:2:end,2:2:end) = my;
my_new(2:2:end,1:2:end,2:2:end) = my;
my_new(2:2:end,2:2:end,2:2:end) = my;

my = cat(3,zeros(nt,nx,1),my,zeros(nt,nx,1));
my = 0.5*(my(:,:,1:end-1)+my(:,:,2:end));
my_new(1:2:end,1:2:end,1:2:end) = my;
my_new(1:2:end,2:2:end,1:2:end) = my;
my_new(2:2:end,1:2:end,1:2:end) = my;
my_new(2:2:end,2:2:end,1:2:end) = my;

end


