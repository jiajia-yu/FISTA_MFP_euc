function mx_new = inter2d_mx(mx)

% nt = size(mx,1);
% nx = size(mx,2)+1;
% ny = size(mx,3);
% 
% nt_new = nt;
% nx_new = nx*2;
% ny_new = ny*2;
% mx_new = zeros(nt_new,nx_new-1,ny_new);
% 
% mx_new(:,2:2:end,1:2:end) = mx;
% mx_new(:,2:2:end,2:2:end) = mx;
% 
% mx = cat(2,zeros(nt,1,ny),mx,zeros(nt,1,ny));
% mx = 0.5*(mx(:,1:end-1,:)+mx(:,2:end,:));
% mx_new(:,1:2:end,1:2:end) = mx;
% mx_new(:,1:2:end,2:2:end) = mx;

%% ------------------------------------------

nt = size(mx,1);
nx = size(mx,2)+1;
ny = size(mx,3);

nt_new = nt*2;
nx_new = nx*2;
ny_new = ny*2;
mx_new = zeros(nt_new,nx_new-1,ny_new);

mx_new(1:2:end,2:2:end,1:2:end) = mx;
mx_new(1:2:end,2:2:end,2:2:end) = mx;
mx_new(2:2:end,2:2:end,1:2:end) = mx;
mx_new(2:2:end,2:2:end,2:2:end) = mx;

mx = cat(2,zeros(nt,1,ny),mx,zeros(nt,1,ny));
mx = 0.5*(mx(:,1:end-1,:)+mx(:,2:end,:));
mx_new(1:2:end,1:2:end,1:2:end) = mx;
mx_new(1:2:end,1:2:end,2:2:end) = mx;
mx_new(2:2:end,1:2:end,1:2:end) = mx;
mx_new(2:2:end,1:2:end,2:2:end) = mx;

end


