clear
clc
close all

fig = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
set(gcf,'unit','points','Position',[0 0 350 120]);

load('results\maze1.mat')
nexttile
eg_illustration(max(rho0,rho1),Qx);

load('results\maze2.mat')
nexttile
eg_illustration(max(rho0,rho1),Qx);

load('results\maze3.mat')
nexttile
eg_illustration(max(rho0,rho1),Qx);

exportgraphics(fig,'results\obstacle_egs.eps','BackgroundColor','none')

%%
clear
clc
close all

num_frame = 5;
vec_dens = 25;
vec_leng = 1;
tol = 1e-1;

load('results\maze1.mat')
mx_mfp = (mx_mfp(:,1:end-1,:) + mx_mfp(:,2:end,:) )/2;
my_mfp = (my_mfp(:,:,1:end-1) + my_mfp(:,:,2:end) )/2;
mx_mfp(mx_mfp < tol) = 0;
my_mfp(my_mfp < tol) = 0;

idx_frame = round(linspace(1,ntp,num_frame+1));
idx_frame = (idx_frame(1:end-1)+idx_frame(2:end))/2;
t_frame = (idx_frame-1)./(ntp-1);
X = repmat( (1:vec_dens:nx) ,length(1:vec_dens:ny),1 );
Y = repmat( (1:vec_dens:ny)',1, length(1:vec_dens:nx));


fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 130 130]);
nexttile
eg_illustration(max(rho0,rho1),Qx);
exportgraphics(fig,'results\obs1.eps','BackgroundColor','none');

close all
fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 80 90]);
nexttile
for k = 1:num_frame
    idx = idx_frame(k);
    eg_illustration(squeeze(rho_mfp(idx,:,:)),Qx);
    hold on

    quiver(X,Y,squeeze(my_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
               squeeze(mx_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
            vec_leng,'r','LineWidth',1);
        
    title(['t=',num2str(t_frame(k))]);
    exportgraphics(fig,['results\obs1_shot',num2str(k),'.eps'],'BackgroundColor','none');
end


%%
clear
clc
close all

num_frame = 5;
vec_dens = 25;
vec_leng = 1;
tol = 1e-1;

load('results\maze2.mat')
mx_mfp = (mx_mfp(:,1:end-1,:) + mx_mfp(:,2:end,:) )/2;
my_mfp = (my_mfp(:,:,1:end-1) + my_mfp(:,:,2:end) )/2;
mx_mfp(mx_mfp < tol) = 0;
my_mfp(my_mfp < tol) = 0;

idx_frame = round(linspace(1,ntp,num_frame+1));
idx_frame = (idx_frame(1:end-1)+idx_frame(2:end))/2;
t_frame = (idx_frame-1)./(ntp-1);
X = repmat( (1:vec_dens:nx) ,length(1:vec_dens:ny),1 );
Y = repmat( (1:vec_dens:ny)',1, length(1:vec_dens:nx));


fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 130 130]);
nexttile
eg_illustration(max(rho0,rho1),Qx);
exportgraphics(fig,'results\obs2.eps','BackgroundColor','none');

close all
fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 80 90]);
nexttile
for k = 1:num_frame
    idx = idx_frame(k);
    eg_illustration(squeeze(rho_mfp(idx,:,:)),Qx);
    hold on

    quiver(X,Y,squeeze(my_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
               squeeze(mx_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
            vec_leng,'r','LineWidth',1);
        
    title(['t=',num2str(t_frame(k))]);
    exportgraphics(fig,['results\obs2_shot',num2str(k),'.eps'],'BackgroundColor','none');
end

%%
clear
clc
close all

num_frame = 5;
vec_dens = 25;
vec_leng = 1;
tol = 1e-1;

load('results\maze3.mat')
mx_mfp = (mx_mfp(:,1:end-1,:) + mx_mfp(:,2:end,:) )/2;
my_mfp = (my_mfp(:,:,1:end-1) + my_mfp(:,:,2:end) )/2;
mx_mfp(mx_mfp < tol) = 0;
my_mfp(my_mfp < tol) = 0;

idx_frame = round(linspace(1,ntp,num_frame+1));
idx_frame = (idx_frame(1:end-1)+idx_frame(2:end))/2;
t_frame = (idx_frame-1)./(ntp-1);
X = repmat( (1:vec_dens:nx) ,length(1:vec_dens:ny),1 );
Y = repmat( (1:vec_dens:ny)',1, length(1:vec_dens:nx));


fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 130 130]);
nexttile
eg_illustration(max(rho0,rho1),Qx);
exportgraphics(fig,'results\obs3.eps','BackgroundColor','none');

close all
fig = tiledlayout(1,1,'TileSpacing','none','Padding','none');
set(gcf,'unit','points','Position',[0 0 80 90]);
nexttile
for k = 1:num_frame
    idx = idx_frame(k);
    eg_illustration(squeeze(rho_mfp(idx,:,:)),Qx);
    hold on

    quiver(X,Y,squeeze(my_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
               squeeze(mx_mfp(idx,1:vec_dens:end,1:vec_dens:end)),...
            vec_leng,'r','LineWidth',1);
        
    title(['t=',num2str(t_frame(k))]);
    exportgraphics(fig,['results\obs3_shot',num2str(k),'.eps'],'BackgroundColor','none');
end

