clear
clc
close all

load('results\image_20.mat')

imshow(rho0,[]);
set(gcf,'unit','points','Position',[100 100 100 80]);
set(gca,'Position',[0.1,0.1,0.65,0.8]);
colorbar('Position',[0.76,0.1,0.05,0.8])
print('results\flb_rho0','-depsc');

imshow(rho1,[]);
set(gcf,'unit','points','Position',[100 100 100 80]);
set(gca,'Position',[0.1,0.1,0.65,0.8]);
colorbar('Position',[0.76,0.1,0.05,0.8])
print('results\flb_rho1','-depsc')

imshow(Qx,[]);
set(gcf,'unit','points','Position',[100 100 100 80]);
set(gca,'Position',[0.1,0.1,0.65,0.8]);
colorbar('Position',[0.76,0.1,0.05,0.8])
print('results\flb_Qx','-depsc')

imshow(Gx,[]);
set(gcf,'unit','points','Position',[100 100 100 80]);
set(gca,'Position',[0.1,0.1,0.65,0.8]);
colorbar('Position',[0.76,0.1,0.05,0.8])
print('results\flb_Gx','-depsc')

%% ot
num_frame = 6;

idx_frame = round(linspace(1,ntp,num_frame));
t_frame = (idx_frame-1)./(ntp-1);

%% ot
close all
for k = 2:num_frame
    idx = idx_frame(k);
    imshow(squeeze(rho_ot(idx,:,:)),[]);
    set(gcf,'unit','centimeters','position',[10 5 2 3])
    set(gca,'Position',[0.1,0.05,0.65,0.8]);
    colorbar('Position',[0.76,0.1,0.05,0.7])
    title(['t=',num2str(t_frame(k))]);
    fig = gcf;
    exportgraphics(fig,['results\flb_ot_shot',num2str(k),'.eps']);
end

%% mfp1
close all
for k = 2:num_frame
    idx = idx_frame(k);
    imshow(squeeze(rho_mfp1(idx,:,:)),[]);
    set(gcf,'unit','centimeters','position',[10 5 2 3])
    set(gca,'Position',[0.1,0.05,0.65,0.8]);
    colorbar('Position',[0.76,0.1,0.05,0.7])
    title(['t=',num2str(t_frame(k))]);
    fig = gcf;
    exportgraphics(fig,['results\flb_mfp1_shot',num2str(k),'.eps']);
end

%% mfp2
close all
for k = 2:num_frame
    idx = idx_frame(k);
    imshow(squeeze(rho_mfp2(idx,:,:)),[]);
    set(gcf,'unit','centimeters','position',[10 5 2 3])
    set(gca,'Position',[0.1,0.05,0.65,0.8]);
    colorbar('Position',[0.76,0.1,0.05,0.7])
    title(['t=',num2str(t_frame(k))]);
    fig = gcf;
    exportgraphics(fig,['results\flb_mfp2_shot',num2str(k),'.eps']);
end

%% mfp3
close all
for k = 2:num_frame
    idx = idx_frame(k);
    imshow(squeeze(rho_mfp3(idx,:,:)),[]);
    set(gcf,'unit','centimeters','position',[10 5 2 3])
    set(gca,'Position',[0.1,0.05,0.65,0.8]);
    colorbar('Position',[0.76,0.1,0.05,0.7])
    title(['t=',num2str(t_frame(k))]);
    fig = gcf;
    exportgraphics(fig,['results\flb_mfp3_shot',num2str(k),'.eps']);
end

%% mfg
close all
for k = 2:num_frame
    idx = idx_frame(k);
    imshow(squeeze(rho_mfg(idx,:,:)),[]);
    set(gcf,'unit','centimeters','position',[10 5 2 3])
    set(gca,'Position',[0.1,0.05,0.65,0.8]);
    colorbar('Position',[0.76,0.1,0.05,0.7])
    title(['t=',num2str(t_frame(k))]);
    fig = gcf;
    exportgraphics(fig,['results\flb_mfg_shot',num2str(k),'.eps']);
end
