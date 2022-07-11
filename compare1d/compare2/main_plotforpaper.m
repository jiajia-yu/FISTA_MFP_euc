%%
clear
close all
clc

load('C:\Users\JIAJI\Dropbox\Jiajia_finished_mfg_euc\codes\_finalversion\compare1d\results\gaussian1d_compall.mat')
x = linspace(1/nx,1-1/nx,nx);
T = repmat(t(1:2:end),1,nx)';
X = repmat(x,length(t(1:2:end)),1)';

fig = tiledlayout(2,3,'TileSpacing','Compact','Padding','none');
% fig.Units = 'inches';
% fig.OuterPosition = [0 0 1.5 1.5];
set(gcf,'unit','points','Position',[0 0 350 240]);
% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize',[12 15])

nexttile
plot3(T, X, rho_fista(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

nexttile
plot3(T, X, rho_alg(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

nexttile
plot3(T, X, rho_gprox(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

nexttile
plot3(T, X, rho_ml(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

nexttile
plot3(T, X, rho_mg5(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

nexttile
plot3(T, X, rho_mg10(1:2:end,:)','Color',[0 0.4470 0.7410],'LineWidth',1)
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square

% set(gcf,'PaperUnits','centimeters');
% set(gcf, 'PaperSize',[12 15])
% exportgraphics(fig,'test.jpg','Resolution',300,'BackgroundColor','white')
exportgraphics(fig,'../results/gaussian1d_all_plot3.eps','Resolution',600)

