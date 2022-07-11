function show_evolution(rho,filename)
if nargin == 1
    filename = 0; 
end
[ntp,nx] = size(rho);
rho0 = rho(1,:);
rho1 = rho(end,:);

t = linspace(0,1,ntp)';
x = linspace(1/nx,1-1/nx,nx);

%% pdf
fig = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile

plot3( repmat(t,1,nx)', repmat( x, ntp,1)', rho');
xlabel('t','FontSize',15,'FontWeight','bold'); 
ylabel('x','FontSize',15,'FontWeight','bold'); 
zlabel('\rho','FontSize',15,'FontWeight','bold');
axis square
if filename 
    exportgraphics(fig,[filename,'_plot3.pdf'],'BackgroundColor','white')
end

% contour(repmat( x, ntp,1), repmat(t,1,nx), rho,'LineWidth',1);
% xlabel('x','FontSize',15,'FontWeight','bold'); 
% ylabel('t','FontSize',15,'FontWeight','bold')
% if filename 
%     print(fig,'-dpdf',[filename,'_contour']); 
% end
% 
% mesh(repmat( x, ntp,1), repmat(t,1,nx), rho);
% xlabel('x','FontSize',15,'FontWeight','bold'); 
% ylabel('t','FontSize',15,'FontWeight','bold'); 
% zlabel('\rho','FontSize',15,'FontWeight','bold');
% if filename 
%     print(fig,'-dpdf',[filename,'_mesh']); 
% end

%% gif
im = cell(ntp,1);
figure('papersize',[5,4],'paperposition',[0,0,5,4]);
set(gcf,'color','w');
for i = 1:ntp
    clf;
    plot(x,rho0,'r','LineWidth',2);hold on;
    plot(x,rho1,'b','LineWidth',2);hold on;
    plot(x,rho(i,:),'LineWidth',2);
    pause(0.01);
    drawnow
    frame = getframe(gcf);
    im{i} = frame2im(frame);
    pause(0.01);
end


for idx = 1:ntp
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0);
    end
end


    
    


end

