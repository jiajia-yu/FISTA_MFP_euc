function show_evolution(rho,showfunc,filename)
if nargin < 3 
    filename = 0;
else
    filename = [filename,'_',showfunc];
end
[ntp,nx,ny] = size(rho);
zmin = min(rho,[],'all');
zmax = max(rho,[],'all');
x = linspace(1/nx,1-1/nx,nx);
y = linspace(1/ny,1-1/ny,ny);
[X,Y] = meshgrid(x,y);
num_frame = 5;
idx_frame = round(linspace(1,ntp,num_frame));
t_frame = idx_frame./(ntp);
%% show
im = cell(ntp,1);
figure('papersize',[5,4],'paperposition',[0,0,5,4]);
set(gcf,'color','w');

switch showfunc
    case 'imshow'
        for t = 1:ntp
            imshow(squeeze(rho(t,:,:)),[]);
%             imshow(squeeze(rho(t,:,end:-1:1)),[]);
            colormap default
            pause(0.1);

            frame = getframe(gcf);
            im{t} = frame2im(frame);  
        end
        
    case 'mesh'
        
        for t = 1:ntp
            mesh(X,Y,squeeze(rho(t,:,:))');
            xlabel('x','FontSize',15,'FontWeight','bold')
            ylabel('y','FontSize',15,'FontWeight','bold')
            zlabel('\rho','FontSize',15,'FontWeight','bold')
            zlim([zmin,zmax])
            pause(0.05);

            frame = getframe(gcf);
            im{t} = frame2im(frame);    
        end
        
    case 'contour'
        for t = 1:ntp
            contour(X,Y,squeeze(rho(t,:,:))','LineWidth',1,'ShowText','on');
            xlabel('x','FontSize',15,'FontWeight','bold')
            ylabel('y','FontSize',15,'FontWeight','bold')
            pause(0.05);

            frame = getframe(gcf);
            im{t} = frame2im(frame);    
        end
end

%%
if filename   
    
    % gif
    for idx = 1:ntp
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    
%     % --- initial density
%     fig = figure('papersize',[5,5],'paperposition',[0,0,5,5]);
%       
%     switch showfunc
%         case 'imshow'
%             imshow(squeeze(rho(1,:,end:-1:1)),[]);
%             colormap default
% 
%         case 'mesh'
%             mesh(X,Y,squeeze(rho(1,:,:)));
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
%             zlabel('\rho','FontSize',15,'FontWeight','bold')
%             zlim([zmin,zmax]);
% 
%         case 'contour'
%             contour(X,Y,squeeze(rho(1,:,:))','LineWidth',1,'ShowText','on');
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
% 
%     end
%     title('t=0');
%     print(fig,'-dpdf',[filename,'_0']); 
%     
%     % --- final density
%     fig = figure('papersize',[5,5],'paperposition',[0,0,5,5]);
%        
%     switch showfunc
%         case 'imshow'
%             imshow(squeeze(rho(end,:,end:-1:1)),[]);
%             colormap default
% 
%         case 'mesh'
%             mesh(X,Y,squeeze(rho(end,:,:)));
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
%             zlabel('\rho','FontSize',15,'FontWeight','bold')
%             zlim([zmin,zmax]);
% 
%         case 'contour'
%             contour(X,Y,squeeze(rho(end,:,:))','LineWidth',1,'ShowText','on');
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
% 
%     end
%     title('t=1');
%     print(fig,'-dpdf',[filename,'_1']); 
%     
    % pdf
    fig = figure('papersize',[5*num_frame,4],'paperposition',[0,0,5*num_frame,4]);
    for k = 1:num_frame
        idx = idx_frame(k);
        subplot(1,num_frame,k);
        
        switch showfunc
            case 'imshow'
                imshow(squeeze(rho(idx,:,:)),[]);
%                 imshow(squeeze(rho(idx,:,end:-1:1)),[]);
%                 colormap default

            case 'mesh'
                mesh(X,Y,squeeze(rho(idx,:,:)));
                xlabel('x','FontSize',15,'FontWeight','bold')
                ylabel('y','FontSize',15,'FontWeight','bold')
                zlabel('\rho','FontSize',15,'FontWeight','bold')
                zlim([zmin,zmax]);
                
            case 'contour'
                contour(X,Y,squeeze(rho(idx,:,:))','LineWidth',1,'ShowText','on');
                xlabel('x','FontSize',15,'FontWeight','bold')
                ylabel('y','FontSize',15,'FontWeight','bold')
            
        end
        
        title(['t=',num2str(t_frame(k))]);
    end

%     print(fig,'-dpdf',filename); 


    
end

end