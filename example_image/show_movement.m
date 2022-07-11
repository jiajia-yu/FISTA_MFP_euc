function show_movement(rho,mx,my,opts,filename)
mx = (mx(:,1:end-1,:) + mx(:,2:end,:) )/2;
my = (my(:,:,1:end-1) + my(:,:,2:end) )/2;

[ntp,nx,ny] = size(rho);
% zmin = min(rho,[],'all');
% zmax = max(rho,[],'all');

if nargin < 4
    opts = [];
end

if nargin < 5 
    filename = 0;
% else
%     filename = [filename,'_',showfunc];
end

if isfield(opts,'num_frame') num_frame = opts.num_frame; else num_frame = 5; end
if isfield(opts,'vec_dens') vec_dens = opts.vec_dens; else vec_dens = 16; end
if isfield(opts,'vec_leng') vec_leng = opts.vec_leng; else vec_leng = 1; end

idx_frame = round(linspace(1,ntp,num_frame));
t_frame = (idx_frame-1)./(ntp-1);

X = repmat( (1:vec_dens:nx) ,floor(ny/vec_dens),1 );
Y = repmat( (1:vec_dens:ny)',1, floor(nx/vec_dens));


%% show
im = cell(ntp,1);
figure('papersize',[7,7],'paperposition',[0,0,7,7]);
set(gcf,'color','w');

% switch showfunc
%     case 'imshow'
        for t = 1:ntp-1
            imshow(squeeze(rho(t,:,:)),[]);
%             imshow(squeeze(rho(t,:,end:-1:1)),[]);
            colormap default
            hold on
            
            quiver(X,Y,squeeze(my(t,1:vec_dens:end,1:vec_dens:end)),squeeze(mx(t,1:vec_dens:end,1:vec_dens:end)),...
                vec_leng,'r','LineWidth',0.75);

            pause(0.2);

            frame = getframe(gcf);
            im{t} = frame2im(frame);  
        end
        imshow(squeeze(rho(ntp,:,:)),[]);
        colormap default
        pause(0.2);
        frame = getframe(gcf);
        im{ntp} = frame2im(frame);
        
%     case 'mesh'
%         
%         for t = 1:ntp
%             mesh(X,Y,squeeze(rho(t,:,:))');
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
%             zlabel('\rho','FontSize',15,'FontWeight','bold')
%             zlim([zmin,zmax])
%             pause(0.05);
% 
%             frame = getframe(gcf);
%             im{t} = frame2im(frame);    
%         end
%         
%     case 'contour'
%         for t = 1:ntp
%             contour(X,Y,squeeze(rho(t,:,:))','LineWidth',1,'ShowText','on');
%             xlabel('x','FontSize',15,'FontWeight','bold')
%             ylabel('y','FontSize',15,'FontWeight','bold')
%             pause(0.05);
% 
%             frame = getframe(gcf);
%             im{t} = frame2im(frame);    
%         end
% end

%%
if filename   
    
    % gif
    for idx = 1:ntp
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0.2);
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
    for k = 1:num_frame-1
        idx = idx_frame(k);
        subplot(1,num_frame,k);
        
%         switch showfunc
%             case 'imshow'
                imshow(squeeze(rho(idx,:,:)),[]);
%                 imshow(squeeze(rho(idx,:,end:-1:1)),[]);
                colormap default
                hold on
                
                quiver(X,Y,squeeze(my(idx,1:vec_dens:end,1:vec_dens:end)),squeeze(mx(idx,1:vec_dens:end,1:vec_dens:end)),...
                vec_leng,'r','LineWidth',1);


%             case 'mesh'
%                 mesh(X,Y,squeeze(rho(idx,:,:)));
%                 xlabel('x','FontSize',15,'FontWeight','bold')
%                 ylabel('y','FontSize',15,'FontWeight','bold')
%                 zlabel('\rho','FontSize',15,'FontWeight','bold')
%                 zlim([zmin,zmax]);
%                 
%             case 'contour'
%                 contour(X,Y,squeeze(rho(idx,:,:))','LineWidth',1,'ShowText','on');
%                 xlabel('x','FontSize',15,'FontWeight','bold')
%                 ylabel('y','FontSize',15,'FontWeight','bold')
%             
%         end
        
        title(['t=',num2str(t_frame(k))]);
    end
    
    idx = idx_frame(num_frame);
    subplot(1,num_frame,num_frame);

    imshow(squeeze(rho(idx,:,:)),[]);
    colormap default
    title(['t=',num2str(t_frame(num_frame))]);

%     print(fig,'-dpdf',filename); 


    
end

end