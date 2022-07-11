function show_movement(rho,mx,my,Qx,opts,filename)
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
end

if isfield(opts,'num_frame') num_frame = opts.num_frame; else num_frame = 5; end
if isfield(opts,'vec_dens') vec_dens = opts.vec_dens; else vec_dens = 16; end
if isfield(opts,'vec_leng') vec_leng = opts.vec_leng; else vec_leng = 1; end

idx_frame = round(linspace(1,ntp,num_frame+1));
idx_frame = (idx_frame(1:end-1)+idx_frame(2:end))/2;
t_frame = (idx_frame-1)./(ntp-1);

X = repmat( (1:vec_dens:nx) ,floor(ny/vec_dens),1 );
Y = repmat( (1:vec_dens:ny)',1, floor(nx/vec_dens));


%% show
im = cell(ntp,1);
figure('papersize',[7,7],'paperposition',[0,0,7,7]);
set(gcf,'color','w');

for t = 1:ntp-1
    eg_illustration(squeeze(rho(t,:,:)),Qx);
%     imshow(squeeze(rho(t,:,:)),[]);
%     colormap default
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
    
    % pdf
    fig = figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
    figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    for k = 1:num_frame
        idx = idx_frame(k);
        eg_illustration(squeeze(rho(idx,:,:)),Qx);
%         imshow(squeeze(rho(idx,:,:)),[]);
%         colormap default
        hold on
%         colorbar
        
        quiver(X,Y,squeeze(my(idx,1:vec_dens:end,1:vec_dens:end)),...
                   squeeze(mx(idx,1:vec_dens:end,1:vec_dens:end)),...
                vec_leng,'r','LineWidth',1);

        
        title(['t=',num2str(t_frame(k))]);
        exportgraphics(t,[filename,'_',num2str(k),'.pdf'],'BackgroundColor','none')
    end
    
   
    

    print(fig,'-dpdf',filename); 


    
end

end