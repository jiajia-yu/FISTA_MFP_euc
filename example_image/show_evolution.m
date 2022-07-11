function show_evolution(rho,filename)
if nargin < 2 
    filename = 0;
end
[ntp,~,~] = size(rho);
num_frame = 6;
idx_frame = round(linspace(1,ntp,num_frame));
t_frame = (idx_frame-1)./(ntp-1);
%% show
im = cell(ntp,1);
figure('papersize',[5,5],'paperposition',[0,0,5,5]);
set(gcf,'color','w');

for t = 1:ntp
    imshow(squeeze(rho(t,:,:)),[]);
%     colormap default
    colorbar
    pause(0.1);

    frame = getframe(gcf);
    im{t} = frame2im(frame);  
end
        


%%
if filename   
    % pdf
%     fig = figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
%     figure('papersize',[5*num_frame,5],'paperposition',[0,0,5*num_frame,5]);
    t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
    nexttile
    for k = 2:num_frame
        idx = idx_frame(k);
        imshow(squeeze(rho(idx,:,:)),[]);
%         colormap default
        colorbar
        title(['t=',num2str(t_frame(k))]);
        exportgraphics(t,[filename,'_',num2str(k),'.pdf'],'BackgroundColor','none')
    end
    
    % gif
    for idx = 1:ntp
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end

end