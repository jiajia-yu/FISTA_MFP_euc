function [a,bx,by] = proj_K(alpha,betax,betay)
%% description
% alpha Nt*Nx(*Ny)
% beta Nt*Nx(*Ny*2)
% entrywise, 
% (a,b) = argmin{(a-alpha)^2+||b-beta||^2} 
% s.t. a+||b||^2/2<=0

% I. if alpha+||beta||^2/2<=0, then a=alpha,b=beta
%II. if alpha+||beta||^2/2> 0, solve ||beta||^2*t^3+2(alpha+1)t-2=0

%%
a = alpha;
bx = betax;
by = betay;

alpha = (alpha(1:end-1,:,:) + alpha(2:end,:,:))/2;
betax= (betax(:,1:end-1,:) + betax(:,2:end,:))/2;
betay= (betay(:,:,1:end-1) + betay(:,:,2:end))/2;
ind = (alpha + (betax.^2+betay.^2)/2 > 0);

polyb = -alpha(ind)-1;
polyd = -(betax(ind).^2+betay(ind).^2)/2;
t = ones(size(alpha));
t(ind) = solve_cubic(1,polyb,0,polyd);
alpha = -t+1+alpha;
betax = betax./t;
betay = betay./t;

alpha  = (alpha(1:end-1,:,:) +alpha(2:end,:,:)) /2;
betax = (betax(:,1:end-1,:)+betax(:,2:end,:))/2;
betay = (betay(:,:,1:end-1)+betay(:,:,2:end))/2;

a(2:end-1,:,:) = alpha;
bx(:,2:end-1,:)= betax;
by(:,:,2:end-1)= betay;
        
end