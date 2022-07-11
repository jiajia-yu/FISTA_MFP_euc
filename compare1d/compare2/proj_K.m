function [a,b] = proj_K(alpha,beta)
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
b = beta;

alpha = (alpha(1:end-1,:)+alpha(2:end,:))/2;
beta = (beta(:,1:end-1) +beta(:,2:end) )/2;
ind = (alpha + beta.^2/2 > 0);

polyb = -alpha(ind)-1;
polyd = -beta(ind).^2/2;
t = ones(size(alpha));
t(ind) = solve_cubic(1,polyb,0,polyd);
alpha = -t+1+alpha;
beta = beta./t;    

alpha = (alpha(1:end-1,:)+alpha(2:end,:))/2;
beta = (beta(:,1:end-1)+beta(:,2:end))/2;

a(2:end-1,:) = alpha;
b(:,2:end-1) = beta;


end