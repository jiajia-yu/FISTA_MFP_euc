function [rho,mx,my,outs] = ot2d_fista_mg(rho0,rho1,level,opts)

if length(size(rho0))==2
    [nxp,nyp] = size(rho0);
    rho0 = reshape(rho0,1,nxp,nyp);
    rho1 = reshape(rho1,1,nxp,nyp);    
end
nx = nxp-1;
ny = nyp-1;
if ~isfield(opts,'nt') opts.nt = 4*2^level; end
if ~isfield(opts,'tol') opts.tol = 1e-6; end
if isfield(opts,'maxit') maxit = opts.maxit; else maxit = 1e2; end
if isfield(opts,'maxit_des') maxit_des = opts.maxit_des; else maxit_des = 10;  end
if isfield(opts,'maxit_aes') maxit_aes = opts.maxit_aes; else maxit_aes = 5e2; end

%%
nts = zeros(level,1);
nxs = zeros(level,1);
nys = zeros(level,1);
nts(level) = opts.nt;
nxs(level) = nx;
nys(level) = ny;

rho0s = cell(level,1);
rho1s = cell(level,1);
rho0s{level} = rho0;
rho1s{level} = rho1;
for l = level-1:-1:1
    nts(l) = nts(l+1)/2;
    nxs(l) = nxs(l+1)/2;
    nys(l) = nys(l+1)/2;
    
    rho0_cur = rho0s{l+1};
    rho0_cur = rho0_cur(:,1:2:end,1:2:end);
    rho0s{l} = rho0_cur;
    
    rho1_cur = rho1s{l+1};
    rho1_cur = rho1_cur(:,1:2:end,1:2:end);
    rho1s{l} = rho1_cur;
end

%%
t_total = 0;
nit_total = 0;

opts.rho = ones(nts(level)-1,nxs(level),nys(level));
opts.mx = ones(nts(level),nxs(level)-1,nys(level));
opts.my = ones(nts(level),nxs(level),nys(level)-1);

rhos = cell(level,1);
mxs = cell(level,1);
mys = cell(level,1);
drhos = cell(level-1,1);
dmxs = cell(level-1,1);
dmys = cell(level-1,1);

for nit = 1:opts.maxit
    opts.maxit = maxit_des;
    for l = level:-1:1
        tic
        opts.nt = nts(l);
        if l == level
            [rho_cur,mx_cur,my_cur,outs] = ot2d_fista(rho0s{l},rho1s{l},opts);
            rhos{l} = rho_cur;
            mxs{l} = mx_cur;
            mys{l} = my_cur;
        else
            opts.rho = restr2d_rho(rho_cur);
            opts.mx = restr2d_mx(mx_cur);
            opts.my = restr2d_my(my_cur);
            [rho_cur,mx_cur,my_cur,outs] = ot2d_fista(rho0s{l},rho1s{l},opts);
            rhos{l} = rho_cur;
            mxs{l} = mx_cur;
            mys{l} = my_cur;
            drhos{l} = rho_cur - opts.rho;
            dmxs{l} = mx_cur - opts.mx;
            dmys{l} = my_cur - opts.my;
        end
        t_l = toc;
        t_total = t_total + t_l;
        nit_total = nit_total + length(outs.objs);
    end   
    
    opts.maxit = maxit_aes;
    for l = 1:1:level
        tic
        opts.nt = nts(l);
        if l == 1
            opts.rho = rhos{1};
            opts.mx = mxs{1};
            opts.my = mys{1};
            [rho_cur,mx_cur,my_cur,outs] = ot2d_fista(rho0s{l},rho1s{l},opts);
            rhos{l} = rho_cur;
            mxs{l} = mx_cur;
            mys{l} = my_cur;
        else
            drhos{l-1} = rho_cur - opts.rho;
            dmxs{l-1} = mx_cur - opts.mx;
            dmys{l-1} = my_cur - opts.my;
            opts.rho = rhos{l} + ...
                inter2d_rho(zeros(nxs(l-1)+1,nys(l-1)+1),zeros(nxs(l-1)+1,nys(l-1)+1),drhos{l-1});
            opts.mx = mxs{l} + inter2d_mx(dmxs{l-1});
            opts.my = mys{l} + inter2d_my(dmys{l-1});
            [rho_cur,mx_cur,my_cur,outs] = ot2d_fista(rho0s{l},rho1s{l},opts);
            rhos{l} = rho_cur;
            mxs{l} = mx_cur;
            mys{l} = my_cur;
        end
        t_l = toc;
        t_total = t_total + t_l;
        nit_total = nit_total + length(outs.objs);
        fprintf('%d iters in %.3f sec\n',length(outs.objs),t_l);
    end
    
end

%%
rho = rhos{level};
mx = mxs{level};
my = mys{level};
outs.t_total = t_total;
outs.nit_total = nit_total;


