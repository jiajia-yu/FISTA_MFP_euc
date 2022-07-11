% function [rho,mx,outs] = ot1d_fista_ml(rho0,rho1,level,opts)
% fprintf('no interpolation in time space\n')
% 
% nxp = length(rho0);
% nx = nxp-1;
% if ~isfield(opts,'nt') opts.nt = 20; end
% 
% %%
% nxs = zeros(level,1);
% nxs(level) = nx;
% 
% rho0s = cell(level,1);
% rho1s = cell(level,1);
% rho0s{level} = rho0;
% rho1s{level} = rho1;
% for l = level-1:-1:1
%     nxs(l) = nxs(l+1)/2;
%     
%     rho0_cur = rho0s{l+1};
%     rho0_cur = rho0_cur(1:2:end);
%     rho0s{l} = rho0_cur;
%     
%     rho1_cur = rho1s{l+1};
%     rho1_cur = rho1_cur(1:2:end);
%     rho1s{l} = rho1_cur;
% end
% opts.rho = ones(opts.nt-1,nxs(1));
% opts.mx = ones(opts.nt,nxs(1)-1);
% 
% %%
% t_total = 0;
% nit_total = 0;
% for l = 1:level
%     if l > 1
%         opts.rho = inter1d_rho(rho0s{l-1},rho1s{l-1},rho_coarse);
%         opts.mx = inter1d_mx(mx_coarse);
%     end
%     tic
%     [rho_coarse,mx_coarse,outs] = ot1d_fista(rho0s{l},rho1s{l},opts);
%     t_l = toc;
%     t_total = t_total + t_l;
%     nit_total = nit_total + length(outs.objs);
%     fprintf('%d iters in %.3f sec\n',length(outs.objs),t_l);
% end
% 
% %%
% rho = rho_coarse;
% mx = mx_coarse;
% outs.t_total = t_total;
% outs.nit_total = nit_total;

%% ---------------------------------------------------------------------
function [rho,mx,outs] = ot1d_fista_ml(rho0,rho1,level,opts)
fprintf('interpolation in time space\n')

nxp = length(rho0);
nx = nxp-1;
if ~isfield(opts,'nt') opts.nt = 4*2^level; end

%%
nts = zeros(level,1);
nxs = zeros(level,1);
nts(level) = opts.nt;
nxs(level) = nx;

rho0s = cell(level,1);
rho1s = cell(level,1);
rho0s{level} = rho0;
rho1s{level} = rho1;
for l = level-1:-1:1
    nts(l) = nts(l+1)/2;
    nxs(l) = nxs(l+1)/2;
    
    rho0_cur = rho0s{l+1};
    rho0_cur = rho0_cur(1:2:end);
    rho0s{l} = rho0_cur;
    
    rho1_cur = rho1s{l+1};
    rho1_cur = rho1_cur(1:2:end);
    rho1s{l} = rho1_cur;
end
opts.rho = ones(nts(1)-1,nxs(1));
opts.mx = ones(nts(1),nxs(1)-1);

%%
t_total = 0;
nit_total = 0;
for l = 1:level
    opts.nt = nts(l);
    if l > 1
        opts.rho = inter1d_rho(rho0s{l-1},rho1s{l-1},rho_coarse);
        opts.mx = inter1d_mx(mx_coarse);
    end
    tic
    [rho_coarse,mx_coarse,outs] = ot1d_fista(rho0s{l},rho1s{l},opts);
    t_l = toc;
    t_total = t_total + t_l;
    nit_total = nit_total + length(outs.objs);
    fprintf('%d iters in %.3f sec\n',length(outs.objs),t_l);
end

%%
rho = rho_coarse;
mx = mx_coarse;
outs.t_total = t_total;
outs.nit_total = nit_total;


