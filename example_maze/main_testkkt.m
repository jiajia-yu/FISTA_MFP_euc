% kkt
clear
% load('results\maze_20');
load('results\maze2');

[res1,res2] = kkt_mfp(rho_mfp,mx_mfp,my_mfp,opts);
disp(res1);
disp(res2');

%% lambda_P = 20
% fprintf('\nlambda_P = 20\n')
% 
% [res1,res2] = kkt_mfp(rho_mfp1,mx_mfp1,my_mfp1,opts_mfp1);
% disp(res1);
% disp(res2');
% 
% %% lambda_P = 200
% fprintf('\nlambda_P = 200\n')
% 
% [res1,res2] = kkt_mfp(rho_mfp2,mx_mfp2,my_mfp2,opts_mfp2);
% disp(res1);
% disp(res2');
% 
% %% lambda_P = 2000
% fprintf('\nlambda_P = 2000\n')
% 
% [res1,res2] = kkt_mfp(rho_mfp3,mx_mfp3,my_mfp3,opts_mfp3);
% disp(res1);
% disp(res2');
% 
% %%
% Qx = opts.Qx/2000;
% 
% opts.Qx = Qx*20;
% [res_stat, res_feas] = kkt_mfp(rho_mfp1,mx_mfp1,my_mfp1,opts);
% disp(res_stat);
% disp(res_feas');
% 
% opts.Qx = Qx*200;
% [res_stat, res_feas] = kkt_mfp(rho_mfp2,mx_mfp2,my_mfp2,opts);
% disp(res_stat);
% disp(res_feas');
% 
% opts.Qx = Qx*2000;
% [res_stat, res_feas] = kkt_mfp(rho_mfp3,mx_mfp3,my_mfp3,opts);
% disp(res_stat);
% disp(res_feas');

