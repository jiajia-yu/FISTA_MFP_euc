clear
load('results\image_20');

%% MFG
fprintf('\nmean field game\n')

[res1,res2] = kkt_mfg(rho_mfg,mx_mfg,my_mfg,opts_mfg);
disp(res1);
disp(res2');

%% OT
fprintf('\noptimal transport\n')
[res1,res2] = kkt_ot(rho_ot,mx_ot,my_ot);
disp(res1);
disp(res2');

%% MFP1 interaction \rho*log(\rho)
fprintf('\nmean field planning 1\n')

[res1,res2] = kkt_mfp1(rho_mfp1,mx_mfp1,my_mfp1,opts_mfp1);
disp(res1);
disp(res2');

%% MFP2 interaction \rho*\rho/2
fprintf('\nmean field planning 2\n')

[res1,res2] = kkt_mfp2(rho_mfp2,mx_mfp2,my_mfp2,opts_mfp2);
disp(res1);
disp(res2');

%% MFP3 interaction 1/\rho
fprintf('\nmean field planning3\n')

[res1,res2] = kkt_mfp3(rho_mfp3,mx_mfp3,my_mfp3,opts_mfp3);
disp(res1);
disp(res2');

