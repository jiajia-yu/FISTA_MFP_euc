clear
clc

opts_display = [];
opts_display.num_frame = 5;
opts_display.vec_dens = 16;
opts_display.vec_leng = 1;
%%
close all
result_name = 'maze1';
load(result_name);

show_movement(rho_ot,mx_ot,my_ot,opts_display,[result_name,'_ot'])
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,[result_name,'_mfp1'])
show_movement(rho_mfp2,mx_mfp2,my_mfp2,opts_display,[result_name,'_mfp2'])
show_movement(rho_mfg,mx_mfg,my_mfg,opts_display,[result_name,'_mfg'])

%%
close all
result_name = 'maze2';
load(result_name);

show_movement(rho_ot,mx_ot,my_ot,opts_display,[result_name,'_ot'])
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,[result_name,'_mfp1'])
show_movement(rho_mfp2,mx_mfp2,my_mfp2,opts_display,[result_name,'_mfp2'])
show_movement(rho_mfg,mx_mfg,my_mfg,opts_display,[result_name,'_mfg'])

%%
close all
result_name = 'maze3';
load(result_name);

show_movement(rho_ot,mx_ot,my_ot,opts_display,[result_name,'_ot'])
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,[result_name,'_mfp1'])
show_movement(rho_mfp2,mx_mfp2,my_mfp2,opts_display,[result_name,'_mfp2'])
show_movement(rho_mfg,mx_mfg,my_mfg,opts_display,[result_name,'_mfg'])

%%
close all
result_name = 'maze4';
load(result_name);

show_movement(rho_ot,mx_ot,my_ot,opts_display,[result_name,'_ot'])
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,[result_name,'_mfp1'])
show_movement(rho_mfp2,mx_mfp2,my_mfp2,opts_display,[result_name,'_mfp2'])
show_movement(rho_mfg,mx_mfg,my_mfg,opts_display,[result_name,'_mfg'])

%%
close all
result_name = 'mfp1_10_500';
load(result_name);
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,result_name);

%%
close all
result_name = 'mfp1_10_2000';
load(result_name);
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,result_name);

%%
close all
result_name = 'mfp1_1_2000';
load(result_name);
show_movement(rho_mfp1,mx_mfp1,my_mfp1,opts_display,result_name);


