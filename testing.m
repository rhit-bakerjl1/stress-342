close all;
clear;
clc;

vf  = 1;
v_guess  = 60;
[vel_optim]     = func_MDsecant(@(v0)func_vel_resid(v0, vf), v_guess, 1e-3, 100000, 1);
