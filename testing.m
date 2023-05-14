close all;
clear;
clc;

vf  = 1;
v_guess  = 60;
numBlocks   = 100;
vel_block   = ones(1,numBlocks+1)*5;
vel_air     = ones(1,numBlocks+1)*5;
crk_len     = 0.5;      % in

for bloq   = 1:numBlocks
    if (exist("crk_len", 'var'))
        vel_block(bloq+1)  = func_MDsecant(@(v0)func_vel_resid(v0, vel_air(bloq), crk_len), vel_air(bloq), 1e-3, 100000, 1);
    else
        vel_block(bloq+1)  = func_MDsecant(@(v0)func_vel_resid(v0, vel_air(bloq)), vel_block(bloq)*2, 1e-3, 100000, 1);
    end
    vel_air(bloq+1)     = func_MDsecant(@(v0)func_vel_air_resid(v0, vel_block(bloq+1)), vel_block(bloq+1)+1, 1e-3, 100000,1);
    disp(" ");
    disp("---------------------------------------------");
    disp("I have cleared Block " + bloq + "!");
    disp("---------------------------------------------");
    disp(" ");
end

vel_block(1)     = 0;

figure(1);
clf;
subplot(2,1,1);
plot(0:numBlocks, vel_air, "o");
xlabel("Block #");
ylabel("Speed Needed to Clear (mph)");

vel_lost_block  = vel_block(2:end) - vel_air(1:end-1);
vel_lost_air    = vel_air(2:end) - vel_block(2:end);
subplot(2,1,2);
plot(numBlocks:-1:1, vel_lost_block, "o");
hold on;
plot(numBlocks:-1:1, vel_lost_air, "o")
xlabel("Block #");
ylabel("Speed Lost Across Block (mph)");
legend("From Block", "From Air");

% Percentage from each one
vel_lost    = vel_lost_air + vel_lost_block;
lost_percent_bloq   = (vel_lost_block./vel_lost)*100;
lost_percent_air    = (vel_lost_air./vel_lost)*100;

% Plot lost velocities
figure(2);
clf;
plot(numBlocks:-1:1, lost_percent_bloq, "o");
hold on;
plot(numBlocks:-1:1, lost_percent_air, "o");
xlabel("Block #");
ylabel("Percentage of Speed Lost");
legend("From Block", "From Air");

disp("Velocity needed to pass through seven blocks: " + vel_air(8) + " mph");
