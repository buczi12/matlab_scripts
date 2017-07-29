%% Visualization
clear all, close all %#ok<CLSCR>
resDir = './simResults/';

% List all mat files in resDir
fileNames = dir([resDir, '*.mat']);

figure(1)
 for m = 1:length(fileNames)
     % plot just data subset
     if(strcmp(fileNames(m).name(5:6), 'C5'))
         load([resDir, fileNames(m).name])
         plot(result.time, result.signals.values)
         hold on
     end
 end
 
title('RLC circuit V_o_u_t')
grid on
xlabel('time [s]'); ylabel('Voltage [V]')
