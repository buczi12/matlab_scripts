% This script presents how to generate dataset for simulation purposes
% when there are many variables in the (Simulink) model
% as an example this script will prepare dataset,
% run model simulation, collect and save output data.

% The model is simple series RLC electrical circuit and 
% we want to get step response for various R, L and C parameters

% Simulink model is not available via github
%% Workspace parameters
% folder for output data
resDir = './simResults/';
addpath(resDir)

%% Simulation parameters
% sampling frequency
fs = 1e5;                                         % [Hz]
% simulation time
tsim = 0.01;                                      % [s]

%% Model parameters
% constants
Vin = 5;                                          % [V]
% variables
Rs = [5, 9];                                      % [ohm]
Ls = [6, 9];                                      % [mH]
Cs = [1, 2, 5];                                   % [uF]
%... and so on

%% Cases Generation
% get total number of datasets
caseNum = length(Rs) * length(Ls) * length(Cs);
% preallocate structure array to save time
casesData(caseNum).Vin = Vin;
k = 1;
for R = Rs
    for L = Ls
        for C = Cs
            casesData(k).tsim = tsim;
            casesData(k).fs = fs;
            casesData(k).Vin = Vin;
            casesData(k).R = R;
            casesData(k).L = L;
            casesData(k).C = C;
            %... and so on
            % create filename for output data 
            % from k-th simulation
            casesData(k).filename = ['R', num2str(R),...
                                    'L', num2str(L),...
                                    'C', num2str(C)];
            % output namefile format will be eg. "R2L1C4.mat"
            
            k = k + 1;
        end
    end
end
save('casesData', 'casesData')
