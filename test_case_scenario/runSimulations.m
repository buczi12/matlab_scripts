%% Run simulations
clear all %#ok<CLSCR>
resDir = './simResults/';
load casesData.mat
% model variables' names should be formatted as simCase.R'

for n = 1:length(casesData)
    simCase = casesData(n);
    % Check if result file already exists
    if(~exist([resDir, simCase.filename, '.mat'], 'file'))
        sim RLCmodel

        % model produces n-th output data as 'result' structure
        % save data as filename in resDir flder
        save([resDir, simCase.filename], 'result')

        disp([num2str(n), ' simulations out of ',...
             num2str(length(casesData)), ' done.'])
    else
        disp([num2str(simCase.filename),...
        ' simulation result already exists.'])
    end
end

