function [setup, conditions, spectra] = processBeamformingConfig(folder)

%% Read configuration
configFile = [folder filesep 'config.xlsx'];

if ~exist('config.xlsx', 'file')
    error('Cannot run. Missing file "config.xlsx" in this folder.');
else
    [action, setup] = readBeamformingConfig(configFile);
end

%% Perform specified action

% Inspect a raw data sample?
if contains(action, 'inspect', 'IgnoreCase', true)
    [setup, conditions, spectra] = inspectRawSample(setup);
        
% Batch process raw data?
elseif contains(action, 'process', 'IgnoreCase', true)
    
    
% Convert to database entries?
elseif contains(action, 'convert', 'IgnoreCase', true)
    % TODO to be completed
    
else
    error('Unknown action specified');
end

