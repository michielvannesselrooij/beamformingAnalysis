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
if contains(action, 'inspect', 'IgnoreCase', true) &&...
       contains(action, 'raw', 'IgnoreCase', true)
   
    [setup, conditions, spectra] = inspectSample(setup, 'raw');
    
% Inspect a raw data sample?
elseif contains(action, 'inspect', 'IgnoreCase', true) &&...
       contains(action, 'processed', 'IgnoreCase', true)
   
    [setup, conditions, spectra] = inspectSample(setup, 'processed');
        
% Batch process raw data?
elseif contains(action, 'batch', 'IgnoreCase', true)
    
    
% Convert to database entries?
elseif contains(action, 'convert', 'IgnoreCase', true)
    % TODO to be completed
    
else
    error('Unknown action specified');
end

