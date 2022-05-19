function [setup, conditions, spectra] = processBeamformingConfig(folder)
% MAIN FUNCTION DETERMINING WHICH FUNCTIONS TO RUN, BASED ON ACTION
% AND CONFIGURATION RETRIEVED FROM CONFIG.xlsx

% INPUTS
% -------
% folder        string      path to folder with the config file to consider
% -------

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
    
    [setup, conditions, spectra] = batchBeamforming(setup);
    
% Convert to database entries?
elseif contains(action, 'convert', 'IgnoreCase', true)
    
    offset = input('Database id after which to add results: ');
    
    dataFolder = 'batchProcessResults';
    if exist(dataFolder, 'dir') ~= 7
        dataFolder = pwd;
        fprintf(['Cannot find folder %s, looking for .mat files in ',...
            ' current folder \n'], dataFolder);
    end
    
    convertToDatabaseBF(dataFolder, offset)
    
else
    error('Unknown action specified');
end

