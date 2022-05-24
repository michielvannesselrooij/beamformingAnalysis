function processBeamformingConfig(folder)
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
   
    inspectSample(setup, 'raw');
    
% Inspect a raw data sample?
elseif contains(action, 'inspect', 'IgnoreCase', true) &&...
       contains(action, 'processed', 'IgnoreCase', true)
   
    inspectSample(setup, 'processed');
        
% Batch process raw data?
elseif contains(action, 'batch', 'IgnoreCase', true)
    
    batchBeamforming(setup);
    
elseif contains(action, 'adjust', 'IgnoreCase', true)
    adjustIntegrationWindow(setup)
    
% Convert to database entries?
elseif contains(action, 'convert', 'IgnoreCase', true)
    
    dataFolder = 'batchProcessResults';
        
    if exist(dataFolder, 'dir') ~= 7
        fprintf(['Cannot find folder %s, looking for .mat files in ',...
            ' reprocessedResults \n'], dataFolder);
        dataFolder = 'reprocessedResults';
        
    elseif exist(dataFolder, 'dir') ~= 7
        fprintf(['Cannot find folder %s, looking for .mat files in ',...
            ' current folder \n'], dataFolder);
        dataFolder = pwd;
    end
    
    convertToDatabaseBF(dataFolder)
    
else
    error('Unknown action specified');
end

