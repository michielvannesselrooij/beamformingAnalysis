function [action, setup] = readBeamformingConfig(file)
% EXTRACTS ANALYSIS CONFIGURATION SETTINGS TO "SETUP" FROM THE CONFIG.xlsx
% FILE AND DETERMINES THE ACTIONS TO PERFORM
%
% INPUTS
% -------
% file        string      relative or absolute path to config.xlsx
% -------

%% Load beamforming configuration file

data                        = importdata(file);
dataIdx                     = ~isnan(data.data.config);
config.data                 = nan(size(data.data.config));
config.data(dataIdx)        = data.data.config(dataIdx);
config.textdata             = data.textdata.config;

clear data

%% Extract data that is always relevant

action                      = config.textdata{1,2};    % Analysis action

setup.h0                    = config.data(2,   1);     % Target distance (nominal)
setup.fRange                = config.data(1,   4:5);   % Beamforming frequency range
setup.fs                    = config.data(2,   4);     % Acquisition frequency
setup.dataPortion           = config.data(3,   4)/100; % Portion of data to use (centered around middle)
setup.timeChunk             = config.data(4,   4);     % Data chunk size in seconds
setup.scanPlaneResolution   = config.data(5,   4);     % Beamforming resolution
setup.diagonalRemoval       = config.data(6,   4);     % Remove CSM diagonal?
setup.brokenMics            = config.data(7:12,4);     % Indices of broken mics

setup.fileGrouper           = config.textdata{12, 10}; % String to look for in filenames to group them

%% Specify source of data variables

sourceSettings = config.textdata([4, 6, 8, 10], 10);
overrides      = config.data(    [2, 4, 6, 8], 9);

for i=1:length(sourceSettings)
    if contains(sourceSettings{i}, 'Override')
        if isnan(overrides(i))
            error('Override selected, but no override value specified')
        end
    end
end

setup.T_source     = sourceSettings{1};
setup.Re_source    = sourceSettings{2};
setup.mu_source    = sourceSettings{3};
setup.AoA_source   = sourceSettings{4};

setup.T_override   = overrides(1);
setup.Re_override  = overrides(2);
setup.mu_override  = overrides(3);
setup.AoA_override = overrides(4);

%% Checks
if setup.dataPortion > 1
    error('Data portion to be used should be 100% at most');
end

if setup.diagonalRemoval ~= 0 && setup.diagonalRemoval ~= 1 && ~isnan(setup.diagonalRemoval)
    error('CSM diagonal removal should be set to 0 or 1');
end

if isempty(setup.fileGrouper) && ~isempty(config.data(9,9))
    error(['String in name to group files set to ' num2str(config.data(9,9)) ...
        ', but should be a string (probably starting with "_")']);
end

% Clean up broken mic selection
setup.brokenMics(isnan(setup.brokenMics)) = [];
setup.brokenMics = sort(unique(setup.brokenMics));

% Delete empty fields
Sfields                     = fieldnames(setup);
tc                          = struct2cell(setup);
cn0                         = cellfun(@isnan, tc, 'UniformOutput', false);
cn                          = cellfun(@sum, cn0);
mask                        = any(cn, 2:ndims(cn));
tc(any(cn, 2:ndims(cn)),:)  = [];
setup                       = reshape(cell2struct(tc, Sfields(~mask), 1), size(setup));

%% Load additional data based on measurement type
facility = config.textdata{4,2};
setup.facility = facility;

if contains(facility, 'windtunnel', 'IgnoreCase', true) ||...
   contains(facility, 'wind tunnel', 'IgnoreCase', true)
    
    setup.hShear        = config.data(8,1);
    
    setup.wing.chord    = config.data(5,1);
    setup.wing.TEpos    = config.data(6,1);
    setup.wing.hRot     = config.data(7,1);
    setup.wing.visible  = config.textdata{12,2};
    setup.intPlane_rel  = config.data(11:14,1);
    
    if setup.wing.hRot < 0
        setup.wing.hRot = abs(setup.wing.hRot);
        warning(['Negative distance TE to center of rotation specified,'...
            ' assuming absolute value']);
    end
    
    if setup.hShear < 0
        setup.hShear = abs(setup.hShear);
        warning(['Negative distance mics to shear layer specified,'...
            ' assuming absolute value']);
    end
   
else
    % If not in a windtunnel, no change to target distance expected
    setup.h = setup.h0;
    
end

%% Load facility-specific configurations
setup = loadFacilityDetails(setup);

