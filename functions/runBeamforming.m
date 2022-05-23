function [setup, conditions, spectra] = runBeamforming(dataFiles, setup)
% CALLS BEAMFORMING.m USING FILENAMES. FIRST FETCHES CONDITIONS FROM
% SOURCES SPECIFIED IN "SETUP" AND COMPLETES "SETUP" BASED ON CONDITIONS
%
% INPUTS
% -------
% dataFiles             string or cell  Path to h5 data file(s)
% setup.T_source        string          Source for temperature value
% setup.Re_source       string          Source for Reynolds number
% setup.mu_source       string          Source for viscosity value
%
% Note: Calls beamforming.m which has own input requirements on "setup"!
% -------

%% Get conditions and specific setup details for file
% Select file to extract conditions from & define name
if ischar(dataFiles)
    
    % Single file
    dataFile1 = dataFiles;
    [~, nameBase, ~] = fileparts(dataFile1);
    
elseif iscell(dataFiles)
    
    % Multiple files, extract name base
    dataFile1 = dataFiles{1};
    [~, nameBase, ~] = fileparts(dataFile1);
    
    if isfield(setup, 'fileGrouper')
        % Remove file grouper from config name
        if contains(nameBase, setup.fileGrouper)
            idx = strfind(nameBase, setup.fileGrouper)-1;
            nameBase = nameBase(1:idx);
        end
    end
    
else
    error('input "dataFiles" should be a string or a set of strings');
    
end
setup.name = nameBase;

% Temperature
T_default = 20;
if contains(setup.T_source, 'data', 'IgnoreCase', true)                         % From data file
    T = h5read(dataFile1,'/Conditions/Temperature');
    fprintf('Temperature from data file: %0.1f \n', T);
    
elseif contains(setup.T_source, 'override', 'IgnoreCase', true)                 % Override value
    if isfield(setup, 'T_override')
        T = setup.T_override;
        fprintf('Explicitly setting temperature: %0.1f \n', T);
        
    else
        warning(['Temperature mode set to override, but no override '...
            'value specified. Using default: ' num2str(T_default) ' C']);
        T = T_default;
    end
    
else                                                                            % Default value
    fprintf('Using default value for temperature (%i C) \n', T_default);
    T = T_default;
    
end
conditions.T = T;


% Speed of sound [m/s]
conditions.c     = 331.5 + (0.6*conditions.T);     


% Reynolds number
Re_default = 0;
if contains(setup.Re_source, 'data', 'IgnoreCase', true)                        % From data file
    Re = h5read(dataFile1,'/Conditions/Reynolds');
    fprintf('Reynolds number from data file: %0.2E \n', Re);
    
elseif contains(setup.Re_source, 'override', 'IgnoreCase', true)                % Override value
    if isfield(setup, 'Re_override')
        Re = setup.Re_override;
        fprintf('Explicitly setting Reynolds number: %0.2E \n', Re);
        
    else
        warning(['Reynolds mode set to override, but no override '...
            'value specified. Using default: ' num2str(Re_default)]);
        Re = Re_default;
    end
    
elseif contains(setup.Re_source, 'name', 'IgnoreCase', true)                    % From file name
    [~, str, ~] = fileparts(dataFile1);
    idx = strfind(str, 'Re') + 2;
    
    if isempty(idx)
        error('Trying to extract Re from file name, but "Re" not found');
    else
        seps = strfind(str, '_');
        dIdx = 0;
        if strcmp(str(idx), '_')
            dIdx = 1;
        end
        next = find(seps > idx+dIdx, 1, 'first');
        
        if isempty(next)
            idxEnd = length(str); % Value at end of name
        else
            idxEnd = seps(next)-1;
        end

        Re = 1e6 * str2double(str(idx+dIdx : idxEnd));

        fprintf('Reynolds number extracted from file name: %0.2E \n', Re);
    end
      
    if isnan(Re)
        error('Trying to extract Reynolds from file name but got NaN');
    end
    
else                                                                            % Default value
    warning(['Using fallback value for Re: ' num2str(Re_default)]);
    Re = Re_default;
    
end
    
conditions.Re = Re;

% Viscosity
mu_default = 1.825e-5;
if contains(setup.mu_source, 'data', 'IgnoreCase', true)                        % From data file
    mu = h5read(dataFile1,'/Conditions/Viscosity');
    fprintf('Viscosity from data file: %0.2E \n', mu);
    
elseif contains(setup.mu_source, 'override', 'IgnoreCase', true)                % Override value
    if isfield(setup, 'mu_override')
        mu = setup.mu_override;
        fprintf('Explicitly setting viscosity: %0.4E \n', mu);
        
    else
        warning(['Viscosity mode set to override, but no override '...
            'value specified. Using default: ' num2str(mu_default) ' kg/ms']);
        mu = mu_default;
    end
    
else                                                                            % Default value
    fprintf('Using default value for viscosity (%0.4E kg/ms) \n', mu_default);
    mu = mu_default;
    
end
conditions.mu = mu;


% Air density
dataFields = h5info(dataFile1, '/Conditions');
dataFields = dataFields.Datasets;
rho = 0;
for i=1:length(dataFields)
    if strcmp(dataFields(i).Name, 'Density')
        rho = h5read(dataFile1,'/Conditions/Density');
        fprintf('Air density from data file: %0.2f \n', rho);
    end
end
if rho == 0
    rho = 1.2;
    fprintf('Using default value for air desnity: %0.3f kg/m3', rho);
end
conditions.rho = rho;


% Settings for windtunnel environments
if isfield(setup, 'wing')
    
    % Determine AoA
    AoA_default = 0;
    if contains(setup.AoA_source, 'data', 'IgnoreCase', true)                   % From data file
        AoA = h5read(dataFile1,'/Conditions/AngleOfAttack');
        fprintf('AoA from data file: %0.1f \n', AoA);
        
    elseif contains(setup.AoA_source, 'override', 'IgnoreCase', true)           % Override value
        if isfield(setup, 'AoA_override')
            AoA = setup.AoA_override;
            fprintf('Explicitly setting AoA: %0.1f \n', AoA);

        else
            warning(['AoA mode set to override, but no override '...
                'value specified. Using default: ' num2str(T_default) ' C']);
            AoA = AoA_default;
        end
        
    elseif contains(setup.AoA_source, 'name', 'IgnoreCase', true)               % From file name
        [~, str, ~] = fileparts(dataFile1);
        idx = strfind(str, 'AoA') + 3;
        
        if isempty(idx)
            error('Trying to extract AoA from file name, but "AoA" not found');
        else
        
            dIdx = 0;   % offset between 'AoA' and the value in the name
            C = 1;      % Sign of AoA
            if strcmp(str(idx), 'n')
                C = -1;
                dIdx = 1;
            elseif strcmp(str(idx), 'p')
                dIdx = 1;
            end

            if strcmp(str(idx+dIdx), '_')
                dIdx = dIdx+1;
            end

            seps = strfind(str, '_');
            next = find(seps > idx+dIdx, 1, 'first');
            
            if isempty(next)
                idxEnd = length(str); % Value at end of name
            else
                idxEnd = seps(next)-1;
            end

            AoA = C * str2double(str(idx+dIdx : idxEnd));
            fprintf('AoA extracted from file name: %0.1f deg\n', AoA);
            
        end
        
        if isnan(AoA)
            error('Trying to extract AoA from file name but got NaN');
        end
        
    else                                                                        % Default value
        warning('Using default value for AoA: 0 deg');
        AoA = AoA_default;
        
    end
    conditions.AoA = AoA;
    
    % Correct distance based on AoA
    if contains(setup.wing.visible, 'suction', 'IgnoreCase', true)
        setup.h = setup.h0 + setup.wing.hRot*sind(conditions.AoA);
    else 
        setup.h = setup.h0 - setup.wing.hRot*sind(conditions.AoA);
    end
    
    % Determine velocity
    conditions.U   = Re / setup.wing.chord * mu / rho;
    
    % Store airfoil location
    if strcmp(setup.wing.rotAxis, 'vertical')
        setup.wing.loc = [...
            setup.wing.TEpos - setup.flowVector(1) * (setup.wing.chord + (setup.wing.chord - setup.wing.hRot)*(1-cosd(conditions.AoA))),...
            setup.scanPlaneLimits(3),...
            setup.wing.TEpos - setup.flowVector(1) * (setup.wing.hRot*(1-cosd(conditions.AoA))),...
            setup.scanPlaneLimits(4)];
            
         
    elseif strcmp(setup.wing.rotAxis, 'horizontal')
        setup.wing.loc = [...
            setup.scanPlaneLimits(1),...
            setup.wing.TEpos - setup.flowVector(2) * (setup.wing.chord + (setup.wing.chord - setup.wing.hRot)*(1-cosd(conditions.AoA))),...
            setup.scanPlaneLimits(2),...
            setup.wing.TEpos - setup.flowVector(2) * (setup.wing.hRot*(1-cosd(conditions.AoA)))];
        
    else
        error('Airfoil rotation axis should be set to horizontal or vertical')
    end
    
end


% Effective Mach number
M                = conditions.U / conditions.c;
conditions.M_eff = M * (setup.h - setup.hShear) / setup.h;

%% Run beamforming
spectra = beamforming(dataFiles, conditions, setup);
