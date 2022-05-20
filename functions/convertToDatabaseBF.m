function convertToDatabaseBF(dataFolder)
% READ AND COMBINE PROCESSED BEAMFORMING DATA (.mat) AND STORE SPECTRA IN
% MUTESKIN DATABASE FORMAT
%
% INPUTS
% -------
% dataFolder        string      relative or absolute path to .mat files
% -------

%% Create output location
outputFolder = 'databaseNoiseFiles';
mkdir(outputFolder);

%% Output
dataFiles = dir([dataFolder filesep '*.mat']);

% Find file sets
name = {};
for i=1:length(dataFiles)
    
    % Check file name structure
    idx1  = strfind(dataFiles(i).name, '_Re');
    idx2  = strfind(dataFiles(i).name, '_AoA');
    
    load([dataFolder filesep dataFiles(i).name],'setup');
    if isfield(setup, 'fileGrouper')
        idx3 = strfind(dataFiles(i).name, setup.fileGrouper);
    end
    
    idx = min([idx1 idx2 idx3]);
    
    % Crop to name base
    nameBase = dataFiles(i).name(1:idx-1);
    
    if sum(contains(name, nameBase)) == 0
        name = [name; nameBase];
    end
    
end
N = length(name);

% Open files and store in database format
for i=1:N
    
    disp(['Processing ' name{i}]);
    
    % Start file
    outputFile = [outputFolder filesep name{i} '.csv'];
    fid = fopen(outputFile, 'w');
    fprintf(fid, '%s\n', 'rho,T,V,alpha,alpha_eff,f,dB');
    
    % Find all data sets for this configuration
    setFiles = dir([dataFolder filesep name{i} '*.mat']);
    
    for j=1:length(setFiles)
        
        fprintf(['  %s\n' setFiles(j).name]);
        
        % Extract data
        load([setFiles(j).folder filesep setFiles(j).name],...
            'spectra', 'conditions');
        
        % Write file
        f           = spectra.f;
        SPL         = spectra.SPL;
        V           = conditions.U;
        alpha       = conditions.AoA;
        
        if isfield(conditions, 'AoA_eff')
            alpha_eff = conditions.AoA_eff;
        else
            alpha_eff = alpha;
            warning(['No effective AoA found. Assuming effective AoA'...
                ' is geometric AoA']);
        end
        
        fprintf(fid, '%s%f%s%f%s%f%s%f%s%f\n', ',,', V, ',', alpha, ',',...
            alpha_eff, ',', f(1), ',', SPL(1));

        for ii=2:length(f)
            if ~isnan(f(ii)) && ~isnan(SPL(ii))
                fprintf(fid, '%s%f%s%f\n', ',,,,,', f(ii), ',', SPL(ii));
            end
        end
            
    end 
    
    % Close file
    fclose(fid);

end

fprintf('Conversion to database ready and stored in %s \n\n', outputFolder);