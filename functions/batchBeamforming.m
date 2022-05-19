function [setup, conditions, spectra] = batchBeamforming(setup)
%% Select files

dataFiles = dir('*.h5');

if length(dataFiles) < 1
    error('No data files of correct type found. Place main script in folder with data');
end

%% Create output location
outputFolder = 'batchProcessResults';
mkdir(outputFolder);

%% Do beamforming and integration

while length(dataFiles) >= 1
    
    % Get first remaining file and remove from list
    dataFile = dataFiles(1);
    dataFiles(1) = [];
    fprintf('\nInspecting file %s \n', dataFile.name);
    filePath = [dataFile.folder filesep dataFile.name];
    
    % Check if file is part of a set
    if isfield(setup, 'fileGrouper')
        if contains(dataFile.name, setup.fileGrouper)
            
            % Look for other files in set
            idx = strfind(dataFile.name, setup.fileGrouper)-1;
            nameBase = dataFile.name(1:idx);
            k=1;
            jj = [];
            for j=1:length(dataFiles)
                if contains(dataFiles(j).name, nameBase)
                    
                    % Create set of filepaths and add original first
                    if k == 1
                        filePaths{1} = filePath;
                    end
                    k=k+1;
                    
                    % Add new file to set and put on remove list 
                    newPath = [dataFile.folder filesep dataFiles(j).name];
                    filePaths{k} =  newPath;
                    jj = [jj, j];
                    
                    fprintf('              & %s \n', dataFiles(j).name);
                    
                end
            end
            dataFiles(jj) = []; % Remove clustered files from list
        end
    end
    filePath = filePaths;
    
    % Beamforming
    [setup, conditions, spectra] = runBeamforming(filePath, setup);
    
    % Integration
    setup       = defineIntegrationWindow(setup);
    spectra.SPL = selectIntegrationWindow(setup, conditions, spectra);

    % Save result
    if length(filePath) > 1
        idx = strfind(dataFile.name, setup.fileGrouper)-1;
        name = dataFile.name(1:idx);
    else
        name = dataFile.name;
    end
    save([outputFolder filesep name '.mat'],...
        'conditions', 'setup', 'spectra');
end