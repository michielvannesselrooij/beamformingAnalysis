function [setup, conditions, spectra] = batchBeamforming(setup)
% DO BEAMFORMING ON ALL DATA IN CURRENT FOLDER AND STORE AS .mat,
% OPTIONALLY COMBINING FILES BASED ON SPECIFIED GROUPER STRING IN NAME
%
% INPUTS
% -------
% setup            structure   Passed to runBeamforming.m, see requirements
% -------


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
                
                % Create set of filepaths and add original first
                if k == 1
                    filePaths{1} = filePath;
                end
                    
                if contains(dataFiles(j).name, nameBase)
                    
                    k=k+1;
                    
                    % Add new file to set and put on remove list 
                    newPath = [dataFile.folder filesep dataFiles(j).name];
                    filePaths{k} =  newPath;
                    jj = [jj, j];
                    
                    fprintf('              & %s \n', dataFiles(j).name);
                    
                end
            end
            dataFiles(jj) = []; % Remove clustered files from list
            filePath = filePaths;
        end
    end
    
    % Determine output file
    if iscell(filePath)
        idx = strfind(dataFile.name, setup.fileGrouper)-1;
        name = dataFile.name(1:idx);
    else
        [~, name, ~] = fileparts(dataFile.name);
    end

    outputFileName = [outputFolder filesep name '.mat'];
    if exist(outputFileName, 'file') == 2
        fprintf('Skipping, .mat file already exists');
        
    else
        % Beamforming
        [setup, conditions, spectra] = runBeamforming(filePath, setup);

        % Integration
        setup       = defineIntegrationWindow(setup);
        spectra.SPL = selectIntegrationWindow(setup, conditions, spectra);

        % Save result
        save(outputFileName, 'conditions', 'setup', 'spectra');
        
        % Create and store figures
        nPlots = 10;
        beamformingSummary(setup, spectra, nPlots)
        
        figureFolder = [outputFolder filesep  'figures_' name];
        mkdir(figureFolder);
        
        figList = handle( sort( double(findall(0, 'type', 'figure') ) ));
        savefig(figList, [figureFolder filesep 'figures'], 'compact');
        for i=1:length(figList)
            set(figList(i),'WindowStyle','normal');
            set(figList(i), 'Units', 'pixels', 'Position', [10 10 800 600]);
            saveas(figList(i), [figureFolder filesep 'figure_' num2str(i) '.png']);
        end
        for i=1:length(figList)
            close;
        end
        
    end

    fprintf('\nDone. Files in queue: %i \n' ,length(dataFiles));
end