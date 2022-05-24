function adjustIntegrationWindow(setup)
% LOAD PROCESSED BEAMFORMING DATA AND ADJUST THE INTEGRATION WINDOW
%
% INPUTS
% -------
% setup            structure   Passed to  defineIntegrationWindow and selectIntegrationWindow
% -------


%% Select files

% Look in current folder
dataFiles = dir('*.mat');

if isempty(dataFiles)
    
    % Look in default batch processed folder
    dataFiles = dir(['batchProcessResults' filesep '*.mat']);
    
    if isempty(dataFiles)
        error(['No .mat files found. Place main script in folder with data' ...
            ' or have a folder called "batchProcessResults"']);
    end
end

fprintf('Reprocessing data with new integration window settings...\n');

%% Create output location
outputFolder = 'reprocessedResults';
mkdir(outputFolder);

%% Do beamforming and integration

while length(dataFiles) >= 1
    
    % Get first remaining file and remove from list
    dataFile = dataFiles(1);
    dataFiles(1) = [];
    fprintf('\nInspecting file %s \n', dataFile.name);
    filePath = [dataFile.folder filesep dataFile.name];
    
    % Determine output file
    [~, name, ~] = fileparts(dataFile.name);

    outputFileName = [outputFolder filesep name '.mat'];
    if exist(outputFileName, 'file') == 2
        fprintf('Skipping, .mat file already exists');
        
    else
        
        % Load data and substitute new integration window
        intPlane_rel_new   = setup.intPlane_rel;
        load(filePath, 'setup', 'conditions', 'spectra');
        setup.intPlane_rel = intPlane_rel_new;

        % Integration on window
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