function [setup, conditions, spectra] = inspectSample(setup, inputType)
% PROCESS (PART OF) FIRST FILE IN FOLDER AND SHOW RESULTS SUCH AS
% MICROPHONE LAYOUT, BEAMFORMING + GEOMETRY + INTEGRATION WINDOW, AND
% FINAL SPECTRUM. OPTIONALLY STORE RESULTS
%
% INPUTS
% -------
% inputType        string      File to load: 'raw' (.h5) or 'processed' (.mat)
% setup            structure   Passed to runBeamforming.m, see requirements
% -------

%% Select file

if strcmp(inputType, 'raw')
    dataFiles = dir('*.h5');
elseif strcmp(inputType, 'processed')
    
    % Look in inspectSample
    if exist('inspectSample', 'dir') == 7         
        dataFiles  = dir(['inspectSample' filesep '*sample*.mat']);
        
    else
        % Look in inspect_*
        if exist('inspect_*', 'dir') == 7         
            folders    = dir('inspect_*');
            dataFiles  = dir([folders(1).name filesep '*sample*.mat']);
            
        % Look in current folder
        elseif exist('*sample*.mat', 'file') == 2 
            dataFiles  = dir('*sample*.mat');
            
        else
            error('Cannot find the processed sample file');
        end
    end
end

if length(dataFiles) < 1
    error('No data files of correct type found. Place main script in folder with data');
else
    dataFile = dataFiles(1);
    fprintf('Inspecting file %s \n', dataFile.name);
    filePath = [dataFile.folder filesep dataFile.name];
end


%% Settings

if strcmp(inputType, 'raw')
    setup.dataPortion = input('Portion of data to use in this check (%): ')/100;
end

nPlots = input('Number of frequencies to plot:            ');


%% Pass to beamforming
if strcmp(inputType, 'raw')
    [setup, conditions, spectra] = runBeamforming(filePath, setup);
else
    load(filePath, 'setup', 'conditions', 'spectra');
    fprintf('Loading .mat file. Skipping beamforming.\n');
end


%% Inspect integration window
setup       = defineIntegrationWindow(setup);
spectra.SPL = selectIntegrationWindow(setup, conditions, spectra);


%% Show detailed output

% Show mic positions
showMicLayout(setup.micPos(1,:), setup.micPos(2,:));

% Show scan plane at various frequencies
if ~isfield(setup, 'fRange')
    setup.fRange = [500 4000];
end

fShow = round(linspace(setup.fRange(1), setup.fRange(2), nPlots));
for i=1:nPlots
    showBeamforming(spectra, setup, fShow(i));
end

% Show final spectrum
figure;
plot(spectra.f, spectra.SPL, 'b-', 'LineWidth', 2);
xlabel('f [Hz]');
ylabel('SPL [dB]');
box on;
set(gca, 'XScale', 'log')


%% Store result
storeResult = input('Save result? (y/n): ','s');
if strcmp(storeResult, 'y')
    
    if isfield(setup, 'name')
        outputFolder = ['inspect_' setup.name];
    else
        outputFolder = 'inspectSample';
    end
    fprintf('Storing results in %s\n', outputFolder);
    mkdir(outputFolder);
    
    % Save data
    save([outputFolder filesep 'sample.mat'], 'conditions', 'setup', 'spectra');
    
    % Save all figures
    figList = handle( sort( double(findall(0, 'type', 'figure') ) ));
    savefig(figList, [outputFolder filesep 'figures'], 'compact');
    for i=1:length(figList)
        set(figList(i),'WindowStyle','normal');
        set(figList(i), 'Units', 'pixels', 'Position', [10 10 800 600]);
        saveas(figList(i), [outputFolder filesep 'figure_' num2str(i) '.png']);
        set(figList(i),'WindowStyle','docked');
    end
    
end