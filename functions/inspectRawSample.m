function [setup, conditions, spectra] = inspectRawSample(setup)
%% Select file
dataFiles = dir('*.h5');

if length(dataFiles) < 1
    error('No data files found. Place main script in folder with data');
else
    dataFile = dataFiles(1);
    fprintf('Inspecting file %s \n', dataFile.name);
    filePath = [dataFile.folder filesep dataFile.name];
end

%% Settings
setup.dataPortion = input('Portion of data to use in this check (%): ')/100;
nPlots = input('Number of frequencies to plot: ');

%% Pass to beamforming
[setup, conditions, spectra] = runBeamforming(filePath, setup);

%% Determine integration window
setup = defineIntegrationWindow(setup);

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