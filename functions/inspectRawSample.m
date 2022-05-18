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
nPlots = input('Number of frequencies to plot:            ');


%% Pass to beamforming
[setup, conditions, spectra] = runBeamforming(filePath, setup);


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