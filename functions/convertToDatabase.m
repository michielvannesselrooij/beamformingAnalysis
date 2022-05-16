clc; close all; clear;

%% Setting

outputFolder = 'noise-files';

%% Output

dataFolder = 'data-DU93';
dataFiles = dir([dataFolder filesep '*.mat']);

offset    = 161;

% Find file sets
name = {};
for i=1:length(dataFiles)
    
    idx = strfind(dataFiles(i).name, '_AoA');
    str = dataFiles(i).name(1:idx-1);
    
    if sum(contains(name, str)) == 0
        name = [name; str];
    end
    
end
N = length(name);

mkdir(outputFolder);
for i=1:N
    
    disp(['Processing ' name{i}]);
    
    % Start file
    fid = fopen([outputFolder filesep 'noise-' num2str(i+offset) '.csv'],'w');
    fprintf(fid, '%s\n', 'rho,T,V,alpha,alpha_eff,f,dB');
    
    % Find all data sets for this configuration
    setFiles = dir([dataFolder filesep name{i} '*.mat']);
    
    for j=1:length(setFiles)
        
        disp(['  ' setFiles(j).name]);
        
        % Extract data
        load([setFiles(j).folder filesep setFiles(j).name]);
        
        % Write file
        V           = conditions.U;
        alpha       = conditions.AoA;
        alpha_eff   = conditions.AoA * 2/3; % Assumption!
        f           = spectra.f_ind_c_averaged;
        SPL         = spectra.SPL_exp;
        fprintf(fid, '%s%f%s%f%s%f%s%f%s%f\n', ',,', V, ',', alpha, ',', alpha_eff, ',', f(1), ',', SPL(1));

        for ii=2:length(f)
            if ~isnan(f(ii)) && ~isnan(SPL(ii))
                fprintf(fid, '%s%f%s%f\n', ',,,,,', f(ii), ',', SPL(ii));
            end
        end
            
    end 
    
    % Close file
    fclose(fid);

end

disp(['Data output ready!']);