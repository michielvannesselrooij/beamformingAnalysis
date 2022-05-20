function showSpectra(f, SPL, name)
% PLOTS ONE OR MORE BEAMFORMING SPECTRA WITH OR WITHOUT A LEGEND
%
% INPUTS
% -------
% f            (cell array of) 1D       Frequency range
% SPL          (cell array of) 1D       Sound pressure levels for f
% name         (cell array of) string   Legend names
% -------

%% Checks
if iscell(f)
    if ~iscell(SPL)
        error('if f is a cell array, then so should SPL');
    end
    
    if length(f) ~= length(SPL)
        error('Number of items in f and SPL should be equal');
    end
    
end

if exist('name', 'var')
    if iscell(name) && length(name) ~= length(f)
        error('Number of items in name, f, and SPL should be equal');
    elseif ~ischar(name) && ~iscell(f)
        error('Name should be a string');
    end
end

%% Prepare

if ~iscell(f)
    % Wrap for consistent behavior
    f    = {f};
    SPL  = {SPL};
    if exist('name', 'var')
        name = {name};
    end
end

N = length(f);

% Make colors
c0 = 0.6;
if iscell(f)
    c = [ repmat([c0, c0], N, 1), linspace(c0, 1, N)' ];
end

%% Plot
figure; 
hold on;
box on;

for i=1:N
    if exist('name', 'var')
        plot(f{i}, SPL{i}, '-', 'LineWidth', 2, 'Color', c(i,:),...
            'displayName', name{i});
    else
        plot(f{i}, SPL{i}, '-', 'LineWidth', 2, 'Color', c(i,:));
    end
end

set(gca, 'XScale', 'log')

xlabel('f [Hz]');
ylabel('SPL [dB]');

if exist('name', 'var')
    legend(name, 'Location', 'EastOutside');
end

