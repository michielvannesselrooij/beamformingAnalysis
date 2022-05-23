function beamformingSummary(setup, spectra, nPlots)
%% PLOTS BEAMFROMING RESULTS IN VARIOUS WAYS. OPTIONALLY SAVES THEM
%
% INPUT
% ------
% See requirements:
%  showMicLayout.m
%  showBeamforming.m
%  showSpectra.m
% ------

%% Show mic positions
showMicLayout(setup.micPos(1,:), setup.micPos(2,:));

%% Show scan plane at various frequencies
if isfield(setup, 'fRange')
    fRange = setup.fRange;
else
    fRange = [500 4000];
end

fShow = round(linspace(fRange(1), fRange(2), nPlots));
for i=1:nPlots
    showBeamforming(spectra, setup, fShow(i));
end

%% Show final spectrum
showSpectra(spectra.f, spectra.SPL);