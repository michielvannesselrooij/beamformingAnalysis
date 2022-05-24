function showIntermediateMap(setup, spectra)
% PLOTS ONE OR MORE INTERMEDIATE BEAMFORMING MAPS
%
% INPUTS
% -------
% f            (cell array of) 1D       Frequency range
% SPL          (cell array of) 1D       Sound pressure levels for f
% name         (cell array of) string   Legend names
% -------

%% Prepare data

pref = 20e-6;

% Make g2 map
B_g2  = reshape(B2_ind, [size(scan_plane_x),length(f_ind_c_averaged)]);
B_sum = sum(B_g2,3);

% Convert to dB
scanPlaneB = 20*log10(sqrt(B_sum)/(h*pref));

% Normalise
scanPlaneB = scanPlaneB - max(max(scanPlaneB));


%% Plot 
figure;
hold on;

set(gca,'FontSize',24);
set(gca,'LineWidth',2.5)
set(gca,'fontname','times');
set(gcf,'color','w');

surf(scan_plane_x, scan_plane_y, scanPlaneB, 'EdgeColor', 'none');

hc = colorbar('EastOutside');
hc.Label.String = 'SPL [dB]';
hc.LineWidth = 2;
colormap copper

caxis([max(max(scanPlaneB))-dynamic_range, max(max(scanPlaneB))])

xlabel('x [m]')
ylabel('y [m]')
title({strcat('CFDBF at ',num2str(f_lower/1000),' - ',num2str(f_upper/1000),' kHz. Intermediate map before SPI')})

grid on
grid minor
view(2);
pbaspect([scan_plane_limits(2)-scan_plane_limits(1) scan_plane_limits(4)-scan_plane_limits(3) 1])