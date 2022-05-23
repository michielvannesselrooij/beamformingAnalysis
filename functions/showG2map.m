function showG2map(setup, spectra)

% Make g2 map
B_g2  = reshape(B2_ind, [size(scanPlaneX),length(f)]);
B_sum = sum(B_g2,3);

% Convert to dB
scanPlaneB = 20*log10(sqrt(B_sum)/(h*pref));

% Normalise
scanPlaneB = scanPlaneB - max(max(scanPlaneB));


%% Plot 
set(gca,'FontSize',24); set(gca,'LineWidth',2.5); hold on;
set(gca,'fontname','times');  set(gcf,'color','w');
surf(scan_plane_x-x_TE,scan_plane_y-y_TE,scan_plane_B, 'EdgeColor', 'none');
col_bar = colorbar('eastoutside');
caxis([max(max(scan_plane_B))-dynamic_range, max(max(scan_plane_B))])
col_bar.Label.String = 'SPL [dB]';
col_bar.LineWidth = 2;
colormap jet
hold on    
alpha 0.75
xlabel('x [m]')
ylabel('y [m]')
title({strcat('CFDBF at ',num2str(f_lower/1000),' - ',num2str(f_upper/1000),' kHz. Intermediate map before SPI')})

grid on
grid minor
view(2);
pbaspect([scan_plane_limits(2)-scan_plane_limits(1) scan_plane_limits(4)-scan_plane_limits(3) 1])