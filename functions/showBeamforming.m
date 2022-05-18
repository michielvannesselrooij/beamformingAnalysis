function showBeamforming(spectra, setup, f_plot)
fprintf('Making beamforming plot \n')

%% Check inputs

% Beamforming field
if ~isfield(spectra, 'B')
    error('Nothing to show: "spectra.B" not defined');
end

% Reference pressure
if isfield(setup, 'p_ref')
    p_ref = setup.p_ref;
else
    p_ref = 20e-6;
    fprintf('  Using default reference pressure: %.0E Pa \n', p_ref)
end

% Scan plane information
if isfield(spectra, 'scan_plane_x')
    scan_plane_x = spectra.scan_plane_x;
else
    error('Missing coordinates: spectra.scan_plane_x"');
end

if isfield(spectra, 'scan_plane_y')
    scan_plane_y = spectra.scan_plane_y;
else
    error('Missing coordinates: spectra.scan_plane_y"');
end

%% Prepare data to show
if ~exist('f_plot', 'var')
    B_select = sum(spectra.B,3);
    
elseif isscalar(f_plot)
    
    % Look up nearest frequency available
    [~,idx] = min(abs(spectra.f_ind_c_averaged - f_plot));
    
    f_plot_final = spectra.f_ind_c_averaged(idx);
    B_select = spectra.B(:,:,idx);
    
else
    error('Input "f_plot" should be a scalar, or left empty');
end
    
SPL = 10*log10(B_select/p_ref^2);

%% Create figure
figure;

% Show beamforming
contourf(spectra.scan_plane_x, spectra.scan_plane_y, SPL, 20, 'LineColor', 'none');

% Show geometry (optional)
hold on
if isfield(setup, 'wing')
    af_loc_x = [setup.wing.loc(1), setup.wing.loc(3)];
    af_loc_y = [setup.wing.loc(2), setup.wing.loc(4)];

    x = [af_loc_x(1), af_loc_x(1), af_loc_x(2), af_loc_x(2), af_loc_x(1)];
    y = [af_loc_y(1), af_loc_y(2), af_loc_y(2), af_loc_y(1), af_loc_y(1)];

    plot(x, y, 'w-', 'LineWidth', 3);
end

% Show integration window
intPlaneX = setup.intPlane([1,3]);
intPlaneY = setup.intPlane([2,4]);
x = [intPlaneX(1), intPlaneX(1), intPlaneX(2), intPlaneX(2), intPlaneX(1)];
y = [intPlaneY(1), intPlaneY(2), intPlaneY(2), intPlaneY(1), intPlaneY(1)];
plot(x, y, 'w--', 'LineWidth', 2);

% Format
colormap copper
hc = colorbar;
title(hc, 'dB');
axis equal
xlabel('x [m]');
ylabel('y [m]');

if exist('f_plot_final', 'var')
    title([num2str(f_plot_final) ' Hz']);
end

%% Wrap up
fprintf('\n');