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
    fprintf('  Using specified reference pressure: %.0E Pa \n', p_ref)
else
    p_ref = 20e-6;
end

% Scan plane information
if isfield(spectra, 'scanPlaneX')
    scanPlaneX = spectra.scanPlaneX;
else
    error('Missing coordinates: spectra.scanPlaneX"');
end

if isfield(spectra, 'scanPlaneY')
    scanPlaneY = spectra.scanPlaneY;
else
    error('Missing coordinates: spectra.scanPlaneY"');
end

% Integration plane
if isfield(setup, 'intPlane')
    intPlaneX = setup.intPlane([1,3]);
    intPlaneY = setup.intPlane([2,4]);
end

%% Prepare data to show

diagnostic = false;
abort = false;
if ~exist('f_plot', 'var')
    
    B_select = sum(spectra.B,3);
    SPL = 10*log10(B_select/p_ref^2);
    
    fprintf('No frequency specified for plotting. Using sum of all frequencies\n');
    
elseif isscalar(f_plot)
    
    % Look up nearest frequency available
    [~,idx] = min(abs(spectra.f - f_plot));
    
    f_plot_final = spectra.f(idx);
    B_select = spectra.B(:,:,idx);
    SPL = 10*log10(B_select/p_ref^2);
    
elseif ischar(f_plot) && strcmp(f_plot, 'diagnostic')
    diagnostic = true;
    
    if isfield(spectra, 'scanPlaneB')
        scanPlaneB = spectra.scanPlaneB;
    else
        warning(['Trying to show intermediate map,'....
            ' but missing "spectra.scanPlaneB". Aborting.']);
        abort = true;
    end
    
else
    error('Input "f_plot" should be a scalar, or left empty');
end
    

%% Create figure
dynRange = 10; % dB

if ~abort

    figure;
    hold on
    
    % Prepare contour data
    if diagnostic
        plotVal = scanPlaneB;
    else
        plotVal = SPL;
    end
    
    % Find max level to show
    if exist('intPlaneX', 'var')
        idx      = (scanPlaneX >= intPlaneX(1)) &...
                   (scanPlaneX <= intPlaneX(2)) &...
                   (scanPlaneY >= intPlaneY(1)) &...
                   (scanPlaneY <= intPlaneY(2));
               
        plotValWithin = plotVal .* repmat(idx, [1,1,size(plotVal,3)]);
        maxLevel = max(plotValWithin(:));
    else
        maxLevel = max(plotVal(:));
    end
    minLevel = maxLevel-dynRange;
    levels = linspace(minLevel, maxLevel, 20);
    
    % Show beamforming
    contourf(scanPlaneX, scanPlaneY, plotVal, levels, 'LineColor', 'none');

    % Show geometry (optional)
    if isfield(setup, 'wing')
        af_loc_x = [setup.wing.loc(1), setup.wing.loc(3)];
        af_loc_y = [setup.wing.loc(2), setup.wing.loc(4)];

        x = [af_loc_x(1), af_loc_x(1), af_loc_x(2), af_loc_x(2), af_loc_x(1)];
        y = [af_loc_y(1), af_loc_y(2), af_loc_y(2), af_loc_y(1), af_loc_y(1)];

        plot(x, y, 'w-', 'LineWidth', 3);
    end
    
    % Show integration window
    if exist('intPlaneX', 'var')
        x = [intPlaneX(1), intPlaneX(1), intPlaneX(2), intPlaneX(2), intPlaneX(1)];
        y = [intPlaneY(1), intPlaneY(2), intPlaneY(2), intPlaneY(1), intPlaneY(1)];
        plot(x, y, 'w--', 'LineWidth', 2);
    end

    % Format
    axis equal
    colormap copper
    caxis([minLevel, maxLevel]);
    
    hc = colorbar;
    title(hc, 'dB');
    xlabel('x [m]');
    ylabel('y [m]');

    if diagnostic
        title('Intermediate (g2) map');
    elseif exist('f_plot_final', 'var')
        title([num2str(f_plot_final) ' Hz']);
    end

    % Wrap up
    fprintf('\n');
end