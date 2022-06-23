function SPL = selectIntegrationWindow(setup, conditions, spectra)
% SOURCES SPECIFIED IN "SETUP" AND COMPLETES "SETUP" BASED ON CONDITIONS
%
% INPUTS
% -------
% conditions.c          double  Speed of sound                          [m/s]
% spectra.f             1D      Frequencies                             [Hz]
% spectra.B             1D      Sound pressure levels for frequencies f [dB]
% spectra.scanPlaneX    2D      x-positions of total beamforming window [m]
% spectra.scanPlaneY    2D      y-positions of total beamforming window [m]
% setup.h               double  Distance target from center microphone  [m]
% setup.intPlane        4x1     Integration window: Xmin Xmax Ymin Ymax [m]    
% setup.micPos          4xN     positions + ID of microphones           [m, -]
%
% (optional)
% setup.p_ref           double  Reference pressure for calculating dB   [Pa]
% setup.resolution      double  Window resolution                       [m]
% -------

%% Inputs
fprintf('Integrating selected window...\n');

c          = conditions.c;

f          = spectra.f;
B          = spectra.B;
scanPlaneX = spectra.scanPlaneX;
scanPlaneY = spectra.scanPlaneY;

h          = setup.h;
intPlane   = setup.intPlane;
micPos     = setup.micPos(1:3,:);

%% Checks

if size(micPos,1) < 3
    error('Input micPos should have 3 at least rows (x;y;z)');
end

if length(f) ~= size(B,3)
    error('Inputs B and f should be vectors of the same length');
end

if size(scanPlaneX) ~= size(scanPlaneY)
    error(['Inputs scanPlaneX and scanPlaneY should be 2D matrices '...
        'of the same length']);
end


%% Settings

% Dynamic range [+dB] considered in the ROI of SPI
dynRange = 6;

% Reference pressure
if isfield(setup, 'p_ref')
    p_ref = setup.p_ref;
    fprintf('  Using specified reference pressure: %.0E Pa \n', p_ref)
else
    p_ref = 20e-6;
end

% Resolution
if isfield(setup, 'resolution')
    resolution = setup.resolution;
else
    resolution = scanPlaneX(1,2) - scanPlaneX(1,1);
end

%% Apply source power integration

% Numerator, sum the source powers from the experiment
idx      = (scanPlaneX >= intPlane(1)) &...
           (scanPlaneX <= intPlane(3)) &...
           (scanPlaneY >= intPlane(2)) &...
           (scanPlaneY <= intPlane(4));
       
B_within = B .* repmat(idx, [1,1,size(B,3)]);


minLevel = repmat( 20*log10(sqrt(max(max(B_within))) / (h*p_ref)) - dynRange,...
             [size(B,1), size(B,2),1]);
         
B_within(20*log10(sqrt(B_within)/(h*p_ref)) < minLevel) = 0;

integrated_B_f = squeeze(sum(sum(B_within)));

% Get PSF of the simulated source in the integration area
% Build scan grid for the integration
x_steps_int = intPlane(1) : resolution : intPlane(3);
y_steps_int = intPlane(2) : resolution : intPlane(4);

scan_plane_x_int = repmat(x_steps_int, length(y_steps_int), 1);
scan_plane_y_int = repmat(flipud(y_steps_int'), 1, length(x_steps_int));
scan_plane_z_int = h .* ones(size(scan_plane_x_int));

% find radii for the integration region
radii_int    = sqrt((repmat(scan_plane_x_int, [1,1,size(micPos,2)]) - repmat(reshape(micPos(1,:), [1,1,size(micPos,2)]), [size(scan_plane_x_int), 1])).^2+...
                 (repmat(scan_plane_y_int, [1,1,size(micPos,2)]) - repmat(reshape(micPos(2,:), [1,1,size(micPos,2)]), [size(scan_plane_y_int), 1])).^2+...
                 (repmat(scan_plane_z_int, [1,1,size(micPos,2)]) - repmat(reshape(micPos(3,:), [1,1,size(micPos,2)]), [size(scan_plane_x_int), 1])).^2);
delt_all_int = radii_int ./ c;


% Beamforming of the point source's PSF
Ix      = rectpulse(mean([intPlane(1) intPlane(3)]),size(micPos,2))';
Iy      = rectpulse(mean([intPlane(2) intPlane(4)]),size(micPos,2))';
Iz      = rectpulse(h,size(micPos,2))';
radius  = sqrt( (Ix-micPos(1,:)).^2+(Iy-micPos(2,:)).^2+(Iz-micPos(3,:)).^2);
delt    = radius / c;
PSF     = zeros([size(radii_int,1)*size(radii_int,2), length(f)]);

for i = 1:length(f)
    fk = f(i);
    % Denominator, sum the normalized PSF of a point source
    % make the steering vector to the center of the integration area, g_z
    % assuming the source is a point source in the center of the
    % integration area (ROI)
    gz = -exp(-2*pi*1i*fk*delt)./radius;
    gk = reshape(-exp(-2*pi*1i*fk*delt_all_int) ./ radii_int,...
                 [size(radii_int,1)*size(radii_int,2), size(radii_int,3)]);
             
    for j = 1:size(gk,1)
        PSF(j,i) = norm(conj(gk(j,:)) * (gz.'*conj(gz)) * gk(j,:).') / (norm(gk(j,:)))^4;
    end
end

PSF         = reshape(PSF, [size(radii_int,1), size(radii_int,2), length(f)]);

% Sum the PSF within the ROI (exclude PSF below predefined level)
minLevel    = repmat(20*log10(sqrt(max(max(PSF)))) - dynRange, [size(PSF,1), size(PSF,2), 1]);
PSF_curr    = PSF;
PSF_curr( 20*log10(sqrt(PSF_curr)) < minLevel ) = 0;
PSF_sum_f   = squeeze(sum(sum(PSF_curr)));

% P_exp per frequency
SPL         = 20*log10(sqrt((integrated_B_f./PSF_sum_f))./(h*p_ref));

