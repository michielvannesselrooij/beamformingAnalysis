function SPL_exp = selectIntegrationWindow(setup, conditions, spectra)
%% Inputs
f         = spectra.f_ind_c_averaged;
B         = spectra.B;
h         = setup.h;
c_sound   = 331.5 + (0.6*conditions.T); % speed of sound [m/s]
int_plane = setup.int_plane;


%% Settings
scan_plane_limits       = setup.scan_plane_limits;
scan_plane_resolution   = setup.scan_plane_resolution;
pref                    = 20e-6;
integration_lower_bound = 6; % Dynamic range [+dB] considered in the ROI of SPI
load mic_poses.mat 

%% Prepare

% Build scan grid
x_steps      = scan_plane_limits(1) : scan_plane_resolution : scan_plane_limits(2);
y_steps      = scan_plane_limits(3) : scan_plane_resolution : scan_plane_limits(4);
scan_plane_x = repmat(x_steps, length(y_steps), 1); 
scan_plane_y = repmat(flipud(y_steps'), 1, length(x_steps));
scan_plane_z = h .* ones(size(scan_plane_x));


%% Apply source power integration

% Numerator, sum the source powers from the experiment
idx      = (scan_plane_x >= int_plane(1)) &...
           (scan_plane_x <= int_plane(3)) &...
           (scan_plane_y >= int_plane(2)) &...
           (scan_plane_y <= int_plane(4));
       
B_within = B .* repmat(idx, [1,1,size(B,3)]);

Lcond    = repmat( 20*log10(sqrt(max(max(B_within))) / (h*pref)) - integration_lower_bound,...
             [size(B,1), size(B,2),1]);
         
B_within(20*log10(sqrt(B_within)/(h*pref)) < Lcond) = 0;

integrated_B_f = squeeze(sum(sum(B_within)));

% Get PSF of the simulated source in the integration area
% Build scan grid for the integration
x_steps_int = int_plane(1) : scan_plane_resolution : int_plane(3);
y_steps_int = int_plane(2) : scan_plane_resolution : int_plane(4);

scan_plane_x_int = repmat(x_steps_int, length(y_steps_int), 1);
scan_plane_y_int = repmat(flipud(y_steps_int'), 1, length(x_steps_int));
scan_plane_z_int = h .* ones(size(scan_plane_x_int));

% find radii for the integration region
% radii_int = zeros(size(scan_plane_x_int,1),size(scan_plane_x_int,2),size(mic_poses,2));

radii_int = sqrt((repmat(scan_plane_x_int, [1,1,size(mic_poses,2)]) - repmat(reshape(mic_poses(1,:), [1,1,size(mic_poses,2)]), [size(scan_plane_x_int), 1])).^2+...
                 (repmat(scan_plane_y_int, [1,1,size(mic_poses,2)]) - repmat(reshape(mic_poses(2,:), [1,1,size(mic_poses,2)]), [size(scan_plane_y_int), 1])).^2+...
                 (repmat(scan_plane_z_int, [1,1,size(mic_poses,2)]) - repmat(reshape(mic_poses(3,:), [1,1,size(mic_poses,2)]), [size(scan_plane_x_int), 1])).^2);
delt_all_int = radii_int ./ c_sound;


% Beamforming of the point source's PSF
Ix     = rectpulse(mean([int_plane(1) int_plane(3)]),size(mic_poses,2))';
Iy     = rectpulse(mean([int_plane(2) int_plane(4)]),size(mic_poses,2))';
Iz     = rectpulse(h,size(mic_poses,2))';
radius = sqrt( (Ix-mic_poses(1,:)).^2+(Iy-mic_poses(2,:)).^2+(Iz-mic_poses(3,:)).^2);
delt   = radius / c_sound;
PSF    = zeros([size(radii_int,1)*size(radii_int,2), length(f)]);

for i = 1:length(f)
    fk = f(i);
    % Denominator, sum the normalized PSF of a point source
    % make the steering vector to the center of the integration area, g_z
    % assuming the source is a point source in the center of the
    % integration area (ROI)
    gz = -exp(-2*pi*1i*fk*delt)./radius;
    gk = reshape(-exp(-2*pi*1i*fk*delt_all_int)./radii_int,[size(radii_int,1)*size(radii_int,2),size(radii_int,3)]);
    for j = 1:size(gk,1)
        PSF(j,i) = norm(conj(gk(j,:)) * (gz.'*conj(gz)) * gk(j,:).') / (norm(gk(j,:)))^4;
    end
end

PSF = reshape(PSF, [size(radii_int,1), size(radii_int,2), length(f)]);

% sum the PSF within the ROI,
% but exclude the PSF lower than a predefined level
Lcond2 = repmat(20*log10(sqrt(max(max(PSF)))) - integration_lower_bound, [size(PSF,1), size(PSF,2), 1]);
PSF_curr = PSF;
PSF_curr(20*log10(sqrt(PSF_curr))<Lcond2)=0;
PSF_sum_f = squeeze(sum(sum(PSF_curr)));

% P_exp per frequency
SPL_exp = 20*log10(sqrt((integrated_B_f./PSF_sum_f))./(h*pref));

