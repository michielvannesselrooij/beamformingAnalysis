function spectra = beamforming_SPI(dataFile, conditions, setup)
% v7.0 - Stripped to beamforming core task (MvN - 2022-05-12)
tic;

%% Check inputs

% Frequency range
if isfield(setup, 'f_select')
    f_select = setup.f_select;
else
    f_select = [500 4000];
end
fprintf('Frequency range: %i - %i Hz \n', f_select(1), f_select(2));

% Distance
if isfield(setup, 'h')
    h = setup.h;
else
    error('No distance to target specified in input "setup.h"');
end

% Sampling frequency [Hz]
if isfield(setup, 'fs')
    fs = setup.fs;
else
    error('No sampling frequency specified in input "setup.fs"');
end

% Temperature
if isfield(conditions, 'T')
    T = conditions.T;
else
    error('No temperature specified in input "setup.T"');
end

% Mach number
if isfield(conditions, 'M_eff')
    M_eff = conditions_M_eff;
else
    error('No effective Mach number specified in "conditions.M_eff"');
end

% Microphone positions
if isfield(setup, 'mic_poses')
    mic_poses = setup.mic_poses;
    n_mics = size(mic_poses,2);
    fprintf('Loaded %i microphones \n', n_mics);
else
    error('No mic positions specified in "setup.mic_poses"');
end

% Indices of mics to be removed
if isfield(setup, 'brokenMics')
    removeMics = setup.brokenMics;
    fprintf('Ignoring following mics: \n');
    disp(removeMics);
else
    fprintf('Using all mics \n');
end

% Scan plane limits
if isfield(setup, 'scan_plane_limits')
    scan_plane_limits = setup.scan_plane_limits;
else
    error('No scan plane limits specified in "setup.scan_plane_limits"');
end

% Detailed settings ------------
fprintf('Detailed settings: \n');

    % Time chunk [s]
    if isfield(setup, 'time_chunk')
        time_chunk = setup.time_chunk;
    else
        time_chunk = 0.1;
    end
    fprintf('    Time chunk: 0.2%f s \n', time_chunk);

    % Data overlap in fraction
    if isfield(setup, 'overlap')
        overlap = setup.overlap;
    else
        overlap = 0.5;
    end
    fprintf('    Data overlap fraction: 0.2%f \n', overlap);

    % Resolution [m]
    if isfield(setup, 'scan_plane_resolution')
        scan_plane_resolution = setup.scan_plane_resolution;
    else
        scan_plane_resolution = 0.01;
    end
    fprintf('    Scan plane resolution: 0.3%f m \n\n', scan_plane_resolution);
    
    % Diagonal removal of the CSM
    if isfield(setup, 'diagonal_removal')
        diagonal_removal = setup.diagonal_removal;
    else
        diagonal_removal = 0;
    end
    fprintf('    Diagonal removal: %i \n', diagonal_removal);
    
    % Portion of data [%] to be used, centered at the half-time
    if isfield(setup, 'data_portion')
        data_portion = setup.data_portion;
    else
        data_portion = 1;
    end
    fprintf('    Portion of data used: %i \n', data_portion*100);
        
    % Flow-corrected steering vector
    if isfield(setup, 'flow_dir_vector')
        flow_corrected_steering_vector = 1;
        flow_dir_vector = setup.flow_dir_vector;
        
        fprintf('    Flow-corrected steering vector: \n');   
        disp(flow_corrected_steering_vector);
        
    else
        flow_corrected_steering_vector = 0;
        fprintf('    Flow-corrected steering vector: no \n');   
    end
fprintf('\n');

%% Load data file
fprintf('Loading data...\n');

datafileload  = strcat(dataFile.path, dataFile.fileName);
dataMic       = h5read(datafileload,'/AcousticData/Measurement');
data          = dataMic(:,1:n_mics);

%% Environment preparation
fprintf('Beamforming...\n');

% Speed of sound [m/s]
c     = 331.5 + (0.6*T);         

% Generate frequency bands, format f(1,number of bands)
if length(f_select) == 1
    f_cen = 10.^(f_select./10);
    f_upper = 2^(1/6).*f_cen;
    f_lower = f_upper./2^(1/3);
else
    f_lower = f_select(1);
    f_upper = f_select(2);
end

% Microphones and data preparation
data                    = data(...
                                max(1,round(size(data,1)/2-(0.5*data_portion)*size(data,1))) : ...
                                round(size(data,1)/2+(0.5*data_portion)*size(data,1)),...
                                :);
data(:,removeMics)      = [];
mic_poses(3,:)          = 0;
mic_poses(:,removeMics) = [];

% Windowing
ind_start_chunk         = 1;
chunk_length            = floor( time_chunk*fs );
overlap_chunk_length    = floor( overlap*chunk_length );
array_us                = ind_start_chunk :  chunk_length-overlap_chunk_length : size(data,1)-overlap_chunk_length;
usable_signal           = [ array_us(1:end-1); array_us(1:end-1)+chunk_length-1 ];
   
%% Develop CSM
% f axis for small chunks
df_chunk = fs/chunk_length;
f_chunk  = (0:chunk_length/2-1)*df_chunk;

% indices of frequency to be taken from fft of small chunks
ind_f_lower = floor(f_lower.*chunk_length./fs + 1);
ind_f_upper = floor(f_upper.*chunk_length./fs + 1);
data_chunk  = zeros(chunk_length,n_mics,size(usable_signal,2));

for k = 1:size(usable_signal,2)
    data_chunk(:,:,k) = data(usable_signal(1,k):usable_signal(2,k),:)...
                        - repmat(mean(data(usable_signal(1,k):usable_signal(2,k),:),1),chunk_length,1);
end

win_hann        = repmat(hann(chunk_length),[1,size(data,2),size(usable_signal,2)]);
data_chunk_hann = data_chunk.*win_hann;
    
% With Hanning window + overall sound pressure level (dB)
X_chunk_hann              = fft(data_chunk_hann);
X_chunk_hann              = X_chunk_hann(1:length(f_chunk),:,:)/(sqrt(sum(hann(chunk_length).^2)));
X_chunk_hann(2:end-1,:,:) = sqrt(2)*X_chunk_hann(2:end-1,:,:);
X_chunk_hann_collect      = X_chunk_hann;

% Averaging chunks
X_chunk_ff = permute(sqrt(1/(chunk_length*df_chunk)).*X_chunk_hann_collect(ind_f_lower:ind_f_upper,:,:),[2,3,1]);
f_ind_c_averaged = f_chunk(ind_f_lower:ind_f_upper);

C_av=zeros(n_mics^2,ind_f_upper-ind_f_lower+1);

for j = 1:size(X_chunk_ff,3)
    CC=rectpulse(conj(X_chunk_ff(:,:,j)),size(X_chunk_ff,1)).*repmat(X_chunk_ff(:,:,j),[size(X_chunk_ff,1),1]);
    C_av(:,j)=mean(real(CC)-1i*imag(CC),2);
end

C_averaged=reshape(C_av,[n_mics,n_mics,ind_f_upper-ind_f_lower+1]);

if diagonal_removal == 1
    for i = 1:size(C_averaged,3)
        C_averaged(:,:,i) = C_averaged(:,:,i) - diag(diag(C_averaged(:,:,i)));
    end
end

C_averaged = C_averaged.*df_chunk;

%% Build scan plane

x_steps = scan_plane_limits(1) : scan_plane_resolution : scan_plane_limits(2);
y_steps = scan_plane_limits(3) : scan_plane_resolution : scan_plane_limits(4);

scan_plane_x = repmat(x_steps,length(y_steps),1); 
scan_plane_y = repmat(flipud(y_steps'),1,length(x_steps));
scan_plane_z = h.*ones(size(scan_plane_x));

%% Beamforming

if flow_corrected_steering_vector == 1
    
    % Steering vector for flow correction
    M_vec           = M_eff.*flow_dir_vector;
    beta_square     = 1-M_eff^2;
    psi_j_vec       = [scan_plane_x(:), scan_plane_y(:), -scan_plane_z(:)];
    x_n_vec         = mic_poses';
    r_vec           = repmat(x_n_vec,size(psi_j_vec,1),1)-rectpulse(psi_j_vec,size(x_n_vec,1));
    radii           = sqrt(sum(r_vec.^2,2));
    radii           = permute(reshape(radii,[n_mics,size(scan_plane_x)]),[2,3,1]);
    
    M_vec_mat       = repmat(M_vec,size(r_vec,1),1);
    delt_all_m      = (-dot(M_vec_mat,r_vec,2)+sqrt((dot(M_vec_mat,r_vec,2)).^2+beta_square*radii.^2))./(c*beta_square);
    delt_all        = permute(reshape(delt_all_m,[n_mics,size(scan_plane_x)]),[2,3,1]);

else
    
    % Basic steering vector, ignore the flow
    radii = sqrt((repmat(scan_plane_x,[1,1,n_mics]) - repmat(reshape(mic_poses(1,:),[1,1,n_mics]), [size(scan_plane_x),1])).^2 ...
               + (repmat(scan_plane_y,[1,1,n_mics]) - repmat(reshape(mic_poses(2,:),[1,1,n_mics]), [size(scan_plane_y),1])).^2 ...
               + (repmat(scan_plane_z,[1,1,n_mics]) - repmat(reshape(mic_poses(3,:),[1,1,n_mics]), [size(scan_plane_x),1])).^2);

    delt_all = radii./c;
    
end

radius = repmat(radii, [1,1,1]);
delt   = repmat(delt_all, [1,1,1]);

B_ind  = zeros([numel(scan_plane_x), length(f_ind_c_averaged)]);
B2_ind = zeros([numel(scan_plane_x), length(f_ind_c_averaged)]);
bar_bf = waitbar(0, 'CFDBF: Performing beamforming ...');

for f_index = 1:length(f_ind_c_averaged)
    waitbar(f_index/length(f_ind_c_averaged));
    
    g     = ((-exp(-2*pi*1i*f_ind_c_averaged(f_index).*delt)) ./ radius);
    gsel  = reshape(g, [numel(scan_plane_x), length(mic_poses)]);
    C_ind = C_averaged(:, :, f_index);
    
    for s = 1:size(gsel,1)
        Norm1 = norm(conj(gsel(s,:))*C_ind.'*gsel(s,:).');
        gnorm = norm(gsel(s,:));
        B_ind(s,f_index)  = Norm1/gnorm^4;
        B2_ind(s,f_index) = Norm1/gnorm^2;
    end
end

close(bar_bf);

B = reshape(B_ind,[size(scan_plane_x),length(f_ind_c_averaged)]);

%% Structure output

spectra.B                          = B;
spectra.f_ind_c_averaged           = f_ind_c_averaged;
spectra.scan_plane_x               = scan_plane_x;
spectra.scan_plane_y               = scan_plane_y;

%% Round up

fprintf('\n')
toc;
fprintf('\n\n')
