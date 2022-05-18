function spectra = beamforming(dataFiles, conditions, setup)
tic;

%% Check inputs

% File selection
if ischar(dataFiles)
    
    % Single file
    dataFile1 = dataFiles;
    dataFiles = {dataFiles}; % Pack into cell
    
elseif iscell(dataFiles)
    
    % Multiple files
    dataFile1 = dataFiles{1};
    
else
    error('input "dataFiles" should be a string or a set of strings');
    
end

dataFile=1;

% Frequency range
if isfield(setup, 'fRange')
    fRange = setup.fRange;
else
    fRange = [500 4000];
end

if length(fRange) == 2
    fprintf('Frequency range: %i - %i Hz \n', fRange(1), fRange(2));
elseif length(fRange) == 1
    fprintf('Number of frequency bands: %i', fRange)
else
    error('Number of elements in frequency range "f_select" must be 2');
end


% Distance
if isfield(setup, 'h')
    h = setup.h;
else
    error('No distance to target specified in input "setup.h"');
end


% Temperature
if isfield(conditions, 'T')
    T = conditions.T;
else
    T = h5read(dataFile1,'/Conditions/Temperature');
end
fprintf('Temperature: %0.1f \n', T);


% Speed of sound
if isfield(conditions, 'c')
    c = conditions.c;
else
    error('No speed of sound specified in "conditions.c"');
end


% Mach number
if isfield(conditions, 'M_eff')
    M_eff = conditions.M_eff;
else
    error('No effective Mach number specified in "conditions.M_eff"');
end


% Microphone positions
if isfield(setup, 'micPos')
    micPos = setup.micPos(1:3, :);
    nMics = size(micPos,2);
    fprintf('Loaded %i microphones \n', nMics);
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
if isfield(setup, 'scanPlaneLimits')
    scanPlaneLimits = setup.scanPlaneLimits;
else
    error('No scan plane limits specified in "setup.scan_plane_limits"');
end


% Detailed settings ------------
fprintf('Detailed settings: \n');

    % Time chunk [s]
    if isfield(setup, 'timeChunk')
        timeChunk = setup.timeChunk;
    else
        timeChunk = 0.1;
    end
    fprintf('    Time chunk: %0.2f s \n', timeChunk);

    
    % Data overlap in fraction
    if isfield(setup, 'overlap')
        overlap = setup.overlap;
    else
        overlap = 0.5;
    end
    fprintf('    Data overlap fraction: %0.2f \n', overlap);

    
    % Resolution [m]
    if isfield(setup, 'scanPlaneResolution')
        scanPlaneResolution = setup.scanPlaneResolution;
    else
        scanPlaneResolution = 0.01;
    end
    fprintf('    Scan plane resolution: %0.2f m \n\n', scanPlaneResolution);
    
    
    % Diagonal removal of the CSM
    if isfield(setup, 'diagonalRemoval')
        diagonalRemoval = setup.diagonalRemoval;
    else
        diagonalRemoval = 0;
    end
    fprintf('    Diagonal removal: %i \n', diagonalRemoval);
    
    
    % Portion of data [%] to be used, centered at the half-time
    if isfield(setup, 'dataPortion')
        dataPortion = setup.dataPortion;
    else
        dataPortion = 1;
    end
    fprintf('    Portion of data used: %i \n', dataPortion*100);
        
    
    % Flow-corrected steering vector
    if isfield(setup, 'flowVector')
        flowCorrectedSteeringVector = 1;
        flowVector = setup.flowVector;
        
        fprintf('    Flow-corrected steering vector: \n');   
        disp(flowVector);
        
    else
        flowCorrectedSteeringVector = 0;
        fprintf('    Flow-corrected steering vector: no \n');   
    end
fprintf('\n');

%% Load data file
fprintf('Loading data...\n');

% Pre-allocate
N = length(dataFiles);
dataSize = zeros(N,1);
micCount = zeros(N,1);

for i=1:N
    fileInfo = h5info(dataFiles{i},'/AcousticData/Measurement');
    dataSize(i) = fileInfo.ChunkSize(1);
    micCount(i) = fileInfo.ChunkSize(2);
end

if numel(unique(micCount)) > 1
    error('Number of microphones not consistent between data files in set ');
end

dataMic = zeros(sum(dataSize), micCount(1));

% Fill data
IdxSteps = [1; dataSize];
for i=1:N
    idx1 = sum(IdxSteps(1:i));
    idx2 = idx1 + IdxSteps(i+1)-1;
    dataMic(idx1:idx2, :) = h5read(dataFiles{i}, '/AcousticData/Measurement');
end

% Filter and sort mics
data          = dataMic(:,1:nMics);

% Sampling frequency [Hz]
if isfield(setup, 'fs')
    fs = setup.fs;
    fprintf('Explicitly setting acquisition frequency: %i \n', fs);
else
    fs = h5read(dataFile1,'/Acquisition/SampleRate');
end

%% Environment preparation
fprintf('Extracting data...\n');

% Generate frequency bands
if length(fRange) == 1
    f_cen   = 10.^(fRange./10);
    f_upper = 2^(1/6).*f_cen;
    f_lower = f_upper./2^(1/3);
else
    f_lower = fRange(1);
    f_upper = fRange(2);
end

% Microphones and data preparation
data                    = data(...
                                max(1,round(size(data,1)/2 - (0.5*dataPortion)*size(data,1))) : ...
                                round(size(data,1)/2       + (0.5*dataPortion)*size(data,1)),...
                                :);
data(:,removeMics)      = [];
micPos(3,:)             = 0;
micPos(:,removeMics)    = [];
nMics = length(micPos);

% Windowing
ind_start_chunk         = 1;
chunk_length            = floor( timeChunk*fs );
overlap_chunk_length    = floor( overlap*chunk_length );
array_us                = ind_start_chunk :  chunk_length-overlap_chunk_length : size(data,1)-overlap_chunk_length;
usable_signal           = [ array_us(1:end-1); array_us(1:end-1)+chunk_length-1 ];
   
%% Develop CSM
fprintf('Developing CSM...\n');

% f axis for small chunks
df_chunk = fs/chunk_length;
f_chunk  = (0:chunk_length/2-1)*df_chunk;

% indices of frequency to be taken from fft of small chunks
ind_f_lower = floor(f_lower.*chunk_length./fs + 1);
ind_f_upper = floor(f_upper.*chunk_length./fs + 1);
data_chunk  = zeros(chunk_length,nMics,size(usable_signal,2));

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

C_av=zeros(nMics^2,ind_f_upper-ind_f_lower+1);

for j = 1:size(X_chunk_ff,3)
    CC=rectpulse(conj(X_chunk_ff(:,:,j)),size(X_chunk_ff,1)).*repmat(X_chunk_ff(:,:,j),[size(X_chunk_ff,1),1]);
    C_av(:,j)=mean(real(CC)-1i*imag(CC),2);
end

C_averaged=reshape(C_av,[nMics,nMics,ind_f_upper-ind_f_lower+1]);

if diagonalRemoval == 1
    for i = 1:size(C_averaged,3)
        C_averaged(:,:,i) = C_averaged(:,:,i) - diag(diag(C_averaged(:,:,i)));
    end
end

C_averaged = C_averaged.*df_chunk;

%% Build scan plane

x_steps = scanPlaneLimits(1) : scanPlaneResolution : scanPlaneLimits(2);
y_steps = scanPlaneLimits(3) : scanPlaneResolution : scanPlaneLimits(4);

scan_plane_x = repmat(x_steps,length(y_steps),1); 
scan_plane_y = repmat(flipud(y_steps'),1,length(x_steps));
scan_plane_z = h.*ones(size(scan_plane_x));

%% Beamforming
fprintf('Calculating radii...\n');

if flowCorrectedSteeringVector == 1
    
    % Steering vector for flow correction
    M_vec           = M_eff .* flowVector;
    beta_square     = 1-M_eff^2;
    psi_j_vec       = [scan_plane_x(:), scan_plane_y(:), -scan_plane_z(:)];
    x_n_vec         = micPos';
    r_vec           = repmat(x_n_vec, size(psi_j_vec, 1), 1) - rectpulse(psi_j_vec, size(x_n_vec, 1));
    radii0          = sqrt(sum(r_vec.^2, 2));
    
    M_vec_mat       = repmat(M_vec,size(r_vec,1),1);
    delt_all_m      = (-dot(M_vec_mat,r_vec,2) + sqrt((dot(M_vec_mat,r_vec,2)).^2 + beta_square*radii0.^2)) ./ (c*beta_square);
    delt_all        = permute(reshape(delt_all_m, [nMics,size(scan_plane_x)]), [2,3,1]);
    radii           = permute(reshape(radii0, [nMics, size(scan_plane_x)]), [2,3,1]);

else
    
    % Basic steering vector, ignore the flow
    radii = sqrt((repmat(scan_plane_x,[1,1,nMics]) - repmat(reshape(micPos(1,:),[1,1,nMics]), [size(scan_plane_x),1])).^2 ...
               + (repmat(scan_plane_y,[1,1,nMics]) - repmat(reshape(micPos(2,:),[1,1,nMics]), [size(scan_plane_y),1])).^2 ...
               + (repmat(scan_plane_z,[1,1,nMics]) - repmat(reshape(micPos(3,:),[1,1,nMics]), [size(scan_plane_x),1])).^2);

    delt_all = radii./c;
    
end

radius = repmat(radii, [1,1,1]);
delt   = repmat(delt_all, [1,1,1]);

fprintf('Beamforming...\n');

B_ind  = zeros([numel(scan_plane_x), length(f_ind_c_averaged)]);
B2_ind = zeros([numel(scan_plane_x), length(f_ind_c_averaged)]);
bar_bf = waitbar(0, 'CFDBF: Performing beamforming ...');

for f_index = 1:length(f_ind_c_averaged)
    waitbar(f_index/length(f_ind_c_averaged));
    
    g     = ((-exp(-2*pi*1i*f_ind_c_averaged(f_index).*delt)) ./ radius);
    gsel  = reshape(g, [numel(scan_plane_x), nMics]);
    C_ind = C_averaged(:, :, f_index);
    
    for s = 1:size(gsel,1)
        Norm1 = norm(conj(gsel(s,:))*C_ind.'*gsel(s,:).');
        gnorm = norm(gsel(s,:));
        B_ind(s,f_index)  = Norm1/gnorm^4;
        B2_ind(s,f_index) = Norm1/gnorm^2;
    end
end

close(bar_bf);

B = reshape(B_ind, [size(scan_plane_x), length(f_ind_c_averaged)]);

%% Structure output

spectra.B                  = B;
spectra.f                  = f_ind_c_averaged;
spectra.scanPlaneX         = scan_plane_x;
spectra.scanPlaneY         = scan_plane_y;

%% Round up
fprintf('\nDone. It took %0.0f seconds \n\n', toc)
