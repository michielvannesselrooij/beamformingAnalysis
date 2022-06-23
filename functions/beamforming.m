function spectra = beamforming(dataFiles, conditions, setup)
% PERFORM BEAMFORMING ANALYSIS ON MULTI-MICROPHONE PRESSURE READINGS
%
% INPUTS
% -------
% dataFiles                     string or cell  Path to h5 data file(s)
% setup.h                       Distance target from center microphone [m]
% setup.micPos                  positions + ID of microphones [m, -]
% setup.scanPlaneLimits         Scan plane Xmin, Xmax, Ymin, Ymax
% conditions.c                  Speed of sound [m/s]
% conditions.M_eff              Effective Mach number
%
% Optional
% setup.fRange                  2x1                 Frequency range
% setup.brokenMics              1D                  Broken mic indices
% setup.overlap                 double              Data overlap (0<->1)
% setup.scanPlaneResolution     double              Resolution [m]
% setup.diagonalRemoval         1/0                 Toggle
% setup.dataPortion             double              Data to use [0<->1]
% setup.flowVector              3x1                 Direction of flow      
% -------

%% Check inputs
tic;

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

% Frequency range
if isfield(setup, 'fRange')
    fRange = setup.fRange;
else
    fRange = [500 4000];
end

if numel(fRange) == 2
    fprintf('Frequency range: %i - %i Hz \n', fRange(1), fRange(2));
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
data = dataMic(:,1:nMics);

% Sampling frequency [Hz]
if isfield(setup, 'fs')
    fs = setup.fs;
    fprintf('Explicitly setting acquisition frequency: %i \n', fs);
else
    fs = h5read(dataFile1,'/Acquisition/SampleRate');
end

%% Environment preparation
fprintf('Extracting data...\n');

% Select portion of data to use
data = data(...
            max(1,round(size(data,1)/2 - (0.5*dataPortion)*size(data,1))) : ...
            round(size(data,1)/2       + (0.5*dataPortion)*size(data,1)),...
            :);

% Prepare mics
data(:,removeMics)   = [];
micPos(3,:)          = 0;
micPos(:,removeMics) = [];
nMics                = length(micPos);

% Prepare data chunks
chunkLength        = floor( timeChunk*fs );
chunkOverlapLength = floor( overlap*chunkLength );
chunkStartIdx      = 1 :  chunkLength-chunkOverlapLength : size(data,1)-chunkOverlapLength;
chunks             = [ chunkStartIdx(1:end-1); chunkStartIdx(1:end-1)+chunkLength-1 ];
   
%% Develop CSM
fprintf('Developing CSM...\n');

% Create data chunks
dataChunk = zeros(chunkLength, nMics, size(chunks,2));
for k = 1:size(chunks,2)
   
    chunkStart = chunks(1,k);
    chunkEnd   = chunks(2,k);
    
    dataSelect = data(chunkStart:chunkEnd, :);
    
    dataChunk(:,:,k) = dataSelect - repmat(mean(dataSelect, 1), chunkLength, 1);
    
end

% Frequencies for data chunks
dfChunk = fs/chunkLength;
fChunk  = (0:chunkLength/2-1)*dfChunk;

% Indices of frequency to be taken from fft of data chunks
fLowerIdx = floor(fRange(1).*chunkLength./fs + 1);
fUpperIdx = floor(fRange(2).*chunkLength./fs + 1);
f         = fChunk(fLowerIdx:fUpperIdx);

% Hanning window + overall sound pressure level (dB)
winHann       = repmat(hann(chunkLength), [1, size(data,2), size(chunks,2)]);
dataChunkHann = dataChunk .* winHann;
    
xChunkHann              = fft(dataChunkHann);
xChunkHann              = xChunkHann(1:length(fChunk),:,:)/(sqrt(sum(hann(chunkLength).^2)));
xChunkHann(2:end-1,:,:) = sqrt(2)*xChunkHann(2:end-1,:,:);
xChunkHannCollect       = xChunkHann;

xChunkff = permute(...
                sqrt(1/(chunkLength*dfChunk)) .* xChunkHannCollect(fLowerIdx:fUpperIdx,:,:), ...
                [2,3,1]);
       
% Average chunks
Cavg0 = zeros(nMics^2, fUpperIdx-fLowerIdx+1);
for j = 1:size(xChunkff,3)
    CC = rectpulse(conj(xChunkff(:,:,j)), size(xChunkff,1)) ...
                    .* repmat(xChunkff(:,:,j), [size(xChunkff,1),1]);
                
    Cavg0(:,j) = mean(real(CC) - 1i*imag(CC), 2);
end

Cavg = reshape(Cavg0, [nMics, nMics, fUpperIdx-fLowerIdx+1]);

% Diagonal removal
if diagonalRemoval == 1
    for i = 1:size(Cavg,3)
        Cavg(:,:,i) = Cavg(:,:,i) - diag(diag(Cavg(:,:,i)));
    end
end

Cavg = Cavg.*dfChunk;

%% Build scan plane

dx = scanPlaneLimits(1) : scanPlaneResolution : scanPlaneLimits(2);
dy = scanPlaneLimits(3) : scanPlaneResolution : scanPlaneLimits(4);

scanX = repmat(dx, length(dy), 1); 
scanY = repmat(flipud(dy'), 1, length(dx));
scanZ = h .* ones(size(scanX));

%% Beamforming
fprintf('Calculating radii...\n');

if flowCorrectedSteeringVector == 1
    
    % Steering vector for flow correction
    Mvec      = M_eff .* flowVector;
    beta2     = 1-M_eff^2;
    scanPlane = [scanX(:), scanY(:), -scanZ(:)];
    rvec      = repmat(micPos', size(scanPlane, 1), 1) - rectpulse(scanPlane, nMics);
    radii0    = sqrt(sum(rvec.^2, 2));
    
    Mvec2     = repmat(Mvec, size(rvec,1), 1);
    deltAllM  = (-dot(Mvec2, rvec, 2) + sqrt((dot(Mvec2, rvec, 2)).^2 + beta2*radii0.^2)) ...
                ./ (c*beta2);
                
    deltAll   = permute(reshape(deltAllM, [nMics,size(scanX)]), [2,3,1]);
    radii     = permute(reshape(radii0,   [nMics, size(scanX)]), [2,3,1]);

else
    
    % Basic steering vector, ignore the flow
    radii = sqrt((repmat(scanX,[1,1,nMics]) - repmat(reshape(micPos(1,:),[1,1,nMics]), [size(scanX),1])).^2 ...
               + (repmat(scanY,[1,1,nMics]) - repmat(reshape(micPos(2,:),[1,1,nMics]), [size(scanY),1])).^2 ...
               + (repmat(scanZ,[1,1,nMics]) - repmat(reshape(micPos(3,:),[1,1,nMics]), [size(scanX),1])).^2);

    deltAll = radii./c;
    
end

radius = repmat(radii, [1,1,1]);
delt   = repmat(deltAll, [1,1,1]);

fprintf('Beamforming...\n');
bar = waitbar(0, 'CFDBF: Performing beamforming ...');

B0   = zeros([numel(scanX), length(f)]);
B0g2 = zeros([numel(scanX), length(f)]); % intermediate map

for f_index = 1:length(f)
    waitbar(f_index/length(f));
    
    g     = ((-exp(-2*pi*1i*f(f_index).*delt)) ./ radius);
    gsel  = reshape(g, [numel(scanX), nMics]);
    C_ind = Cavg(:, :, f_index);
    
    for s = 1:size(gsel,1)
        Norm1 = norm(conj(gsel(s,:))*C_ind.'*gsel(s,:).');
        gnorm = norm(gsel(s,:));
        B0(s,f_index)  = Norm1/gnorm^4;
        B0g2(s,f_index) = Norm1/gnorm^2; % intermediate map
    end
end

% Final map
B          = reshape(B0,  [size(scanX), length(f)]);

% Store intermediate map for diagnostics
Bg2       = reshape(B0g2, [size(scanX), length(f)]);
B_sum      = sum(Bg2, 3);

pref       = 20e-6;
scanPlaneB = 20*log10(sqrt(B_sum)/(h*pref));
scanPlaneB = scanPlaneB - max(max(scanPlaneB));

% Finish
close(bar);

%% Structure output

spectra.B                  = B;
spectra.scanPlaneB         = scanPlaneB;
spectra.f                  = f;
spectra.scanPlaneX         = scanX;
spectra.scanPlaneY         = scanY;

%% Round up
fprintf('\nDone. It took %0.0f seconds \n\n', toc)
