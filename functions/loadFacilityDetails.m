function setup = loadFacilityDetails(setup)
% Augments setup structure with facility-specific settings

% LTT
if contains(setup.facility, 'LTT')
    
    % Microphone positions
    load('mic_positions_LTT.mat', 'micPos')    
    setup.micPos = micPos;
    
    % Settings
    setup.flowVector      = [-1, 0, 0];                 % Flow vector
    setup.scanPlaneLimits = [-1, 1, -1.25/2, 1.25/2];   % Scan plane
    setup.anchor          = [setup.wing.TEpos, 0];      % Anchor point integration window
    setup.wing.rotAxis    = 'vertical';                 % Axis of rotation

    
% A-tunnel
elseif contains(setup.facility, '(A)') || contains(setup.facility, 'A-tunnel')
    
    % Microphone positions
    load('mic_positions_A_vertical.mat', 'micPos')    
    setup.micPos = micPos;
    
    % Settings
    setup.flowVector      = [0, 1, 0];                  % Flow vector
    setup.scanPlaneLimits = [-0.2, 0.2, -1, 1];         % Scan plane
    setup.anchor          = [0, setup.wing.TEpos];      % Anchor point integration window
    setup.wing.rotAxis    = 'horizontal';               % Axis of rotation
    
    
% Unfinished facilities
elseif strcmp(setup.facility, 'External')
    error('Programming to be completed for this facility');
    
elseif contains(setup.facility, 'ROSI')
    error('Programming to be completed for this facility');
  
    
% Unsupported
else
    error('No suitable facility specified in configuration');
    
end
