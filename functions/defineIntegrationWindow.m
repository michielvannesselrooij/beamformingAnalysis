function setup = defineIntegrationWindow(setup)

%% Check inputs
if ~isfield(setup, 'intPlane_rel')
    error('Expecting a field setup.int_plane_rel');
end
     
%% Calculate integration window

% Shift window to relevant anchor point
if isfield(setup, 'anchor')
    
    setup.intPlane([1,3]) = setup.intPlane_rel([1,2]) + setup.anchor(1); % Shift x
    setup.intPlane([2,4]) = setup.intPlane_rel([3,4]) + setup.anchor(2); % Shift y
    
else
    
    warning('No anchor point for integration window defined. Centering at center mic');
    setup.intPlane = setup.intPlane_rel;
    
end

