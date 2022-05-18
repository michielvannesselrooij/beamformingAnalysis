function setup = defineIntegrationWindow(setup, conditions)

%% Check inputs
if ~isfield(setup, 'int_plane_rel')
    error('Expecting a field setup.int_plane_rel');
end

if ~isfield(setup, 'c')
    error('Expecting a field setup.c');
end

if ~isfield(setup, 'x_rot')
    error('Expecting a field setup.x_rot');
end

if ~isfield(conditions, 'AoA')
    error('Expecting a field conditions.AoA');
end

%% Location of the airfoil corners [m] [Xlowerleft Ylowerleft Xupperright Yupperright]

if strcmp(setup.tunnel, 'A')
    af_loc = [0.3,...
             -setup.c + (setup.c - setup.x_rot)*(1-cosd(conditions.AoA)) + setup.d_ver,...
              0.7,...
             -setup.x_rot*(1-cosd(conditions.AoA)) + setup.d_ver]; 
         
elseif strcmp(setup.tunnel,'LTT')
    af_loc = [setup.x_rot*(1-cosd(conditions.AoA)) + setup.d_ver,...
              -1.25/2,...
               setup.c - (setup.c - setup.x_rot)*(1-cosd(conditions.AoA)) + setup.d_ver,...
               1.25/2];
        
else
    error('No valid tunnel specified (A or LTT)');
    
end

setup.af_loc = af_loc;
     
%% Calculate integration window
setup.int_plane = [];

setup.int_plane(1)              = 0.5*(af_loc(1) + af_loc(3))        +setup.int_plane_rel(1)*setup.c;
setup.int_plane(2)              = af_loc(4)      + setup.int_plane_rel(2)*setup.c;
setup.int_plane(3)              = 0.5*(af_loc(1) + af_loc(3))        +setup.int_plane_rel(3)*setup.c;
setup.int_plane(4)              = af_loc(4)      + setup.int_plane_rel(4)*setup.c;

