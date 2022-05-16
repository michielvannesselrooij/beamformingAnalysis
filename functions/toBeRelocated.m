h        = setup.d_hor + setup.x_rot*sind(conditions.AoA);                 % Distance from array to scan plane [m]


fs       = 51200;    



flowWidth = 0.35; % [m]
flowWidth = 0.91; % [m]
shear_layer_width = flowWidth + setup.x_rot*sind(conditions.AoA); 
M     = conditions.U/c;                     % Mach number
M_eff = M*shear_layer_width/h;

% Correct integration boundaries to absolute coordinates
setup.int_plane(1)              = 0.5*(setup.af_loc(1) + setup.af_loc(3))        +setup.int_plane_rel(1)*setup.c;
setup.int_plane(2)              = setup.af_loc(4)      + setup.int_plane_rel(2)*setup.c;
setup.int_plane(3)              = 0.5*(setup.af_loc(1) + setup.af_loc(3))        +setup.int_plane_rel(3)*setup.c;
setup.int_plane(4)              = setup.af_loc(4)      + setup.int_plane_rel(4)*setup.c;
x_TE                            = 0.5*(setup.af_loc(1) + setup.af_loc(3));
y_TE                            = setup.af_loc(4);

scan_plane_limits = [0, 1, -0.5, 0.5];                  
scan_plane_limits               = [-.5, 0.7, -1.25/2, 1.25/2];
