% Shear an "original" Mikado network
% --> for plotting timeseries of deformation (steady state validation)
% Brady Berg, 2024-05-15

set simul system
{
    time_step = 0.01
    viscosity = 0.01
    verbose = 0
}

set space cell
{
    shape = rectangle
}

[[D = [5,10,15,20]]]
new cell
{
    length = [[D]],[[D]] 
}

[[fiber_len = 1]]
set fiber actin
{
    rigidity       = 20
    segmentation   = [[fiber_len / 5]]
    confine        = off %inside,200
    min_length     = 0.5  
    activity       = none
}

% force will be applied to this fiber type
set fiber rforced_fiber
{
    rigidity = 1e3
    segmentation = 0.5
    confine = off %inside,200
    min_length = [[D]]
    activity = none
}

% fiber type to be anchored to the bottom of simulation space
set fiber base_fiber
{
    rigidity = 1e3
    segmentation = 0.5
    confine = off %inside,200
    min_length = [[D]]
    activity = none
}

% define anchoring object
set hand strong_hand
{
    unbinding_rate = 0
    unbinding_force = 1e6
    display = ( size = 10; color=blue)
}

set single pivot
{
    hand = strong_hand
    stiffness = 1000
    activity = fixed
}

% for visualizing center of asters
set solid core
{
    display = ( style=4 ; color=red; )
}

set aster actinNode
{
    stiffness = 1000, 500   % I'm not sure what these do, but they match the default in `aster.cym`
}

% initialize asters
[[dens = 7.5]]
% probably exceed threshold: dens = 10 for D=5, dens = 7.5 for D = 10, 15
% target density N*fiber_len / A: [[dens]]
[[nFil = round(dens * D**2 / fiber_len)]]
% target number of filaments: [[nFil]]
[[nFilPerAster = 1]]
[[nAsters = round(nFil/nFilPerAster)]]
new [[nAsters]] actinNode
{
    solid = core
    radius = 0.3
    point1 = center, 0.3
    fibers = [[nFilPerAster]], actin, ( length = [[fiber_len]];  end_state = green, green;)
}

% initialize fibers at top and bottom of simulation space
new 1 rforced_fiber
{
    length = [[D]]
    position = 0 [[D/2 - 0.01]] 0
    orientation = 1 0 0     % "point" fiber to the right
}

new 1 base_fiber
{
    length = [[D]]
    position = 0 [[-D/2 + 0.01]] 0
    orientation = -1 0 0    % "point" fiber to the left
    % syntax below is `attach# = [object] [distance from reference to other end] [reference]
    attach1 = pivot, 0, minus_end
    attach2 = pivot, [[D/10]], minus_end
    attach3 = pivot, [[2*D/10]], minus_end
    attach4 = pivot, [[3*D/10]], minus_end
    attach5 = pivot, [[4*D/10]], minus_end
    attach6 = pivot, [[5*D/10]], minus_end
    attach7 = pivot, [[6*D/10]], minus_end
    attach8 = pivot, [[7*D/10]], minus_end
    attach9 = pivot, [[8*D/10]], minus_end
    attach10 = pivot, [[9*D/10]], minus_end
    attach11 = pivot, 0, plus_end
}

% define crosslinkers
set hand binder
{
	binding_rate = 10
	binding_range = 0.01
	unbinding_rate = 0
    display = (size = 6; color = yellow)
}
set couple crosslinker
{
	hand1 = binder
	hand2 = binder
	stiffness = 100
	diffusion = 10
}
new [[30*nFil]] crosslinker

% place crosslinkers near top and bottom to ensure
% they connect to the top and bottom filaments
new [[10*D]] crosslinker 
{
    position = surface
    placement = surface,, (y > [[D/2 - 0.01]])
}

new [[10*D]] crosslinker 
{
    position = surface
    placement = surface,, (y < [[-D/2 + 0.01]])
}

% run system with no force; allows crosslinkers to bind
run system
{
    nb_steps = 500
    nb_frames = 5
}

[[F = 20]]
% [[F]]
% activate applied force!
change fiber rforced_fiber
{
    plus_end_force = [[F]] 0        % rightward force applied to plus end of fiber
}

run system
{
    nb_steps = 10000
    nb_frames = 100
}

% report, framewise, network position over entire simulation -> needs to be done outside of .cym
