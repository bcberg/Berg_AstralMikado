% Apply tension to an astral network (force sweep w/ manually set random seeds)
% Brady Berg, 2025-10-14
% (systematic rescaling based on fil_len, k_bend, viscosity)
set simul system
{
    time_step = 0.1
    viscosity = 0.5
    kT = [[0.0042/200]] % reduce default value by 200-fold (pN*um scale factor rel. to old params)
    verbose = 0
    % number of replicates (different seeds) specified in run_rep_full_shear_aster.sh
    random_seed = [[random.randint(0,10**9)]]
}

set space cell
{
    shape = rectangle
}

[[s = 1]]
new cell
{
    length = [[s]],[[s]] 
}

[[fiber_len = 0.1]]
set fiber actin
{
    rigidity       = 0.01
    segmentation   = [[fiber_len/5]]
    confine        = off %inside,200
    min_length     = [[fiber_len]]  
    activity       = none
}

% force will be applied to this fiber type
set fiber forced_fiber
{
    rigidity = 500
    segmentation = 0.05
    confine = off %inside,200
    min_length = [[s]]
    activity = none
    display = (line_width = 10; color=blue)
}

% fiber type to be anchored to the bottom of simulation space
set fiber base_fiber
{
    rigidity = 500
    segmentation = 0.05
    confine = off %inside,200
    min_length = [[s]]
    activity = none
    display = (line_width = 10; color=blue)
}

% define anchoring object
set hand strong_hand
{
    unbinding_rate = 0
    unbinding_force = inf
    display = ( size = 10; color=blue)
}

set single pivot
{
    hand = strong_hand
    stiffness = 500
    activity = fixed
}

set single temp_pivot
{
    hand = strong_hand
    stiffness = 500
    activity = fixed
}

% for visualizing center of asters
set solid core
{
    display = ( style=1; color=(1 0 0 0.33); size=1;) % (r g b a) where a \in [0,1] is transparency
}

set aster actinNode
{
    stiffness = 500, 250   % stiffness1 (pins filament ends at center), stiffness2 (provides torque)

}

% initialize asters
[[dens = 75]]
% target density N*fiber_len / A: [[dens]]
[[nFil = round(dens * s**2 / fiber_len)]]
% target number of filaments: [[nFil]]
[[nFilPerAster = list(range(1,25))]]
[[nAsters = round(nFil/nFilPerAster)]]
new [[nAsters]] actinNode
{
    type = astral % fibers are anchored at random positions near the center, pointing outward
    % avoids constraint from type 'radial', requiring sep of 25 nm by default
    solid = core
    radius = 0.03
    % point1 = center, 0.01
    fibers = [[nFilPerAster]], actin, ( length = [[fiber_len]];  end_state = static,static;)
}

% initialize fibers at top and bottom of simulation space
new 1 forced_fiber
{
    length = [[s]]
    position = 0 [[s/2 - 0.01]] 0
    orientation = 1 0 0     % "point" fiber to the right
    % Temporarily anchor forced_fiber while network crosslinks
    % syntax: `attach# = [object], [distance from reference to other end], [reference]
    attach1 = temp_pivot, 0, minus_end
    attach2 = temp_pivot, [[s/2]], minus_end
    attach3 = temp_pivot, 0, plus_end
}

new 1 base_fiber
{
    length = [[s]]
    position = 0 [[-s/2 + 0.01]] 0
    orientation = -1 0 0    % "point" fiber to the left
    % syntax: `attach# = [object], [distance from reference to other end], [reference]
    attach1 = pivot, 0, minus_end
    attach2 = pivot, [[s/10]], minus_end
    attach3 = pivot, [[2*s/10]], minus_end
    attach4 = pivot, [[3*s/10]], minus_end
    attach5 = pivot, [[4*s/10]], minus_end
    attach6 = pivot, [[5*s/10]], minus_end
    attach7 = pivot, [[6*s/10]], minus_end
    attach8 = pivot, [[7*s/10]], minus_end
    attach9 = pivot, [[8*s/10]], minus_end
    attach10 = pivot, [[9*s/10]], minus_end
    attach11 = pivot, 0, plus_end
}

% define crosslinkers
set hand binder
{
	binding_rate = 1
	binding_range = 0.001
	unbinding_rate = 0
    display = {width = 10; color = (1 1 0)}
}
set couple crosslinker
{
	hand1 = binder
	hand2 = binder
	stiffness = 50
	diffusion = 1
}
new [[30*nFil]] crosslinker {attach1 = actin}

% place crosslinkers near top and bottom to ensure
% they connect to the top and bottom filaments
new [[round(100*s)]] crosslinker {attach1 = forced_fiber}
% {
%     position = surface
%     placement = surface,, (y > [[s/2 - 0.01]])
% }

new [[round(100*s)]] crosslinker {attach1 = base_fiber}
% {
%     position = surface
%     placement = surface,, (y < [[-s/2 + 0.01]])
% }

% run system with no force; allows crosslinkers to bind
run system
{
    nb_steps = 500
    nb_frames = 5
}

[[F = [0.2 * x for x in range(6)]]]
% [[F]]
% activate applied force!
change fiber forced_fiber
{
    plus_end_force = 0 [[F/2]]
    minus_end_force = 0 [[F/2]]
}
% remove anchors on rforced_fiber
delete 3 temp_pivot

% report "initial" network position (after crosslinkers have bound)
report fiber:point initial_pos.txt
% report states of crosslinkers to determine network connectivity
report couple:link links.txt

run system
{
    nb_steps = 2500
    nb_frames = 25
}

% report final network position
report fiber:point final_pos.txt