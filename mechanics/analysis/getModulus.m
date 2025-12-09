function [modulus,linefit,gof] = getModulus(initCoM,finalCoM,forceVals, ...
    type,filtering)
%GETMODULUS estimates a 2D elastic modulus from force-displacement data
%   Inputs:
%       initCoM (n x 2 double): list of center-of-mass coordinates at
%       initial state for each applied force (n = length(forceVals))
%       finalCoM (n x 2 double): list of center-of-mass coordinates at
%       steady-state for each applied force (n = length(forceVals))
%       forceVals (n x 1 double): list of applied shear forces
%       corresponding to each row in finalCoM; must start at 0
%       type (character vector or string): either "shear" or "extension",
%       indicates stress/strain direction
%       filtering (1x1 struct): specifies if (and how) to filter networks,
%       contains 3 true/false fields: 
%           -- 'TF' determines if any filtering is done
%           -- 'spanCheck' reports if a spanning component was detected
%       in the network at the end of the initial crosslinking time window
%           -- 'catchDrops' assigns 0 modulus to networks where
%       displacement decreases from forceVals(1)=0 to forceVals(2)>0
%       (*if unspecified, this field is by default TRUE)
%   Outputs:
%       modulus (1x1 double): coefficient of linear term from linefit
%       linefit (cfit object): result of linear fit to force vs.
%       displacement data
%       gof (gof structure): goodness-of-fit statistics
if strcmp(type,"shear")
    coordIdx = 1;
elseif strcmp(type,"extension")
    coordIdx = 2;
else
    error("type argument must be one of 'shear' or 'extension'")
end
%%%%% OPTION 1 %%%%%
%%% unconstrained coefficient values, default fit options
% fitopts = fitoptions('poly1');

%%%%% OPTION 2 %%%%%
%%% constrain intercept coeff. p2 to be 0
% fitopts = fitoptions('poly1','Lower',[-Inf,0],'Upper',[Inf,0]);

%%%%% OPTION 3 %%%%% (current selection)
%%% use bisquare robust fitting with intercept fixed at 0
fitopts = fitoptions('poly1','Robust','Bisquare','Lower',[-Inf,0], ...
    'Upper',[Inf,0]);

displace = finalCoM(:,coordIdx) - initCoM(:,coordIdx);
% fit with forceVals as the independent variable for better performance
[linefit,gof] = fit(forceVals,displace,'poly1',fitopts);
if filtering.TF
    % try to catch networks that may not contain a spanning component, or
    % may otherwise result in poor fitting behavior, as identified by 
    % at least ONE of the following:
    % (1) lack of spanning component as determined by checkLinks function
    %   ^ Do this check in top-level analysis script & pass result in here
    % (2) adjacent displacements, except between forceVals(1)=0 and 
    % forceVals(2)>0, that differ by at least 75% of max. displacement
    %   ^ Attempt to catch networks with oscillatory displacement data
    % (3) displacement decreases from forceVals(1)=0 to forceVals(2)>0
    %       (condition (3) can be toggled and is by default TRUE)
    
    dispDiffs = diff(displace);
    dispFlag = abs(dispDiffs(2:end)) > (0.75 * max(abs(displace)));
    if ~filtering.spanCheck % (1) if a spanning component was NOT found
        modulus = 0;
    elseif any(dispFlag)    % (2) large oscillations in displacements
        modulus = 0;
    elseif dispDiffs(1)<0    % (3) initial decrease in displacements
        if ~isfield(filtering,'catchDrops')
            filtering.catchDrops = true;
        end
        if filtering.catchDrops
            modulus = 0;
        else
            modulus = 1/linefit.p1;
        end
    else
        modulus = 1/linefit.p1;
    end
else
    % report modulus from "blind" linear regression
    modulus = 1/linefit.p1;
end
end