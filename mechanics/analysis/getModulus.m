function [modulus,linefit,gof] = getModulus(initCoM,finalCoM,forceVals,type)
%GETMODULUS fits a line to force vs. displacement data
%   Inputs:
%       initCoM (n x 2 double): list of center-of-mass coordinates at
%       initial state for each applied force (n = length(forceVals))
%       finalCoM (n x 2 double): list of center-of-mass coordinates at
%       steady-state for each applied force (n = length(forceVals))
%       forceVals (n x 1 double): list of applied shear forces
%       corresponding to each row in finalCoM; must start at 0
%       type (character vector or string): either "shear" or "extension",
%       indicates stress/strain direction
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
% constrain intercept coeff. p2 to be 0
fitopts = fitoptions('poly1','Lower',[-Inf,0],'Upper',[Inf,0]);
displace = finalCoM(:,coordIdx) - initCoM(:,coordIdx);
[linefit,gof] = fit(displace,forceVals,'poly1',fitopts);
modulus = linefit.p1;
end