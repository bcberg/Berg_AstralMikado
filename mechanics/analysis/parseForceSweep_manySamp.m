function [initCoM,finalCoM,allInitCoM,allFinalCoM] = parseForceSweep_manySamp(dir, ...
    runVals,initFrames,finalFrames,cytoparams)
%PARSEFORCESWEEP_MANYSAMP imports data from a force sweep over a network
%and computes initial and final center of mass coordinates at each force;
%determines average initial and final positions from multiple position
%samples
%   Inputs:
%       dir (character vector or string): full filepath to directory which
%       groups all runs of 2D-sweep together
%       runVals (n x 1 double): integer values of the runs which correspond
%       to a force sweep over one network type
%       initFrames (m x 1 double): integer values of the report frames from
%       the "initial position" window (established w/in Cytosim)
%       finalFrames (l x 1 double): integer values of the report frames from
%       the "final position" window (established w/in Cytosim)
%       cytoparams (1x1 struct): parameters describing current network
%       parameters needed as inputs to fiberPointLines function; contains
%       fields nFil, nFilPerAster, dataStartLine, ptsPerFiber
%   Outputs:
%       initCoM (n x 2 double): list of center-of-mass coordinates at
%       initial state for each applied force (n = length(forceVals))
%       finalCoM (n x 2 double): list of center-of-mass coordinates at
%       steady-state for each applied force (n = length(forceVals))
%       allInitCoM (n x 2 x m double): raw center-of-mass coordinates at
%       all sampled frames from "initial position" window;
%       allInitCoM(:,:,k) gives center-of-mass at frame initFrames(k)
%       allFinalCoM (n x 2 x l double): raw center-of-mass coordinates at
%       all sampled frames from "final position" window;
%       allFinalCoM(:,:,k) gives center-of-mass at frame finalFrames(k)
varNames = {'identity', 'posX', 'posY', 'curvature'};
dataLines = fiberPointLines(cytoparams.nFil,cytoparams.nFilPerAster,...
    cytoparams.dataStartLine,cytoparams.ptsPerFiber);
opts = fixedWidthImportOptions('VariableNames',varNames,...
    'VariableWidths',10*ones(1,length(varNames)),'DataLines',dataLines,...
    'SelectedVariableNames',{'posX','posY'});
opts = setvartype(opts,{'posX','posY'},'double');

reportFilePattern = strcat(dir,"run%04d/","pos%04d.txt");
n = length(runVals);
m = length(initFrames);
l = length(finalFrames);
allInitCoM = zeros(n,2,m);
allFinalCoM = zeros(n,2,l);
for idx = 1:n
    thisRun = runVals(idx);
    for kdx = 1:m
        % parse initial position data
        initFilename = sprintf(reportFilePattern,thisRun,initFrames(kdx));
        initPts = readmatrix(initFilename,opts);
        allInitCoM(idx,:,kdx) = mean(initPts,1); % take mean along cols
    end
    
    for kdx = 1:l
        % parse final position data
        finalFilename = sprintf(reportFilePattern,thisRun,finalFrames(kdx));
        finalPts = readmatrix(finalFilename,opts);
        allFinalCoM(idx,:,kdx) = mean(finalPts,1);
    end
end
initCoM = mean(allInitCoM,3);
finalCoM = mean(allFinalCoM,3);
end