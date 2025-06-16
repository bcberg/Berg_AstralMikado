function [initCoM, finalCoM] = parseForceSweep(dir,runVals,cytoparams)
%PARSEFORCESWEEP imports data from a force sweep over a network and
%computes initial and final center of mass coordinates at each force
%   Inputs:
%       dir (character vector or string): full filepath to directory which
%       groups all runs of 2D-sweep together
%       runVals (n x 1 double): integer values of the runs which correspond
%       to a force sweep over one network type
%       cytoparams (1x1 struct): parameters describing current network
%       parameters needed as inputs to fiberPointLines function; contains
%       fields nFil, nFilPerAster, dataStartLine, ptsPerFiber
%   Outputs:
%       initCoM (n x 2 double): list of center-of-mass coordinates at
%       initial state for each applied force (n = length(forceVals))
%       finalCoM (n x 2 double): list of center-of-mass coordinates at
%       steady-state for each applied force (n = length(forceVals))
varNames = {'identity', 'posX', 'posY', 'curvature'};
dataLines = fiberPointLines(cytoparams.nFil,cytoparams.nFilPerAster,...
    cytoparams.dataStartLine,cytoparams.ptsPerFiber);
opts = fixedWidthImportOptions('VariableNames',varNames,...
    'VariableWidths',10*ones(1,length(varNames)),'DataLines',dataLines,...
    'SelectedVariableNames',{'posX','posY'});
opts = setvartype(opts,{'posX','posY'},'double');

% store (x,y) coords of initial and final center of mass for each force
initCoM = zeros(length(runVals),2);
initFilePattern = strcat(dir, "run%04i", "/initial_pos.txt");
finalCoM = zeros(length(runVals),2);
finalFilePattern = strcat(dir, "run%04i", "/final_pos.txt");
for idx = 1:length(runVals)
    thisRun = runVals(idx);
    initFilename = sprintf(initFilePattern,thisRun);
    initPts = readmatrix(initFilename,opts);
    initCoM(idx,:) = mean(initPts,1); % row vector containing mean of each column
    finalFilename = sprintf(finalFilePattern,thisRun);
    finalPts = readmatrix(finalFilename,opts);
    finalCoM(idx,:) = mean(finalPts,1);
end
end