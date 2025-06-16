function [CoM_vs_time] = parseTimeSeries(dir,runVal,frameVals,cytoparams)
%PARSETIMESERIES imports a framewise-report of network position and
%computes center-of-mass for each frame
%   Inputs:
%       dir (character vector or string): full filepath to directory which
%       contains all runs of timeseries_check sweep
%       runVal (1x1 double): integer value of the run to analyze
%       frameVals (n x 1 double): indices of frames to analyze 
%       cytoparams (1x1 struct): parameters describing current network
%       -> needed as inputs to fiberPointLines function; contains
%       fields nFil, nFilPerAster, dataStartLine, ptsPerFiber
%   Outputs:
%       CoM_vs_time (n x 2 double): list of (x,y) coordinates of center of
%       mass corresponding to each specified frame
varNames = {'identity', 'posX', 'posY', 'curvature'};
dataLines = fiberPointLines(cytoparams.nFil,cytoparams.nFilPerAster,...
    cytoparams.dataStartLine,cytoparams.ptsPerFiber);
opts = fixedWidthImportOptions('VariableNames',varNames,...
    'VariableWidths',10*ones(1,length(varNames)),'DataLines',dataLines,...
    'SelectedVariableNames',{'posX','posY'});
opts = setvartype(opts,{'posX','posY'},'double');

numFrames = length(frameVals);
CoM_vs_time = zeros(numFrames,2);
runDir = sprintf(strcat(dir, "run%04i/"), runVal);
reportPattern = strcat(runDir, "pos%04i.txt");
for idx = 1:numFrames
    thisFrame = frameVals(idx);
    thisReport = sprintf(reportPattern,thisFrame);
    points = readmatrix(thisReport,opts);
    CoM_vs_time(idx,:) = mean(points,1); % computes means of columns
end
end