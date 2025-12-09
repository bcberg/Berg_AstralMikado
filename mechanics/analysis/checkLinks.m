function [percTF,percStats,filCross] = checkLinks(runDir,runIdx,cytoparams)
%CHECKLINKS analyzes the percolation state of a network from Cytosim
%   Inputs:
%       runDir (char vector or string): full filepath to directory which
%       groups all runs of parameter sweep together
%       runIdx (scalar integer): integer index of the simulation to
%       analyze (conventionally, the network at 0 applied force, though
%       links.txt is generated before applied force is switched on)
%       cytoparams (1x1 struct): parameters describing current network
%       parameters (needed to calculate filament quantities/indices); 
%       contains fields nFil, nFilPerAster, dataStartLine, ptsPerFiber
%   Outputs:
%       percTF (1 x 2 logical): true/false values from the following tests:
%           (1): spanning check; do the base and forced filament belong to
%           the same connected component?
%           (2): connectivity check; is there a single (i.e., unique)
%           connected component in the actin filament network? (This may or
%           may not include the base/forced filaments).
%       percStats (struct): info re: connected components in network,
%       has fields 'bins', 'binsizes'
%           'bins': bins(idx) gives index of connected component that
%           contains the filament specified by idx
%           'binsizes': binsizes(jdx) gives the number of filaments in
%           component jdx
%       filCross (n x 2 double): list of pairs of filament indices
%       corresponding to the filaments that are connected by a particular
%       crosslinker (may contain duplicates)
howManyAsters = round(cytoparams.nFil/cytoparams.nFilPerAster, ...
    TieBreaker="tozero");  
%   ^ attempt to match Python rounding behavior
howManyFibers = howManyAsters * cytoparams.nFilPerAster;
forcedFiberIdx = howManyFibers + 1;
baseFiberIdx = howManyFibers + 2;

% aster idx 1 groups filaments 1,2,...,nFilPerAster; and so on
centerCross = transpose(reshape(1:howManyFibers, ...
    [cytoparams.nFilPerAster,howManyAsters]));

% parse run%04d/links.txt
runCode = sprintf('run%04d',runIdx);
filename = fullfile(runDir,runCode,'links.txt');
% find which lines contain the crosslinker data
firstCrosslinkerLine = 5;   % structural property of report files
lines = readlines(filename);
lastCrosslinkerLine = length(lines) - 2;    % structural property
varNames = {'class','identity','fiber1','abscissa1','fiber2','abscissa2',...
    'force','cos_angle'};
opts = fixedWidthImportOptions('VariableNames',varNames, ...
    'VariableWidths',10*ones(1,length(varNames)),'DataLines', ...
    [firstCrosslinkerLine,lastCrosslinkerLine],'SelectedVariableNames', ...
    {'fiber1','fiber2'});
opts = setvartype(opts,{'fiber1','fiber2'},'double');
filCross = readmatrix(filename,opts);

if cytoparams.nFilPerAster == 1
    G = graph(filCross(:,1),filCross(:,2));
elseif cytoparams.nFilPerAster >= 2
    edgesPerCenter = nchoosek(cytoparams.nFilPerAster,2);
    centerEdges = zeros([howManyAsters * edgesPerCenter, 2]);
    for idx = 1:howManyAsters
        centerEdges((1+(idx-1)*edgesPerCenter):idx*edgesPerCenter,:) = ...
            nchoosek(centerCross(idx,:),2);
    end
    edges = [filCross; centerEdges];
    G = graph(edges(:,1),edges(:,2));
end
if ismultigraph(G)
    % ignore when multiple crosslinkers connect the same two filaments
    G = simplify(G);
end

[bins,binsizes] = conncomp(G);

% spanning check
if length(bins) < baseFiberIdx
    % i.e., the graph does not contain nodes corresponding to both the 
    % base and forced filament
    spanTF = false;
elseif bins(forcedFiberIdx) == bins(baseFiberIdx)
    spanTF = true;
else
    spanTF = false;
end
% connected graph check (single connected component of actin)
actinComponents = unique(bins(1:howManyFibers));
if isscalar(actinComponents)
    connTF = true;
else
    connTF = false;
end

percTF = [spanTF, connTF];
percStats.bins = bins;
percStats.binsizes = binsizes;
end