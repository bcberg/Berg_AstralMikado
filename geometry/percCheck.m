function [percTF,percStats] = percCheck(crossings,nodes,D)
%PERCCHECK checks if an astral network has percolated
%Note: A network is percolated if a single connected component spans 
%       opposing sides of the square with corners (0,0) and (D,D) 
%       (i.e., top to bottom or left to right)
%   Inputs:
%       crossings (struct): has fields 'filCross' and 'centerCross'
%           'filCross' (numNodes x 2 double): list of pairs of filament 
%           indices corresponding to the filaments that cross at a
%           particular node; listed so that filCross(idx,1)<filCross(idx,2)
%           'centerCross': which filaments cross at each astral center;
%           each row represents a center
%       nodes (numNodes x 2 double): list of (x,y) coordinates of filament
%       crossings (EXCLUDING astral centers)
%       D (scalar): domain size
%   Outputs:
%       percTF (1x5 logical): results of various percolation tests:
%           (1): single component spans top to bottom
%           (2): single component spans left to right
%           (3): single component spans top to bottom OR left to right
%           (4): single component spans top to bottom AND left to right
%           (5): single component contains all filaments/asters
%       percStats (struct): info re: connected components in network,
%       has fields 'bins', 'binsizes'
%           'bins': bins(idx) gives index of connected component that
%           contains the filament specified by idx
%           'binsizes': binsizes(jdx) gives the number of filaments in
%           component jdx
[numAsters, astralNum] = size(crossings.centerCross);
numFil = numAsters * astralNum;

if astralNum == 1
    G = graph(crossings.filCross(:,1),crossings.filCross(:,2));
elseif astralNum >= 2
    edgesPerCenter = nchoosek(astralNum,2);
    centerEdges = zeros([numAsters * edgesPerCenter, 2]);
    for idx = 1:numAsters
        centerEdges((1+(idx-1)*edgesPerCenter):idx*edgesPerCenter,:) = ...
            nchoosek(crossings.centerCross(idx,:),2);
    end
    edges = [crossings.filCross; centerEdges];
    G = graph(edges(:,1),edges(:,2));
end
[bins,binsizes] = conncomp(G);

% identify nodes that fall outside of (0,0)--(D,0)--(D,D)--(0,D)--(0,0)
% each node is due to a filCross between two filaments idx and jdx
% --> these two filaments are in the same connected component compidx
%       (use bins(idx) to identify compidx)
% note: a node may satisfy two of the below conditions (e.g., above AND to
% the right)
nodesAbove = find(nodes(:,2) > D);
nodesBelow = find(nodes(:,2) < 0);
nodesLeft = find(nodes(:,1) < 0);
nodesRight = find(nodes(:,1) > D);

% top-to-bottom percolation check
if isempty(nodesAbove) || isempty(nodesBelow)
    percTB = false;
    tbSpanningComps = [];
else
    filsAbove = reshape(crossings.filCross(nodesAbove,:),1,[]);
    compsAbove = unique(bins(filsAbove));
    filsBelow = reshape(crossings.filCross(nodesBelow,:),1,[]);
    compsBelow = unique(bins(filsBelow));
    tbSpanningComps = intersect(compsAbove,compsBelow);
    percTB = any(tbSpanningComps);
end
% left-to-right percolation check
if isempty(nodesLeft) || isempty(nodesRight)
    percLR = false;
    lrSpanningComps = [];
else
    filsLeft = reshape(crossings.filCross(nodesLeft,:),1,[]);
    compsLeft = unique(bins(filsLeft));
    filsRight = reshape(crossings.filCross(nodesRight,:),1,[]);
    compsRight = unique(bins(filsRight));
    lrSpanningComps = intersect(compsLeft,compsRight);
    percLR = any(lrSpanningComps);
end
% all-sided percolation check
if isempty(tbSpanningComps) || isempty(lrSpanningComps)
    percTBLR = false;
else
    tblrSpanningComps = intersect(tbSpanningComps,lrSpanningComps);
    percTBLR = any(tblrSpanningComps);
end
% connected graph check
if max(binsizes) == numFil
    connTF = true;
else
    connTF = false;
end

percTF = [percTB, percLR, (percTB || percLR), percTBLR, connTF];
percStats.bins = bins;
percStats.binsizes = binsizes;