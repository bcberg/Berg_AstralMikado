function [connTF,connStats] = connCheck(crossings)
%CONNCHECK analyzes the connectedness of an astral network
%   Inputs:
%       crossings (struct): has fields 'filCross' and 'centerCross'
%           'filCross' (numNodes x 2 double): list of pairs of filament 
%           indices corresponding to the filaments that cross at a
%           particular node; listed so that filCross(idx,1)<filCross(idx,2)
%           'centerCross': which filaments cross at each astral center;
%           each row represents a center
%   Outputs:
%       connTF (logical): true if graph of filament intersections has one
%       connected component, false otherwise
%       connStats (struct): statistics describing network connectedness,
%       has fields 'bins', 'binsizes'
%           'bins': bins(idx) gives index of connected component that
%           contains the filament specified by idx
%           'binsizes': binsizes(jdx) gives the number of filaments in
%           component jdx
numAsters = size(crossings.centerCross,1);
astralNum = size(crossings.centerCross,2);
numFil = numAsters * astralNum;

if astralNum == 1
    G = graph(crossings.filCross(:,1),crossings.filCross(:,2));
elseif astralNum >= 2
    edgesPerCenter = nchoosek(astralNum,2);
    centerEdges = zeros([numAsters * nchoosek(astralNum,2), 2]);
    for idx = 1:numAsters
        centerEdges((1+(idx-1)*edgesPerCenter):idx*edgesPerCenter,:) = ...
            nchoosek(crossings.centerCross(idx,:),2);
    end
    edges = [crossings.filCross; centerEdges];
    G = graph(edges(:,1),edges(:,2));
end

[bins,binsizes] = conncomp(G);
if max(binsizes) == numFil
    connTF = true;
else
    connTF = false;
end
connStats.bins = bins;
connStats.binsizes = binsizes;