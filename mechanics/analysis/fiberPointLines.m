function dataLines = fiberPointLines(nFil,nFilPerAster,dataStartLine,ptsPerFiber)
%FIBERPOINTLINES returns the line intervals where points on network
%filaments are reported by Cytosim for particular network parameters
% Note: parsing scheme assumes all aster-filaments are listed together
%   Inputs:
%       nFil (1x1 double): total number of filaments (raw parameter used in
%       .cym.tpl file)
%       nFilPerAster (1x1 double): number of filaments per aster
%       dataStartLine (1x1 double): first line of data to read from report
%       ptsPerFiber (1x1 double): number of discrete points used by Cytosim
%       to represent fiber
%   Outputs:
%       dataLines (nFil x 2 double): list of intervals of lines to read
%       from text file (lines containing coordinates of fiber points)
howManyAsters = round(nFil/nFilPerAster,TieBreaker="tozero");  % attempt to match Python rounding behavior
howManyFibers = howManyAsters * nFilPerAster;
lastDataLine = dataStartLine + (ptsPerFiber+1) * howManyFibers - 2;
startLines = (dataStartLine:(ptsPerFiber+1):lastDataLine)';
dataLines = [startLines, startLines + ptsPerFiber - 1];
end