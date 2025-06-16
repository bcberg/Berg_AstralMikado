function samp = percOneSample(numAsters,l,D,astralNum,useMEX)
%PERCONESAMPLE Generates one astral network and checks if it is percolated
%   Inputs:
%       numAsters (scalar): (whole) number of asters to distribute
%       l (scalar): length of individual filament
%       D (scalar): domain size; asters are distributed in the square
%       with corners at (0,0) and (D,D)
%       astralNum (scalar): (whole) number of filaments per aster
%       Nsamp (scalar): integer number of networks to sample
%       useMEX (boolean): passing true allows use of MEX function for
%       network generation, passing false uses pure Matlab function
%   Outputs:
%       samp (1x5 logical): results of various percolation tests:
%           (1): single component spans top to bottom
%           (2): single component spans left to right
%           (3): single component spans top to bottom OR left to right
%           (4): single component spans top to bottom AND left to right
%           (5): single component contains all filaments/asters
numFil = numAsters * astralNum;     % mex gets slow if numFil too high
if useMEX && numFil < 1e4
    [network,crossings,~] = generateAstralNetwork_mex(numAsters,l,D,...
        astralNum, true);
    [samp,~] = percCheck(crossings,network.nodes,D);
else
    [network,crossings,~] = generateAstralNetwork(numAsters,l,D, ...
        astralNum, true);
    [samp,~] = percCheck(crossings,network.nodes,D);
end
end