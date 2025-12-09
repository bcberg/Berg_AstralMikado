function p = estPercProb_eqSpaced(numAsters,l,D,astralNum,Nsamp,useMEX)
%ESTPERCPROB_EQSPACED uses Monte Carlo sampling to estimate percolation
%probabilities of astral networks; filaments are equally spaced about centers
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
%       p (1x5 double): estimated percolation probabilities for astral
%       networks with the given parameters; ordered as probability that
%           (1): single component spans top to bottom
%           (2): single component spans left to right
%           (3): single component spans top to bottom OR left to right
%           (4): single component spans top to bottom AND left to right
%           (5): single component contains all filaments/asters
percCount = zeros(1,5);
numFil = numAsters * astralNum;     % mex gets slow if numFil too high
if useMEX && numFil < 1e4
    for idx = 1:Nsamp
        [network,crossings,~] = generateAstralNetwork_eqSpaced_mex(numAsters,l,D,...
            astralNum, true);
        [percTF,~] = percCheck(crossings,network.nodes,D);
        percCount = percCount + percTF;
    end
else
    for idx = 1:Nsamp
        [network,crossings,~] = generateAstralNetwork_eqSpaced(numAsters,l,D, ...
            astralNum, true);
        [percTF,~] = percCheck(crossings,network.nodes,D);
        percCount = percCount + percTF;
    end
end
p = percCount / Nsamp;
end