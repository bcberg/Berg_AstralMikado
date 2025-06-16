function curve = getPercCurve(l,D,astralNum,densSpec,Nsamp,Ncores, ...
    saveDirectory)
%GETPERCCURVE Computes percolation probability curves for a single astral
% network type
% Percolation definitions are ordered as follows:
%   (1): single component spans top to bottom
%   (2): single component spans left to right
%   (3): single component spans top to bottom OR left to right
%   (4): single component spans top to bottom AND left to right
%   (5): single component contains all filaments/asters
%   Inputs:
%       l (scalar): length of individual filament
%       D (scalar): domain size; asters are distributed in the square
%       with corners at (0,0) and (D,D)
%       astralNum (scalar): (whole) number of filaments per aster
%       densSpec (1x3 double): specifies density sweep range and number of
%       log-spaced points; range is 10.^densSpec(1:2), number of points is
%       densSpec(3)
%       Nsamp (scalar): number of networks to sample per density value
%       Ncores (scalar): number of CPU cores allotted to job
%       saveDirectory (char vector or string): name of folder in which to 
%       save simulation data; pass empty value to skip writing .mat file
%   Outputs:
%       curve (6 x numUniqueDensVals double): curve(1,:) gives line density
%       values, curve(2:6,:) gives percolation probability estimates in the
%       order specified above
% Example call from command line
% matlab -nodesktop -r "getPercCurve(1,20,5,[0,log10(50),50],2000,24,'$TMPDIR')"

filename = sprintf("percProbs_l%02i_D%02i_an%02i",l,D,astralNum);

%%%%%%%%%%%%%% Sampling notes %%%%%%%%%%%%%%
% (a) At smaller l/D, probabilities are 0 until higher densities.
% Manually adjust densSpec to capture the transition
% (b): There may be fewer density values in certain percolation curves due
% to requiring integer numbers of asters (density resolution is coarser at
% higher astral number, larger l/D)
% (c): Once two consecutive densities satisfy ALL percProbs = 1, all
% sampling at higher density values is skipped and remaining probabilities 
% are set to 1

%%%%%%%%%%%%%% Percolation probability estimation %%%%%%%%%%%%%%
densityRange = logspace(densSpec(1),densSpec(2),densSpec(3));
numFilRange = densityRange * D^2 / l;
numAstersRange = round(numFilRange / astralNum);
numAstersUsed = unique(numAstersRange);
numAstersUsed = numAstersUsed(numAstersUsed > 0);
numUniqueDensVals = length(numAstersUsed);
actualDensities = numAstersUsed * (astralNum * l / D^2);
percProbs = zeros(5,numUniqueDensVals);

pool = parpool(Ncores);
useMEX = false;      % adjust depending on machine
addAttachedFiles(pool,{'generateAstralNetwork.m','percCheck.m', ...
    'percOneSample.m'})
F = cell(1,numUniqueDensVals);
for idx = 1:numUniqueDensVals
    F{idx} = parallel.FevalFuture.empty(Nsamp,0);
    for jdx = 1:Nsamp
        F{idx}(jdx) = parfeval(pool,@percOneSample,1, ...
            numAstersUsed(idx),l,D,astralNum,useMEX);
    end
end

samples = fetchOutputs(F{1});
probabilities = sum(samples,1) / Nsamp;
percProbs(1:5,1) = transpose(probabilities);
fprintf('Completed %02i/%02i densities\n',1,numUniqueDensVals)
for idx = 2:numUniqueDensVals
    samples = fetchOutputs(F{idx});
    probabilities = sum(samples,1) / Nsamp;
    percProbs(1:5,idx) = transpose(probabilities);
    fprintf('Completed %02i/%02i densities\n',idx,numUniqueDensVals)
    % See sampling note (c) for lines 74-82
    if idx<numUniqueDensVals && prod(percProbs(1:5,(idx-1):idx),"all") == 1
        fprintf('Skipping remaining %i densities\n',numUniqueDensVals-idx);
        percProbs(1:5,(idx+1):numUniqueDensVals) = 1;
        for killIdx = (idx+1):numUniqueDensVals
            cancel(F{killIdx})
        end
        break
    end
end
clear('F')
delete(pool)

curve = [actualDensities; percProbs];
if ~isempty(saveDirectory)
    fileToSave = fullfile(saveDirectory,filename + ".mat");
    save(fileToSave)    % saves all variables, including function call args
end
end