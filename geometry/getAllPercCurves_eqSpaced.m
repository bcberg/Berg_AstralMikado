function curves = getAllPercCurves_eqSpaced(l,D,astralNumList,densSpec, ...
    Nsamp,Ncores,saveDirectory)
%GETALLPERCCURVES_EQSPACED Computes percolation probability curves for a 
% list of astral network types; filaments are equally spaced about centers
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
%       astralNumList (double vector): (whole) number of filaments per
%       aster, function computes a curve for each value in this list
%       densSpec (1x3 double): specifies density sweep range and number of
%       log-spaced points; range is 10.^densSpec(1:2), number of points is
%       densSpec(3)
%       Nsamp (scalar): number of networks to sample per density value
%       Ncores (scalar): number of CPU cores allotted to job
%       saveDirectory (char vector or string): name of folder in which to 
%       save simulation data; pass empty value to skip writing .mat file
%   Outputs:
%       curves (struct): curves.an%02i(1,:) gives line density values,
%       curves.an%02i(2:6,:) gives percolation probability estimates in the
%       order specified above (replace %02i with desired astral number)
filename = sprintf("percProbs_eqSpaced_l%02i_D%02i",l,D);

%%%%%%%%%%%%%% Sampling notes %%%%%%%%%%%%%%
% (1) At smaller l/D, probabilities are 0 until higher densities.
% Manually adjust densSpec to capture the transition
% (2): There may be fewer density values in certain percolation curves due
% to requiring integer numbers of asters (density resolution is coarser at
% higher astral number, larger l/D)
% (c): Once two consecutive densities (for a given astralNum) satisfy ALL
% percProbs = 1, all sampling at higher density values is skipped and 
% remaining probabilities are set to 1

%%%%%%%%%%%%%% Percolation probability estimation %%%%%%%%%%%%%%
densityRange = logspace(densSpec(1),densSpec(2),densSpec(3));
numFilRange = densityRange * D^2 / l;
numNetTypes = length(astralNumList);
astralNumList = reshape(astralNumList,[numNetTypes,1]);
numAstersArray = round(repmat(numFilRange,[numNetTypes,1]) ./ ...
    repmat(astralNumList,[1,densSpec(3)]));
numAstersUsed = cell(numNetTypes,1);
numUniqueDensVals_byrow = zeros(numNetTypes,1);
actualDensities = cell(numNetTypes,1);
percProbs = cell(numNetTypes,1);
for idx = 1:numNetTypes
    numAstersUsed{idx} = unique(numAstersArray(idx,:));
    numAstersUsed{idx} = numAstersUsed{idx}(numAstersUsed{idx}>0);
    numUniqueDensVals_byrow(idx) = numel(numAstersUsed{idx});
    actualDensities{idx} = astralNumList(idx)*numAstersUsed{idx}*l/D^2;
    percProbs{idx} = zeros(5,numUniqueDensVals_byrow(idx));
end

Njobs = sum(numUniqueDensVals_byrow);
storeHere = zeros(Njobs,2);
pool = parpool(Ncores);
useMEX = false;      % adjust depending on machine
addAttachedFiles(pool,{'generateAstralNetwork_eqSpaced.m','percCheck.m', ...
    'estPercProb_eqSpaced.m'})
F = parallel.FevalFuture.empty(Njobs,0);
job_idx = 1;
for idx = 1:numNetTypes
    for jdx = 1:numUniqueDensVals_byrow(idx)
        F(job_idx) = parfeval(pool,@estPercProb_eqSpaced,1, ...
            numAstersUsed{idx}(jdx),l,D,astralNumList(idx),Nsamp,useMEX);
        storeHere(job_idx,:) = [idx, jdx];
        job_idx = job_idx + 1;
    end
end
doneJobs = 0;
while doneJobs < Njobs
    try
        [completeIdx,p] = fetchNext(F);
        doneJobs = doneJobs + 1;
        p = reshape(p,[5,1]);
        net_idx = storeHere(completeIdx,1);
        dens_idx = storeHere(completeIdx,2);
        percProbs{net_idx}(1:5,dens_idx) = p;
        fprintf('Completed job %04i/%04i (%2.1f%%)\n',doneJobs,Njobs, ...
            doneJobs/Njobs*100)
        if dens_idx > 1 && dens_idx < numUniqueDensVals_byrow(net_idx)
            % note that z=0 if dens_idx-1 job has not been fetched
            z = prod(percProbs{net_idx}(1:5,(dens_idx-1):dens_idx)) == 1;
            if z
                fprintf('Skipping %i densities for astral number %i\n', ...
                    numUniqueDensVals_byrow(net_idx)-dens_idx, ...
                    astralNumList(net_idx));
                percProbs{net_idx}(1:5, ...
                    (dens_idx+1):numUniqueDensVals_byrow(net_idx)) = 1;
                toSkip = (storeHere(:,1) == net_idx).*(storeHere(:,2) > ...
                    dens_idx);
                doneJobs = doneJobs + sum(toSkip);
                cancel(F(logical(toSkip)))
            end
        end
    catch
        warning(['Failed to fetch an output. The corresponding future ' ...
            'may have been canceled.']);
    end
end
clear('F')
delete(pool)

for idx = 1:numNetTypes
    netLabel = sprintf('an%02i',astralNumList(idx));
    curves.(netLabel) = [actualDensities{idx}; percProbs{idx}];
end

if ~isempty(saveDirectory)
    fileToSave = fullfile(saveDirectory,filename + ".mat");
    save(fileToSave)    % saves all variables, including function call args
end
end