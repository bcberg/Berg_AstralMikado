% filLengthQuestions.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% At a given astral number, is there a difference between fewer+larger
% asters and more+smaller asters, fixing overall filament density?

astralNum = 2;
D = 50;
lOverD = [1e-2; 2e-2; 5e-2; 1e-1; 2e-1; 4e-1; 8e-1; 1];
numRatioVals = length(lOverD);
numDensVals = 50;
densityRange = logspace(-1,1,numDensVals);
% Note: values in numAstersArray may not be unique, especially for high l/D
numAstersArray = round(repmat(densityRange*D/astralNum, ...
    [numRatioVals,1]) ./ repmat(lOverD,[1,numDensVals]));
sampPerDensity = 2000;

numAstersUsed = cell(numRatioVals,1);
numUniqueDensVals_byrow = zeros(numRatioVals,1);
actualDensities = cell(numRatioVals,1);
percProbs_oneclass = cell(numRatioVals,1);
for idx = 1:numRatioVals
    numAstersUsed{idx} = unique(numAstersArray(idx,:));
    numUniqueDensVals_byrow(idx) = numel(numAstersUsed{idx});
    actualDensities{idx} = astralNum*numAstersUsed{idx} * lOverD(idx) / D;
    percProbs_oneclass{idx} = zeros(1,numUniqueDensVals_byrow(idx));
end

pool = parpool(8);  % adjust depending on machine
useMEX = true;      % adjust depending on machine
Njobs = sum(numUniqueDensVals_byrow);
F = parallel.FevalFuture.empty(Njobs,0);
F_partition = zeros(Njobs,2);
job_idx = 1;
for idx = 1:numRatioVals
    for jdx = 1:numUniqueDensVals_byrow(idx)
        F(job_idx) = parfeval(pool,@estPercProb,1, ...
            numAstersUsed{idx}(jdx),lOverD(idx)*D,D,astralNum, ...
            sampPerDensity,useMEX);
        F_partition(job_idx,1) = idx;  % record this future's ratio idx
        F_partition(job_idx,2) = jdx;  % record density sub-index, too
        job_idx = job_idx + 1;
    end
end
waitMessage = 'Waiting for FevalFutures to complete... (%3i/%3i)';
h = waitbar(0,sprintf(waitMessage,0,Njobs));
for idx = 1:Njobs
    [completedIdx,p] = fetchNext(F);
    whereToStore = F_partition(completedIdx,:);
    percProbs_oneclass{whereToStore(1)}(whereToStore(2)) = p;
    waitbar(idx/Njobs,h,sprintf(waitMessage,idx,Njobs));
end
delete(pool)
delete(h)
clear('F')
% Ubuntu path
save(['~/Documents/AstralMikadoCYM/data/', ...
    sprintf('lOverD_nAST%02i',astralNum),'.mat']);

% to stop evaluation of parfeval futures, run:
% cancel(F)
% or run:
% cancelAll(p.FevalQueue)

% plotting
fig1 = figure(1);
set(fig1,'defaultLineLineWidth',0.75)
cmap = colormap(turbo(numRatioVals));
hold on
for idx = 1:numRatioVals
    plot(actualDensities{idx},percProbs_oneclass{idx},'DisplayName', ...
        sprintf('$l/D = %1.3f$',lOverD(idx)),'Color',cmap(idx,:));
end
hold off
xscale('log')
xlabel('Filament density [$\mu m^{-1}$]')
ylabel('Percolation probability')
title(sprintf('$a_n = %i$',astralNum))
legend('Location','eastoutside')
% Ubuntu path
exportgraphics(fig1,['~/Documents/AstralMikadoCYM/data/', ...
    sprintf('lOverD_nAST%02i',astralNum),'.png'],'Resolution',300)
% Windows path
% exportgraphics(fig1,['C:\Users\bcber\Documents\AstralMikadoCYM\data\', ...
%     sprintf('lOverD_nAST%02i',astralNum), '.png'],'Resolution',300)

%% How does the ratio l/D matter at each astral number?
% and is there an l << D regime where it doesn't matter much?
% e.g., in Lee & Torquato 1990, used l = 1 and smallest system D = 5 (up to
% D = 40)

% astralNumList = 1:24;