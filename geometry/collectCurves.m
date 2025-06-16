% collectCurves.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Filepaths, parameters

directory = "~/Documents/AstralMikadoCYM/data/";
filepattern = "run240922/percProbs_l%02i_D%02i_an%02i";
collectedFilepattern = "percProbs_l%02i_D%02i";
astralNumList = 1:24;
numNetTypes = length(astralNumList);
l = 1;
D = 50;

%% Collect individual curves

actualDensities_all = cell(numNetTypes,1);
numAstersUsed_all = cell(numNetTypes,1);
numUniqueDensVals_byrow = zeros(numNetTypes,1);
percProbs_all = cell(numNetTypes,1);
for idx = 1:numNetTypes
    thisFile = directory + sprintf(filepattern,l,D,idx);
    load(thisFile,'densityRange','actualDensities','Nsamp','numAstersUsed', ...
        'numUniqueDensVals','percProbs');
    actualDensities_all{idx} = actualDensities;
    numAstersUsed_all{idx} = numAstersUsed;
    numUniqueDensVals_byrow(idx) = numUniqueDensVals;
    percProbs_all{idx} = percProbs;
    clear actualDensities numAstersUsed numUniqueDensVals percProbs
end
actualDensities = actualDensities_all;
numAstersUsed = numAstersUsed_all;
percProbs = percProbs_all;
clear actualDensities_all numAstersUsed_all percProbs_all thisFile idx
fileToSave = directory + sprintf(collectedFilepattern,l,D) + ".mat";
save(fileToSave)