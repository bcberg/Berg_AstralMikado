% percMCRG.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Importing data and fits

% Ubuntu path
% saveDir = "~/Documents/AstralMikadoCYM/data";
% Windows path
saveDir = "C:\Users\bcber\Documents\AstralMikadoCYM\data";
filePattern = "percProbs_l%02i_D%02i";

l = 1;
Dlist = [5,10,15,20];
numD = numel(Dlist);
allCurves = cell(1,numD);
allAstralNums = [];
for idx = 1:numD
    load(fullfile(saveDir,sprintf(filePattern,l,Dlist(idx))),'curves', ...
        'astralNumList')
    allCurves{idx} = curves;
    allAstralNums = union(allAstralNums,astralNumList);
    clear curves
end
numNetTypes = length(allAstralNums);
percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
numPercTypes = length(percTypes); % i.e., = 5

load(fullfile(saveDir,'percMCRG_nonparamFits.mat'),'smsplFits')
load(fullfile(saveDir,'percMCRG_paramFits.mat'),'sigFits')

%% Attempting MCRG (nonparametric fits)

thisPercType = 1;
testAstralNum = 14;
crossings = zeros(numD-1,1);
otherDers = zeros(numD-1,1);
% smallCurve = smsplFits{testAstralNum,thisPercType,1};
smallCurve = sigFits{testAstralNum,thisPercType,1};
for kdx = 1:(numD-1)
    % thisCurve = smsplFits{testAstralNum,thisPercType,kdx+1};
    thisCurve = sigFits{testAstralNum,thisPercType,kdx+1};
    fun = @(x) smallCurve(x) - thisCurve(x);
    [x,fval] = fsolve(fun,5);
    crossings(kdx) = x;
    otherDers(kdx) = differentiate(thisCurve,x);
end
smallDers = differentiate(smallCurve,crossings);
Lambda = otherDers ./ smallDers;

x = log(Dlist(2:end)/Dlist(1))';
y = log(Lambda);
linefit = fit(x,y,'poly1');

fig1 = figure(1);
plot(linefit,'-b',x,y,'*k')
xlabel("$\ln(b/b')$")
ylabel("$\ln \Lambda_{b,b'}$")
title(sprintf("$\\nu\\approx %2.4f$",1/linefit.p1))

%% trying to estimate nu for all astral numbers
nu = zeros(24,1);
for idx = 1:24
    crossings = zeros(numD-1,1);
    otherDers = zeros(numD-1,1);
    % smallCurve = smsplFits{idx,thisPercType,1};
    smallCurve = sigFits{idx,thisPercType,1};
    for kdx = 1:(numD-1)
        % thisCurve = smsplFits{idx,thisPercType,kdx+1};
        thisCurve = sigFits{idx,thisPercType,kdx+1};
        fun = @(x) smallCurve(x) - thisCurve(x);
        [x,fval] = fsolve(fun,5);
        crossings(kdx) = x;
        otherDers(kdx) = differentiate(thisCurve,x);
    end
    smallDers = differentiate(smallCurve,crossings);
    Lambda = otherDers ./ smallDers;

    x = log(Dlist(2:end)/Dlist(1))';
    y = log(Lambda);
    linefit = fit(x,y,'poly1');
    nu(idx) = 1/linefit.p1;
end

fig2 = figure(2);
hold on
plot(1:24,nu,'*b-')
plot([1,24],[4/3,4/3],'--k')
xlabel('Astral number')
ylabel('$\nu$ estimate')