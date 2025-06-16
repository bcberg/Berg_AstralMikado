% oneVsTwoHalves.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% How does N*l/D^2 compare to 2*N*(l/2)/D^2?
% i.e., classical Mikado with length l vs astral number 2 with length l/2

resample = false;
if resample
    numDensVals = 50;
    % Note: there may be fewer density values in certain percolation curves due
    % to requiring integer numbers of asters (density resolution is coarser at
    % higher astral number)
    sampPerDensity = 2000;
    l = 5;
    D = 50;
    densityRange = logspace(-1,1,numDensVals);
    numFilRange = round(densityRange * (D^2 / l));
    actualDensities = numFilRange * l / D^2;  % or, 2*numFilRange*(l/2)/D^2
    percProbs_clas = zeros(1,numDensVals);
    percProbs_astr = zeros(1,numDensVals);

    pool = parpool(4);
    parfor idx = 1:numDensVals
        numFil = numFilRange(idx);
        percCount_clas = 0;
        percCount_astr = 0;
        for jdx = 1:sampPerDensity
            [net_clas,cross_clas,asters_clas] = generateAstralNetwork(numFil, ...
                l, D, 1, true);
            [net_astr,cross_astr,asters_astr] = generateAstralNetwork(numFil, ...
                l/2, D, 2, true);
            [percTF_clas,~] = percCheck(cross_clas,net_clas.nodes,D);
            [percTF_astr,~] = percCheck(cross_astr,net_astr.nodes,D);
            if percTF_clas(1) || percTF_clas(2)
                percCount_clas = percCount_clas + 1;
            end
            if percTF_astr(1) || percTF_astr(2)
                percCount_astr = percCount_astr + 1;
            end
        end
        percProbs_clas(idx) = percCount_clas / sampPerDensity;
        percProbs_astr(idx) = percCount_astr / sampPerDensity;
    end
    delete(pool);
    % Ubuntu path
    % save('~/Documents/AstralMikadoCYM/data/oneVsTwoHalves.mat')
    % Windows path
    % save('C:\Users\bcber\Documents\AstralMikadoCYM\data\oneVsTwoHalves.mat')
end

%% Plotting

load('C:\Users\bcber\Documents\AstralMikadoCYM\data\oneVsTwoHalves.mat')
fig1 = figure(1);
set(fig1,'defaultLineLineWidth',0.75)
label1 = sprintf('$a_n = 1$, $l = %1.1f$',l);
label2 = sprintf('$a_n = 2$, $l = %1.1f$',l/2);
hold on
plot(actualDensities,percProbs_clas,'-b','DisplayName',label1)
plot(actualDensities,percProbs_astr,'--r','DisplayName',label2)
hold off
xscale('log')
xlabel('Filament density [$\mu m^{-1}$]')
ylabel('Percolation probability')
title([sprintf('$D = %2.1f',D), '\;\mu m$'])
legend('Location','southeast')
% Windows path
exportgraphics(fig1,['C:\Users\bcber\Documents\AstralMikadoCYM\data\', ...
    'oneVsTwoHalves.png'],'Resolution',300);