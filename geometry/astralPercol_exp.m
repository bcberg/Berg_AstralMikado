% astralPercol_exp.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Parameters

recompute = input("Recompute percolation probabilities (1=yes, 0=no)? ");
if recompute
    D = input("Enter system size: ");
    l = input("Enter mean filament length: (default = 1)");
    if isempty(l)
        l = 1;
    end
    astralNumList = (1:24)';
    densSpec = [0,2+log10(2),50];
    Nsamp = 2000;
    Ncores = 8;
else
    % load a particular existing dataset
    D = 10;
    l = 1;
end
% Ubuntu path
dataDirectory = "/home/bcberg/Documents/astral-mikado-data/mat_files";
figDirectory = "/home/bcberg/Documents/astral-mikado-data/exploratory_figures";
% Windows path
% saveDirectory = "C:\Users\bcber\Documents\astral-mikado-data\mat_files";
filename = sprintf("percProbs_exp_l%02i_D%02i",l,D);

%% Percolation probability estimation

if recompute
    curves = getAllPercCurves_exp(l,D,astralNumList,densSpec,Nsamp,Ncores, ...
        dataDirectory);
end
load(fullfile(dataDirectory,filename + ".mat"))

%% Curves for all percolation definitions

fig1 = figure(1);
set(fig1,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
    'defaultLineLineWidth',0.6)
fig1 = tiledlayout(3,2,'TileSpacing','compact');

percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
numPercTypes = length(percTypes); % i.e., = 5
cmap = colormap(turbo(numNetTypes));
for pdx = 1:numPercTypes
    nexttile(fig1,mod(pdx,6) + 6 * (mod(pdx,6)==0))
    hold on
    for idx = 1:numNetTypes
        plot(actualDensities{idx},percProbs{idx}(pdx,:),'DisplayName', ...
            sprintf('an%02i',astralNumList(idx)),'Color',cmap(idx,:))
    end
    hold off
    xscale('log')    % 'linear' or 'log'
    xlabel('Filament density [$\mu m^{-1}$]')
    ylabel('Percolation probability')
    title(percTypes{pdx})
    if pdx == 1
        lg = legend('NumColumns',4,'FontSize',12);
        lg.Layout.Tile = 6;
    end
end
title(fig1,sprintf('Domain size: %2.1f, Filament length: %1.1f',D,l))
exportgraphics(fig1,fullfile(figDirectory,filename + "_curves.pdf"), ...
    'Resolution',300)

%% Heatmaps for all percolation definitions

fig2 = figure(2);
set(fig2,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
    'defaultLineLineWidth',0.6)
fig2 = tiledlayout(3,2,'TileSpacing','compact');

percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
numPercTypes = length(percTypes); % i.e., = 5
% cmap = colormap(parula(numNetTypes));
for pdx = 1:numPercTypes
    nexttile(fig2,mod(pdx,6) + 6 * (mod(pdx,6)==0))
    % For heatmap, need to interpolate values at each value of densityRange
    % Currently using cubic splines
    numDensVals = length(densityRange);
    heatmap = zeros(numNetTypes,numDensVals);
    for idx = 1:numNetTypes
        heatmap(idx,:) = clip(spline(log10(actualDensities{idx}), ...
            percProbs{idx}(pdx,:),log10(densityRange)),0,1);
    end

    pcolor(densityRange,astralNumList,heatmap)
    xscale('log')
    xlabel('Filament density [$\mu m^{-1}$]')
    ylabel('Astral number')
    title(percTypes{pdx})
    if pdx == 1
        bar = colorbar;
        bar.Layout.Tile = 6;
    end
end
title(fig2,sprintf('Domain size: %2.1f, Filament length: %1.1f',D,l))
exportgraphics(fig2,fullfile(figDirectory,filename + "_heat.pdf"), ...
    'Resolution',300)