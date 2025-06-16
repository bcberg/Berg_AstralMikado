% timeseries_sysSize.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "timeseries_sysSize241012";

% Ubuntu paths
saveDir = '~/Documents/astral-mikado-data';
% Windows paths
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';

newReports = true;
if newReports
    dir = "~/Documents/AstralMikadoCYM/runs/timeseries_sysSize/";

    % ensure parameters match with timeseries_sysSize.cym.tpl & driver
    dens = 7.5;  % N * len_fil / A
    len_fil = 1;
    sysSizeList = [5; 10; 15; 20]; % side length of square region
    numSizes = length(sysSizeList);
    nFilList = dens * sysSizeList.^2 / len_fil; % may need to round this if numbers 
    nFilPerAster = 1;   % "original" Mikado network
    dataStartLine = 6; % CHECK - USING REPORTF
    ptsPerFiber = 6;
    frameVals = (0:105)';   % CHECK - USING REPORTF
    secsPerFrame = 1;       % CHECK REPORT FILES TO CONFIRM
    forceVal = 20;  % [pN]
    
    for sizeIdx = 1:numSizes
        sizeLabel = sprintf('s%03i',sysSizeList(sizeIdx));
        thisRun = sizeIdx - 1;
        cytoparams = struct('nFil',nFilList(sizeIdx),'nFilPerAster',...
            nFilPerAster,'dataStartLine',dataStartLine,...
            'ptsPerFiber',ptsPerFiber);
        CoM_timeseries.(sizeLabel) = parseTimeSeries(dir,thisRun, ...
            frameVals,cytoparams);
    end
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(fullfile(saveDir,'mat_files',simName + ".mat"), ...
        'nFilPerAster','dens','len_fil','nFilList','sysSizeList', ...
        'numSizes','frameVals','CoM_timeseries','forceVal','secsPerFrame')
    clearvars -except simName saveDir
end
% Ubuntu filepath
load(fullfile(saveDir,'mat_files',simName + ".mat"), ...
    'nFilPerAster','dens','len_fil','nFilList','sysSizeList', ...
    'numSizes','frameVals','CoM_timeseries','forceVal','secsPerFrame')

%% Plot timeseries

fig1 = figure(1);
set(fig1,'units','inches','Position',[1,1,6,8],'defaultTextInterpreter','tex')
fig1 = tiledlayout(fig1,2,1,'TileSpacing','compact');
nexttile(1)
hold on;
for sizeIdx = 1:numSizes
    sizeLabel = sprintf('s%03i',sysSizeList(sizeIdx));
    thisCoMhistory = CoM_timeseries.(sizeLabel);
    plot(frameVals*secsPerFrame, thisCoMhistory(:,1), 'LineWidth',0.75, ...
        'DisplayName', sizeLabel)
end
% indicate where force is switched on
ylimits = ylim(gca);
plot([5,5],ylimits,'--','Color',0.25*[1,1,1],'LineWidth',0.75)
hold off;
xlim('tight')
xlabel('Time [s]')
ylabel('x-center of mass [{\mu}m]')

nexttile(2)
hold on;
for sizeIdx = 1:numSizes
    sizeLabel = sprintf('s%03i',sysSizeList(sizeIdx));
    thisCoMhistory = CoM_timeseries.(sizeLabel);
    plot(frameVals*secsPerFrame, thisCoMhistory(:,2), 'LineWidth',0.75, ...
        'DisplayName', sizeLabel)
end
ylimits = ylim(gca);
plot([5,5],ylimits,'--','Color',0.25*[1,1,1],'LineWidth',0.75, ...
    'DisplayName',sprintf('%3.1f pN applied',forceVal))
hold off;
xlim('tight')
xlabel('Time [s]')
ylabel('y-center of mass [{\mu}m]')
lg = legend;
lg.Layout.Tile = 'east';
titletext = sprintf('a_n=%i, \\rho=%2.1f, L_{fil}=%1.1f, dt=0.01s', ...
    nFilPerAster,dens,len_fil);
title(fig1,titletext)

exportgraphics(fig1,fullfile(saveDir,'exploratory_figures', ...
    simName + ".png"),'Resolution',300)
