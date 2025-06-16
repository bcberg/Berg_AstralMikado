% timeseries_astralNum.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "timeseries_astralNum250403";

% Ubuntu paths
saveDir = '~/Documents/astral-mikado-data';
% Windows paths
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';

newReports = false;
if newReports
    dir = "~/Documents/AstralMikadoCYM/runs/timeseries_astralNum/";

    % ensure parameters match with timeseries_check.cym.tpl & driver
    dens = 7.5;  % N * len_fil / A
    len_fil = 1;
    D = 10;
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = [1,4,8,16];   % "original" Mikado network
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6; % CHECK - USING REPORTF
    ptsPerFiber = 6;
    frameVals = (0:650)';   % CHECK - USING REPORTF
    secsPerFrame = 0.1;       % CHECK REPORT FILES TO CONFIRM
    shearForce = 5; % [pN]
    tensForce = 20; % [pN]

    netLabelPat = 'an%02d';
    % analyze shear runs
    for netIdx = 1:4
        netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
        thisRun = netIdx - 1;
        cytoparams = struct('nFil',nFil,'nFilPerAster', ...
            nFilPerAsterList(netIdx),'dataStartLine',dataStartLine, ...
            'ptsPerFiber',ptsPerFiber);
        shearCoMs.(netLabel) = parseTimeSeries(dir,thisRun, ...
            frameVals,cytoparams);
    end
    % analyze tensile runs
    for netIdx = 1:4
        netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
        thisRun = netIdx - 1 + nNetTypes;
        cytoparams = struct('nFil',nFil,'nFilPerAster', ...
            nFilPerAsterList(netIdx),'dataStartLine',dataStartLine, ...
            'ptsPerFiber',ptsPerFiber);
        tensCoMs.(netLabel) = parseTimeSeries(dir,thisRun, ...
            frameVals,cytoparams);
    end
    
    save(fullfile(saveDir,'mat_files',simName + ".mat"), ...
        'dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
        'frameVals','secsPerFrame','shearForce','tensForce', ...
        'shearCoMs','tensCoMs','netLabelPat')
    clearvars -except simName saveDir
end

load(fullfile(saveDir,'mat_files',simName + ".mat"), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce', ...
    'shearCoMs','tensCoMs','netLabelPat')

%% Plot timeseries

% shear
fig1 = figure(1); clf;
set(fig1,'units','centimeters','Position',[1,1,20,10], ...
    'defaultLineLineWidth',0.75)
fig1 = tiledlayout(fig1,1,2);
nexttile(1)
hold on
for netIdx = 1:nNetTypes
    netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs.(netLabel);
    plot(frameVals*secsPerFrame,theseCoMs(:,1),'DisplayName',netLabel)
end
ylims = ylim(gca);
plot([5,5], ylims, '--k', 'DisplayName',sprintf("%2.1f pN applied", ...
    shearForce))
plot([45,45], ylims, '-k', 'DisplayName', 'Final pos. sample')
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('x-center of mass [{\mu}m]')
lg = legend;
lg.Layout.Tile = 'east';

nexttile(2)
hold on
for netIdx = 1:nNetTypes
    netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs.(netLabel);
    plot(frameVals*secsPerFrame,theseCoMs(:,2))
end
ylims = ylim(gca);
plot([5,5], ylims, '--k')
plot([45,45], ylims, '-k')
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('y-center of mass [{\mu}m]')

exportgraphics(fig1,fullfile(saveDir,'exploratory_figures',simName + ...
    "shear.png"),'Resolution',300)

% tensile
fig2 = figure(2); clf;
set(fig2,'units','centimeters','Position',[1,1,20,10], ...
    'defaultLineLineWidth',0.75)
fig2 = tiledlayout(fig2,1,2);
nexttile(1)
hold on
for netIdx = 1:nNetTypes
    netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs.(netLabel);
    plot(frameVals*secsPerFrame,theseCoMs(:,1),'DisplayName',netLabel)
end
ylims = ylim(gca);
plot([5,5], ylims, '--k', 'DisplayName', sprintf("%2.1f pN applied", ...
    tensForce))
plot([45,45], ylims, '-k', 'DisplayName', 'Final pos. sample')
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('x-center of mass [{\mu}m]')
lg = legend;
lg.Layout.Tile = 'east';

nexttile(2)
hold on
for netIdx = 1:nNetTypes
    netLabel = sprintf(netLabelPat,nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs.(netLabel);
    plot(frameVals*secsPerFrame,theseCoMs(:,2))
end
ylims = ylim(gca);
plot([5,5], ylims, '--k')
plot([45,45], ylims, '-k')
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('y-center of mass [{\mu}m]')

exportgraphics(fig2,fullfile(saveDir,'exploratory_figures',simName + ...
    "tens.png"),'Resolution',300)

