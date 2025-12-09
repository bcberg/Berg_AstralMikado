% figures_manu.m
% Brady Berg, 03/18/2025
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Importing data

% Ubuntu paths
saveDir = '~/Documents/astral-mikado-data';
dataFolder = '~/Documents/astral-mikado-data/mat_files';
figureFolder = '~/Documents/astral-mikado-data/subfigures';
% astralFuncFolder = '~/Documents/astral-mikado/';
% Windows paths
% dataFolder = 'C:\Users\bcber\Documents\astral-mikado-data\mat_files';
% figureFolder = 'C:\Users\bcber\Documents\astral-mikado-data\subfigures';
% astralFuncFolder = 'C:\Users\bcber\Documents\astral-mikado\';

filetype = '.pdf';

% shear modulus sweep
load(fullfile(dataFolder,'rep_shear241105_l01_D10.mat'),'shearModuli', ...
    'forceVals','nFilPerAsterList','nNetTypes','numRep','dens','disps')
sweep_s = struct('forceVals',forceVals,'nFilPerAsterList', ...
    nFilPerAsterList,'nNetTypes',nNetTypes,'numRep',numRep,'dens',dens);
disps_s = disps;
moduli_s = shearModuli;
clear('forceVals','nFilPerAsterList','nNetTypes','numRep','shearModuli',...
    'dens','disps')

% tensile modulus sweep
load(fullfile(dataFolder,'rep_tens241116_l01_D10.mat'),'tensModuli', ...
    'forceVals','nFilPerAsterList','nNetTypes','numRep','dens','disps')
sweep_t = struct('forceVals',forceVals,'nFilPerAsterList', ...
    nFilPerAsterList,'nNetTypes',nNetTypes,'numRep',numRep,'dens',dens);
disps_t = disps;
moduli_t = tensModuli;
clear('forceVals','nFilPerAsterList','nNetTypes','numRep','tensModuli', ...
    'dens','disps')

%%%% negative outliers - currently clipping to 0 %%%%%
moduli_s(23,2) = 0;     % fit yielded -6.267
moduli_t(14,22) = 0;    % fit yielded -289.8, data looks nonpercolated
moduli_t(20,19) = 0;    % fit yielded -246.8
%%%%%

% modulus distribution sampling
load(fullfile(dataFolder,'mod_distr241120_l01_D10.mat'),'shearModuli', ...
    'tensModuli','nFilPerAsterList','nNetTypes','numRep','shearForces', ...
    'tensForces')
moduli_s_hiN = shearModuli;
moduli_t_hiN = tensModuli;
modSampling = struct('nFilPerAsterList',nFilPerAsterList,'nNetTypes', ...
    nNetTypes,'numRep',numRep,'shearForces',shearForces,'tensForces', ...
    tensForces);
clear('shearModuli','tensModuli','nNetTypes','numRep','shearForces', ...
    'tensForces','nFilPerAsterList')

% example astral networks in MATLAB
load(fullfile(dataFolder,'exNetTBSpan.mat'),'astralNum','rho','l','D', ...
    'targetFilNum','numAsters','network','crossings','asters','percTF', ...
    'percStats')
tbSpanParams = struct('astralNum',astralNum,'rho',rho,'l',l,'D',D, ...
    'targetFilNum',targetFilNum,'numAsters',numAsters);
tbNetwork = network;
tbCrossings = crossings;
tbAsters = asters;
tbPercTF = percTF;
tbPercStats = percStats;
clear('astralNum','rho','l','D','targetFilNum','numAsters','network', ...
    'crossings','asters','percTF','percStats')

load(fullfile(dataFolder,'exNetConn.mat'),'astralNum','rho','l','D', ...
    'targetFilNum','numAsters','network','crossings','asters','percTF', ...
    'percStats')
connParams = struct('astralNum',astralNum,'rho',rho,'l',l,'D',D, ...
    'targetFilNum',targetFilNum,'numAsters',numAsters);
connNetwork = network;
connCrossings = crossings;
connAsters = asters;
connPercTF = percTF;
connPercStats = percStats;
clear('astralNum','rho','l','D','targetFilNum','numAsters','network', ...
    'crossings','asters','percTF','percStats')

% percolation data
load(fullfile(dataFolder,sprintf("percProbs_l%02i_D%02i",1,10)), ...
    'curves','densityRange')

% nonparametric (smoothing spline) fits to percolation data
load(fullfile(dataFolder,'percMCRG_nonparamFits_tzOnly.mat'), ...
    'smsplFits','Dlist','percDataUsed')

% dangling ends data
load(fullfile(dataFolder,'danglingEnds_250319.mat'),'numNets', ...
        'nFilPerAsterList','numNetTypes','endsLengths','endsMeans', ...
        'stdOfMeans','endsMeansTrunc','stdOfMeansTrunc','rho','l','D', ...
        'targetFilNum','labelPat','useFraction')
danglingParams = struct('nFilPerAsterList',nFilPerAsterList,'l',l, ...
    'D',D,'rho',rho,'numNets',numNets,'numNetTypes',numNetTypes, ...
    'targetFilNum',targetFilNum,'labelPat',labelPat);
clear('nFilPerAsterList','l','D','rho','numNets','numNetTypes', ...
    'targetFilNum','labelPat')

% filament bending rigidity sweep
load(fullfile(dataFolder,"kbend_shear250228.mat"),'nFilPerAsterList', ...
    'forceVals','numRep','nKbend','nNetTypes','kbend_list','shearModuli')
kbendParams = struct('nFilPerAsterList',nFilPerAsterList,'numRep',numRep, ...
    'forceVals',forceVals,'nNetTypes',nNetTypes,'kbend_list',kbend_list, ...
    'nKbend',nKbend);
moduli_s_kbend = shearModuli;
clear('nFilPerAsterList','numRep','forceVals','nNetTypes','kbend_list', ...
    'nKbend','shearModuli')

% angular stiffness at astral centers sweep
load(fullfile(dataFolder,"krot_shear250226.mat"),'nFilPerAsterList', ...
    'forceVals','numRep','nKrot','nNetTypes','krot_list','shearModuli')
krotParams = struct('nFilPerAsterList',nFilPerAsterList,'numRep',numRep,...
    'forceVals',forceVals,'nNetTypes',nNetTypes,'krot_list',krot_list, ...
    'nKrot',nKrot);
moduli_s_krot = shearModuli;
clear('nFilPerAsterList','numRep','forceVals','nNetTypes','krot_list', ...
    'nKrot','shearModuli')

% crosslinker density sweep
load(fullfile(dataFolder,"crosslinker_dens_shear250320.mat"), ...
    'nFilPerAsterList','forceVals','numRep','nNetTypes', ...
    'crosslinker_dens_list','nCrslnkDens','shearModuli')
crslnkDensParams = struct('nFilPerAsterList',nFilPerAsterList,'numRep',...
    numRep,'forceVals',forceVals,'nNetTypes',nNetTypes,'nCrslnkDens', ...
    nCrslnkDens,'crosslinker_dens_list',crosslinker_dens_list);
moduli_s_crslnkDens = shearModuli;
clear('nFilPerAsterList','numRep','forceVals','nNetTypes','nCrslnkDens',...
    'shearModuli','crosslinker_dens_list')

% crosslinker stiffness sweep
load(fullfile(dataFolder,"crosslinker_stiff_shear250401.mat"), ...
    'nFilPerAsterList','forceVals','numRep','nNetTypes', ...
    'crosslinker_stiff_list','nCrslnkStiff','shearModuli')
crslnkStiffParams = struct('nFilPerAsterList',nFilPerAsterList,'numRep',...
    numRep,'forceVals',forceVals,'nNetTypes',nNetTypes,'nCrslnkStiff', ...
    nCrslnkStiff,'crosslinker_stiff_list',crosslinker_stiff_list);
moduli_s_crslnkStiff = shearModuli;
clear('nFilPerAsterList','numRep','forceVals','nNetTypes','nCrslnkStiff',...
    'shearModuli','crosslinker_stiff_list')

% shear modulus vs. astral number at additional filament densities
load(fullfile(dataFolder,"rep_shear_otherRhos250324.mat"),'dens','len_fil', ...
    'D','nFil','nFilPerAsterList','forceVals','numRep','nNetTypes', ...
    'numDens','shearModuli')
moduli_s_otherRhos = shearModuli;
otherRhoParams = struct('dens',dens,'nFilPerAsterList',nFilPerAsterList, ...
    'nFil',nFil,'len_fil',len_fil,'D',D,'forceVals',forceVals,'numRep', ...
    numRep,'nNetTypes',nNetTypes,'numDens',numDens);
clear('dens','nFilPerAsterList','nFil','len_fil','D','forceVals','numRep', ...
    'nNetTypes','numDens')

% "Mikado rigidity" (Wilhelm & Frey 2003)
load(fullfile(dataFolder,"mikado_rigidity250320.mat"),'D_list','numD',...
        'dens_list','numDens','forceVals','numRep','shearModuli')
moduli_s_mikado = shearModuli;
mikadoParams = struct('D_list',D_list,'numD',numD,'dens_list',dens_list,...
    'numDens',numDens,'forceVals',forceVals,'numRep',numRep);
clear('D_list','numD','dens_list','numDens','forceVals','numRep', ...
    'shearModuli')

% timeseries at various astral numbers
load(fullfile(dataFolder,"timeseries_astralNum250403.mat"), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce', ...
    'shearCoMs','tensCoMs','netLabelPat')
timeserParams = struct('dens',dens,'len_fil',len_fil,'D',D,'nFil',nFil, ...
    'nFilPerAsterList',nFilPerAsterList,'nNetTypes',nNetTypes, ...
    'frameVals',frameVals,'secsPerFrame',secsPerFrame,'shearForce', ...
    shearForce,'tensForce',tensForce,'netLabelPat',netLabelPat);
clear('dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce','netLabelPat')

% counting productive segments per node
load(fullfile(dataFolder,"springPerNode.mat"),'l','D', ...
        'rho','targetFilNum','an_list','nNetTypes','numAsters_list','N',...
        'springPerNodeCounts','springPerNodeMeans','stdev', ...
        'means_by_net_type','ci_95','selected_an','num_selected', ...
        'aggregatedCounts')
springCountParams = struct('l',l,'D',D,'rho',rho,'targetFilNum', ...
    targetFilNum,'an_list',an_list,'nNetTypes',nNetTypes, ...
    'numAsters_list',numAsters_list,'N',N,'selected_an',selected_an, ...
    'num_selected',num_selected);
springPerNode_stdev = stdev;
springPerNode_means_by_net_type = means_by_net_type;
springPerNode_ci_95 = ci_95;
clear('l','D','rho','targetFilNum','an_list','nNetTypes', ...
    'numAsters_list','N','selected_an','num_selected','stdev', ...
    'means_by_net_type','ci_95')

%% General figure settings

%%%%% font specifications %%%%%
fsz = 22;
fsz2 = 18;
insetFsz = 16;
fname = 'CMU Serif';

%%%%% standard colors %%%%%
colors = colororder("gem12");

%%%%% for 95% CIs %%%%%
z = norminv(0.975);

%%%%% CUSTOM COLORS %%%%%
lgray = [197,190,181]/255;
gray = [85,87,89]/255;
uciBlue = [0,100,164]/255;
uciGold = [255,210,0]/255;
uciOrange = [247,141,45]/255;

%% MODULI VS. ASTRAL NUMBER %%

[stdevs_s,means_s] = std(moduli_s,0,2); % compute summary stats across cols
[stdevs_t,means_t] = std(moduli_t,0,2); % compute summary stats across cols

modPlot = figure(1);
set(modPlot,'units','centimeters','Position',[1,1,20,15])
modPlot = tiledlayout(modPlot,2,1,'TileSpacing','loose');
nexttile(1)
hold on
for idx = 1:sweep_s.nNetTypes
    plot(sweep_s.nFilPerAsterList(idx),moduli_s(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
errorbar(sweep_s.nFilPerAsterList, means_s, ...
    z*stdevs_s/sqrt(sweep_s.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,sweep_s.nFilPerAsterList(end)])
xticks(0:4:sweep_s.nFilPerAsterList(end))
% xlabel('Astral number $a_n$')
ylabel({'Shear modulus'; '$G$ [${\rm pN}\cdot\mu{\rm m^{-1}}$]'})

nexttile(2)
hold on
for idx = 1:sweep_t.nNetTypes
    plot(sweep_t.nFilPerAsterList(idx),moduli_t(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
errorbar(sweep_t.nFilPerAsterList, means_t, ...
    z*stdevs_t/sqrt(sweep_t.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,sweep_t.nFilPerAsterList(end)])
xticks(0:4:sweep_t.nFilPerAsterList(end))
xlabel('Astral number $a_n$')
ylim([0,450])   % excludes one high outlier @ a_n = 18
% ylabel({"Young's modulus"; "$Y$ [${\rm pN}\cdot\mu{\rm m^{-1}}$]"})
% exportgraphics(modPlot,fullfile(figureFolder,['moduliVsAstralNum', ...
%     filetype]),'ContentType','vector')

%% MODULUS DISTRIBUTION PLOTS %%

modDistPlot = figure(2); clf;
set(modDistPlot,'units','centimeters','Position',[1,1,15,15])
modDistPlot = tiledlayout(modDistPlot,2,1,'TileSpacing','loose');
% legLabels = arrayfun(@(x) sprintf('$a_n=%i$',x), ...
%     modSampling.nFilPerAsterList,'UniformOutput',false);
% modSampling.nFilPerAsterList is [1;2;3;4;8;12;16;20;24]
% chosenTypes = [1,4,5,7];    % indices of astral nums 1,4,8,16
chosenTypes = [1,2,4,5,7];  % indices of astral nums 1,2,4,8,16
numCurves = length(chosenTypes);
distColors = colormap(parula(numCurves+1));
nexttile(1)
hold on
for idx = 1:numCurves
    l = cdfplot(moduli_s_hiN(chosenTypes(idx),:));
    l.LineWidth = 2;
    l.Color = distColors(idx,:);
end
xl = xlim;
plot([0,xl(2)],0.5*[1,1],'--','LineWidth',1.5,'Color','k')
hold off
xlabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
xlim([0,xl(2)])
ylabel('')
title('')
set(gca,'FontName',fname,'FontSize',fsz)
nexttile(2)
hold on
for idx = 1:numCurves
    l = cdfplot(moduli_t_hiN(chosenTypes(idx),:));
    l.LineWidth = 2;
    l.Color = distColors(idx,:);
end
xl = xlim;
plot([0,xl(2)],0.5*[1,1],'--','LineWidth',1.5,'Color','k')
hold off
xlabel("Young's modulus $Y$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]")
xlim([0,xl(2)])
ylabel('')
title('')
set(gca,'FontName',fname,'FontSize',fsz)
ylabel(modDistPlot,'Empirical CDF','Interpreter','latex','FontSize', ...
    fsz,'FontName',fname)
lg = legend(string(modSampling.nFilPerAsterList(chosenTypes)),'Location','east', ...
    'FontSize',fsz);
title(lg,'$a_n$')
lg.Layout.Tile = 'east';
% exportgraphics(modDistPlot,fullfile(figureFolder,['modECDFs',filetype]),...
%     'ContentType','vector')

%% EXAMPLE PROBABILITY DENSITIES (HISTOGRAMS OF MODULI) %%

% modSampling.nFilPerAsterList is [1;2;3;4;8;12;16;20;24]

shear1hist = figure(3); clf;
set(shear1hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_s_hiN(1,:),'Normalization','probability');
h.FaceColor = distColors(1,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',insetFsz)
% ylabel('Prob. density')
% exportgraphics(shear1hist,fullfile(figureFolder,['shear1hist',filetype]), ...
%     'ContentType','vector')

shear16hist = figure(4); clf;
set(shear16hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_s_hiN(7,:),'Normalization','probability');
h.FaceColor = distColors(5,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',insetFsz)
% ylabel('Prob. density')
% exportgraphics(shear16hist,fullfile(figureFolder,['shear16hist',filetype]), ...
%     'ContentType','vector')

tens1hist = figure(5); clf;
set(tens1hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_t_hiN(1,:),'Normalization','probability');
h.FaceColor = distColors(1,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',insetFsz)
% ylabel('Prob. density')
% exportgraphics(tens1hist,fullfile(figureFolder,['tens1hist',filetype]), ...
%     'ContentType','vector')

tens16hist = figure(6); clf;
set(tens16hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_t_hiN(7,:),'Normalization','probability');
h.FaceColor = distColors(5,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',insetFsz)
% ylabel('Prob. density')
% exportgraphics(tens16hist,fullfile(figureFolder,['tens16hist',filetype]), ...
%     'ContentType','vector')

%% PROBABILITY OF RELATIVE WEAKNESS (<0.5*MEDIAN) %%

alpha = 0.5;
med_s = median(moduli_s_hiN,2);
med_t = median(moduli_t_hiN,2);
probPartition_s = zeros(3,modSampling.nNetTypes);
probPartition_t = zeros(3,modSampling.nNetTypes);
for idx = 1:modSampling.nNetTypes
    probPartition_s(1,idx) = sum(moduli_s_hiN(idx,:) < alpha*med_s(idx)) / ...
        modSampling.numRep;
    probPartition_s(3,idx) = sum(moduli_s_hiN(idx,:) > (1/alpha)*med_s(idx)) /...
        modSampling.numRep;
    probPartition_t(1,idx) = sum(moduli_t_hiN(idx,:) < alpha*med_t(idx)) / ...
        modSampling.numRep;
    probPartition_t(3,idx) = sum(moduli_t_hiN(idx,:) > (1/alpha)*med_t(idx))/...
        modSampling.numRep;
end
probPartition_s(2,:) = ones(1,modSampling.nNetTypes) - ...
    probPartition_s(1,:) - probPartition_s(3,:);
probPartition_t(2,:) = ones(1,modSampling.nNetTypes) - ...
    probPartition_t(1,:) - probPartition_t(3,:);

relativelyWeak = figure(7); clf;
set(relativelyWeak,'units','centimeters','InnerPosition',[1,1,12,16], ...
    'defaultLineLineWidth',2)
relativelyWeak = tiledlayout(relativelyWeak,2,1,'TileSpacing','loose');
nexttile(1)
set(gca,'FontName',fname,'FontSize',fsz)
hold on
plot(modSampling.nFilPerAsterList(1:7),probPartition_s(1,1:7),'o-', ...
    'DisplayName','Weak','Color',colors(1,:),'MarkerFaceColor', ...
    colors(1,:),'MarkerEdgeColor',colors(1,:))
text(0.5,0.36,sprintf('$\\Pr\\{G < %1.1f\\cdot G_{\\rm med}\\}$',alpha), ...
    'Color',colors(1,:),'FontSize',20)
% plot(modSampling.nFilPerAsterList(1:7),probPartition_s(2,1:7),'o:', ...
%     'DisplayName','Near median','Color','k','MarkerFaceColor','k', ...
%     'MarkerEdgeColor','k')
% text(7,0.8,sprintf('$ %1.1f\\cdot G_{50}\\leq G\\leq %1.1f\\cdot G_{50}$', ...
%     alpha,1/alpha),'Color','k','FontSize',16)
% plot(modSampling.nFilPerAsterList(1:7),probPartition_s(3,1:7),'o--', ...
%     'DisplayName','Strong','Color',colors(2,:),'MarkerFaceColor', ...
%     colors(2,:),'MarkerEdgeColor',colors(2,:))
% text(12,0.075,{'Strong';sprintf('($G > %1.1f\\cdot G_{50}$)',1/alpha)},'Color', ...
%     colors(2,:),'FontSize',20,'HorizontalAlignment','center')
hold off
xticks(0:4:16)
xlim([0,16])
% xlabel('Astral number $a_n$')
ylim([0,0.4])
ylabel('Shear')
% legend('Location','eastoutside')

nexttile(2)
set(gca,'FontName',fname,'FontSize',fsz)
hold on
plot(modSampling.nFilPerAsterList(1:7),probPartition_t(1,1:7),'o-', ...
    'DisplayName','Weak','Color',colors(1,:),'MarkerFaceColor', ...
    colors(1,:),'MarkerEdgeColor',colors(1,:))
text(0.5,0.36,sprintf('$\\Pr\\{Y < %1.1f\\cdot Y_{\\rm med}\\}$',alpha),'Color', ...
    colors(1,:),'FontSize',20)
% plot(modSampling.nFilPerAsterList(1:7),probPartition_t(2,1:7),'o:', ...
%     'DisplayName','Near median','Color','k','MarkerFaceColor','k', ...
%     'MarkerEdgeColor','k')
% text(7,0.85,sprintf('$ %1.1f\\cdot Y_{50}\\leq Y\\leq %1.1f\\cdot Y_{50}$', ...
%     alpha,1/alpha),'Color','k','FontSize',16)
% plot(modSampling.nFilPerAsterList(1:7),probPartition_t(3,1:7),'o--', ...
%     'DisplayName','Strong','Color',colors(2,:),'MarkerFaceColor', ...
%     colors(2,:),'MarkerEdgeColor',colors(2,:))
% text(12,0.075,{'Strong';sprintf('($Y > %1.1f\\cdot Y_{50}$)',1/alpha)},'Color', ...
%     colors(2,:),'FontSize',20,'HorizontalAlignment','center')
hold off
xticks(0:4:16)
xlim([0,16])
xlabel('Astral number $a_n$')
ylim([0,0.4])
ylabel("Young's")
ylabel(relativelyWeak,'Proportion weaker than half of median', ...
    'FontName',fname,'FontSize',fsz,'Interpreter','latex')
% lg = legend('NumColumns',2);
% lg.Layout.Tile = 'south';
% exportgraphics(relativelyWeak,fullfile(figureFolder,['relativelyWeak', ...
%     filetype]),'ContentType','vector')

%% DANGLING ENDS DISTRIBUTIONS %%

endsHists = figure(8); clf;
set(endsHists,'units','centimeters','Position',[1,1,15,15],'visible','on')
endsHists = tiledlayout(2,2,'TileSpacing','compact');
tilesPerPage = 4;

selectedAstralNums = [1,4,8,16];
for idx = 1:length(selectedAstralNums)
    % nexttile(endsHists,mod(idx,tilesPerPage) + tilesPerPage * (mod(idx,tilesPerPage)==0))
    nexttile(idx)
    astralNum = selectedAstralNums(idx);
    netLabel = sprintf(danglingParams.labelPat,astralNum);
    theseEnds = endsLengths.(netLabel);
    % hold on
    h = histogram(theseEnds,'Normalization','probability', ...
        'BinWidth',0.1);
    % h.FaceColor = 'k';
    % h.FaceAlpha = 1;
    % hTrunc = histogram(theseEnds(theseEnds<danglingParams.l), ...
    %     'Normalization','count','BinWidth',0.1);
    % hTrunc.FaceColor = uciBlue;
    % hTrunc.FaceAlpha = 1;
    ylim([0,0.4])
    title(sprintf('$a_n = %i$',astralNum))
    set(gca,'FontName',fname,'FontSize',18)
    % hold off
    % ylims for 'count' normalization
    % if astralNum == 1
    %     ylim([0,6000])
    % else
    %     ylim([0,2500])
    % end
end
xlabel(endsHists,'Dangling end length [$\mu\rm{m}$]','Interpreter', ...
    'latex','FontSize',18)
ylabel(endsHists,'Probability','FontName',fname,'FontSize',18)
% lg = legend('All ends','Ends $<\ell$','Interpreter','latex');
% lg.Layout.Tile = 'east';
% exportgraphics(endsHists,fullfile(figureFolder,['endsHists',filetype]),...
%     'ContentType','vector')

%% MEAN AND TOTAL DANGLING ENDS %%

linewidth = 1.5;
% useFraction statistics
[useSTD,useMean] = std(useFraction,0,2);

endsSummary = figure(9); clf;
set(endsSummary,'units','centimeters','Position',[1,1,15,12])
colororder([0 0 0; 0.8500 0.3250 0.0980])
yyaxis left
hold on
l1 = errorbar(danglingParams.nFilPerAsterList,mean(endsMeans,2), ...
    z*stdOfMeans/sqrt(danglingParams.numNets),'LineStyle','-','Color','k', ...
    'LineWidth',linewidth);
l2 = errorbar(danglingParams.nFilPerAsterList,mean(endsMeansTrunc,2), ...
    z*stdOfMeansTrunc/sqrt(danglingParams.numNets),'LineStyle','-','Color', ...
    colors(3,:),'LineWidth',linewidth);
hold off
yl = ylim;
ylim([0,yl(2)])
ylabel('Dangling end length [$\mu\rm{m}$]')
yyaxis right
l3 = errorbar(danglingParams.nFilPerAsterList,useMean, ...
    z*useSTD/sqrt(danglingParams.numNets),'LineWidth',linewidth);
ylim([0,1])
ylabel('Use fraction')
xlabel('Astral number $a_n$')
xlim([danglingParams.nFilPerAsterList(1)-0.25, ...
    danglingParams.nFilPerAsterList(end)+0.25])
xticks([1,4:4:24])
legend([l1,l2,l3],{'All ends','Ends $<\ell$','Use fraction'}, ...
    'Location','south')
set(gca,'FontName',fname,'FontSize',18)
% exportgraphics(endsSummary,fullfile(figureFolder,['endsSummary',filetype]), ...
%     'ContentType','vector')

%% PERCOLATION DEFINITION VISUALS %%

tbSpanEx = figure(10);
set(tbSpanEx,'units','centimeters','Position',[1,1,16,12])
hold on
plot([0,tbSpanParams.D],[0,0],'-','LineWidth',4,'Color',uciBlue)
plot([0,tbSpanParams.D],[tbSpanParams.D,tbSpanParams.D],'-', ...
    'LineWidth',4,'Color',uciBlue)
for idx = 1:tbSpanParams.numAsters
    for jdx = 1:tbSpanParams.astralNum
        plot(tbAsters.centers(idx,1) + tbSpanParams.l * cos(tbAsters.orients(idx,jdx)) * (0:1), ...
            tbAsters.centers(idx,2) + tbSpanParams.l * sin(tbAsters.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',3,'Color','k')
    end
end
for idx = 1:size(tbNetwork.springs,1)
    augNodeA = tbNetwork.springs(idx,1);
    augNodeB = tbNetwork.springs(idx,2);
    coords = [tbNetwork.augNodes(augNodeA,1:2);
        tbNetwork.augNodes(augNodeB,1:2)];
    plot(coords(:,1), coords(:,2),'.-','MarkerSize',18,'LineWidth',3,...
        'Color',uciOrange)
    % if augNodeA <= tbSpanParams.numAsters && tbSpanParams.astralNum >= 2
    %     % i.e., if augNodeA is an astral center
    %     plot(tbNetwork.augNodes(augNodeA,1),tbNetwork.augNodes(augNodeA,2), ...
    %         '.r','MarkerSize',24)
    % elseif augNodeB <= tbSpanParams.numAsters && tbSpanParams.astralNum >= 2
    %     % likewise for augNodeB
    %     plot(tbNetwork.augNodes(augNodeB,1),tbNetwork.augNodes(augNodeB,2), ...
    %         '.r','MarkerSize',24)
    % end
end
if tbSpanParams.astralNum >= 2
    plot(tbAsters.centers(:,1),tbAsters.centers(:,2),'.r','MarkerSize',24)
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - tbSpanParams.D/2)) * [-1,1] + tbSpanParams.D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - tbSpanParams.D/2)) * [-1,1] + tbSpanParams.D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
% exportgraphics(tbSpanEx,fullfile(figureFolder,['tbSpanEx_rnd',filetype]), ...
%     'ContentType','vector')

connEx = figure(11);
set(connEx,'units','centimeters','Position',[1,1,16,12])
hold on
plot([0,connParams.D],[0,0],'-','LineWidth',4,'Color',uciBlue)
plot([0,connParams.D],[connParams.D,connParams.D],'-','LineWidth', ...
    4,'Color',uciBlue)
for idx = 1:connParams.numAsters
    for jdx = 1:connParams.astralNum
        plot(connAsters.centers(idx,1) + connParams.l * cos(connAsters.orients(idx,jdx)) * (0:1), ...
            connAsters.centers(idx,2) + connParams.l * sin(connAsters.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',3,'Color','k')
    end
end
for idx = 1:size(connNetwork.springs,1)
    augNodeA = connNetwork.springs(idx,1);
    augNodeB = connNetwork.springs(idx,2);
    coords = [connNetwork.augNodes(augNodeA,1:2);
        connNetwork.augNodes(augNodeB,1:2)];
    plot(coords(:,1), coords(:,2),'.-','MarkerSize',18,'LineWidth',3, ...
        'Color',uciOrange)
    % if augNodeA <= connParams.numAsters && connParams.astralNum >= 2
    %     % i.e., if augNodeA is an astral center
    %     plot(connNetwork.augNodes(augNodeA,1),connNetwork.augNodes(augNodeA,2), ...
    %         '.r','MarkerSize',24)
    % elseif augNodeB <= connParams.numAsters && connParams.astralNum >= 2
    %     % likewise for augNodeB
    %     plot(connNetwork.augNodes(augNodeB,1),connNetwork.augNodes(augNodeB,2), ...
    %         '.r','MarkerSize',24)
    % end
end
if connParams.astralNum >= 2
    plot(connAsters.centers(:,1),connAsters.centers(:,2),'.r','MarkerSize',24)
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - connParams.D/2)) * [-1,1] + connParams.D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - connParams.D/2)) * [-1,1] + connParams.D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
% exportgraphics(connEx,fullfile(figureFolder,['connEx_rnd',filetype]), ...
%     'ContentType','vector')

%% PERCOLATION PROBABILITY HEATMAPS %%

percFsz = 16;
percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
chosenTypes = [1,5];
numDensVals = length(densityRange);
heatmaps = zeros(24,numDensVals,2);
networkLabel = @(an) sprintf('an%02i',an);
for idx = 1:24
    theseCurves = curves.(networkLabel(idx));
    for pdx = 1:2
        heatmaps(idx,:,pdx) = clip(spline(log10(theseCurves(1,:)), ...
            theseCurves(1+chosenTypes(pdx),:),log10(densityRange)),0,1);
    end
end

percMaps = figure(12);
set(percMaps,'units','centimeters','Position',[1,1,12,18])
percMaps = tiledlayout(percMaps,2,1,'TileSpacing','loose');
nexttile(1)
pcolor(densityRange,1:24,heatmaps(:,:,1))
xscale('log')
yticks([1,6,12,18,24])
title('Spanning percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
nexttile(2)
pcolor(densityRange,1:24,heatmaps(:,:,2))
xscale('log')
yticks([1,6,12,18,24])
title('Connectivity percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
xlabel(percMaps,'Density [$\mu{\rm m^{-1}}$]','Interpreter','latex', ...
    'FontSize',percFsz,'FontName',fname)
ylabel(percMaps,'Astral number $a_n$','Interpreter','latex','FontSize', ...
    percFsz,'FontName',fname)
bar = colorbar;
title(bar,'Probability')
bar.Layout.Tile = 'east';
% exportgraphics(percMaps,fullfile(figureFolder,['percMaps', ...
%     filetype]),'ContentType','vector')

%% PERCOLATION THRESHOLD ESTIMATES %%

percTypes = [1,5];
% Dlist = [5,10,15,20] so use D index 2
numPercTypes = length(percTypes);
% to10 = zeros(numTypes,24);
to50 = zeros(numPercTypes,24);
% to90 = zeros(numTypes,24);
opts = optimoptions("fsolve","Display","off");
for idx = 1:numPercTypes
    for jdx = 1:24
        thisCurve = smsplFits{jdx,percTypes(idx),2};
        % to10(idx,jdx) = fsolve(@(x) thisCurve(x) - 0.1, 5, opts);
        if percTypes(idx) == 5 && jdx ==1
            initGuess = 50;
        else
            initGuess = 10;
        end
        to50(idx,jdx) = fsolve(@(x) thisCurve(x) - 0.5, initGuess, opts);
        % to90(idx,jdx) = fsolve(@(x) thisCurve(x) - 0.9, 5, opts);
    end
end
% approxMiddle = (to10 + to90) / 2;

threshPlot = figure(13); clf;
set(threshPlot,'units','centimeters','Position',[1,1,15,12], ...
    'defaultLineLineWidth',2.5)
hold on
for idx = 1:numPercTypes
    plot(1:24,to50(idx,:),'Color',colors(idx,:),'LineStyle','-')
    % plot(1:24,approxMiddle(idx,:),'Color',colors(idx,:),'LineStyle','-.')
end
yline(6,':','6','Interpreter','latex','LineWidth',2, ...
    'FontSize',insetFsz,'LabelVerticalAlignment','bottom')
yline(7.5,'-','7.5','Interpreter','latex','LineWidth',2, ...
    'FontSize',insetFsz,'LabelVerticalAlignment','middle')
yline(9,'--','9','Interpreter','latex','LineWidth',2, ...
    'FontSize',insetFsz,'LabelVerticalAlignment','top')
% plot(1:24,7.5*ones(1,24),'--k')
hold off
xticks([1,4:4:24])
xlim([1,24])
xlabel('Astral number $a_n$')
yticks(2.^(2:5))
yticklabels(arrayfun(@(x) sprintf("$2^%i$",x),2:5))
ylim([4,51])
yscale('log')
ylabel('Critical density $\rho_c$ [$\mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz,'FontName',fname)
legend({'Spanning','Connectivity'},'Location', ...
    'northeast','FontSize',16,'FontName',fname)
% exportgraphics(threshPlot,fullfile(figureFolder,['threshPlot', ...
%     filetype]),'ContentType','Vector')

%% CORRELATION BETWEEN deltaRHO AND MODULI %%

connRhoc = transpose(to50(2,:));
simDensities = repmat(6:1.5:9,[24,1]);
deltaRho = simDensities - repmat(connRhoc,[1,3]);

moduli_loDens = squeeze(moduli_s_otherRhos(:,1,:));
moduli_hiDens = squeeze(moduli_s_otherRhos(:,2,:));

% compute a correlation coefficient
allModuli = [reshape(moduli_loDens,[],1); reshape(moduli_s,[],1);
    reshape(moduli_hiDens,[],1)];
deltaRho_vector = [repmat(deltaRho(:,1),otherRhoParams.numRep,1); 
    repmat(deltaRho(:,2),sweep_s.numRep,1); 
    repmat(deltaRho(:,3),otherRhoParams.numRep,1)];
filter = logical((allModuli>0) .* (deltaRho_vector<0));
correlationCoeff = -corrcoef(log10(allModuli(filter)), ...
    log10(-deltaRho_vector(filter)));

corrPlot = figure(14); clf;
set(corrPlot,'units','centimeters','Position',[1,1,18,12])
hold on
% first 3 plot calls facilitate legend format
plot(deltaRho(1,1),moduli_loDens(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho(1,2),moduli_s(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho(1,3),moduli_hiDens(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho(:,1),[1,15]),moduli_loDens,'Marker','v', ...
    'MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho(:,2),[1,30]),moduli_s,'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho(:,3),[1,15]),moduli_hiDens,'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-15,-0.65])
xscale('log')
yscale('log')
xlabel('Density rel. to critical density $\rho - \rho^{\rm conn}_c(a_n)$')
ylabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz2,'FontName',fname)
leglabels = {sprintf('%1.1f',otherRhoParams.dens(1)),sprintf('%1.1f', ...
    sweep_s.dens),sprintf('%1.1f',otherRhoParams.dens(2))};
lg = legend(leglabels,'Location','southeast','FontSize',insetFsz);
title(lg,'$\rho$ [${\rm \mu m^{-1}}$]')
% exportgraphics(corrPlot,fullfile(figureFolder,['corrPlot',filetype]), ...
%     'ContentType','vector')

%% SUPP: SHEAR MODULI VS. ASTRAL NUMBER AT 2 OTHER FILAMENT DENSITIES %%

[stdevMod_otherRhos,meanMod_otherRhos] = std(moduli_s_otherRhos,0,3);
markerList = {'v','^'};
otherRhos = figure(15);
set(otherRhos,'units','centimeters','Position',[1,1,18,12])
hold on
for densIdx = 1:otherRhoParams.numDens
    for netIdx = 1:otherRhoParams.nNetTypes
        plot(otherRhoParams.nFilPerAsterList(netIdx), ...
            squeeze(moduli_s_otherRhos(netIdx,densIdx,:)), ...
            'Marker',markerList{densIdx},'MarkerSize',3,'MarkerFaceColor', ...
            lgray,'MarkerEdgeColor',lgray)
    end
end
for idx = 1:sweep_s.nNetTypes
    plot(sweep_s.nFilPerAsterList(idx),moduli_s(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
l1 = errorbar(otherRhoParams.nFilPerAsterList,meanMod_otherRhos(:,1), ...
    z*stdevMod_otherRhos(:,1)/sqrt(otherRhoParams.numRep),'Color','k', ...
    'Marker',markerList{1},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams.dens(1)));
l2 = errorbar(otherRhoParams.nFilPerAsterList,meanMod_otherRhos(:,2), ...
    z*stdevMod_otherRhos(:,2)/sqrt(otherRhoParams.numRep),'Color', ...
    uciOrange,'Marker',markerList{2},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams.dens(2)));
l3 = errorbar(sweep_s.nFilPerAsterList, means_s, ...
    z*stdevs_s/sqrt(sweep_s.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue,'DisplayName',sprintf('%1.1f',sweep_s.dens));
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,otherRhoParams.nFilPerAsterList(end)])
xticks(0:4:otherRhoParams.nFilPerAsterList(end))
xlabel('Astral number $a_n$')
yl = ylim;
ylim([0,yl(2)])
ylabel({'Shear modulus $G$'; '[${\rm pN} \cdot \mu{\rm m}^{-1}$]'})
lg = legend([l1,l3,l2],'Location','northeast');
title(lg,'$\rho$ [${\rm \mu m^{-1}}$]')
% exportgraphics(otherRhos,fullfile(figureFolder,['otherRhos',filetype]), ...
%     'ContentType','vector')

%% SUPP: FILAMENT BENDING RIGIDITY (kbend) %%

[stdevMod_kbend,meanMod_kbend] = std(moduli_s_kbend,0,3);

kbendShear = figure(98);
set(kbendShear,'units','centimeters','Position',[1,1,15,12])
cmap = colormap(parula(kbendParams.nKbend+1));
hold on
for jdx = 1:kbendParams.nKbend
    errorbar(kbendParams.nFilPerAsterList,meanMod_kbend(:,jdx), ...
        (z/sqrt(kbendParams.numRep))*stdevMod_kbend(:,jdx),'Color', ...
        cmap(jdx,:),'LineWidth',0.75)
end
set(gca,'FontName',fname,'FontSize',fsz2)
lg = legend(string(kbendParams.kbend_list),'Location','northeast', ...
    'NumColumns',2);
lg.FontSize = insetFsz;
title(lg,{"Filament bending"; "rigidity [${\rm pN}\cdot \mu {\rm m}^2$]"})
xlabel('Astral number $a_n$')
xticks(kbendParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN} \cdot \mu{\rm m^{-1}}$]')
yl = ylim;
ylim([0,yl(2)])

% exportgraphics(kbendShear,fullfile(figureFolder,['kbendShear',filetype]), ...
%     'ContentType','vector')

%% SUPP: ANGULAR STIFFNESS AT ASTRAL CENTERS (krot) %%

[stdevMod_krot,meanMod_krot] = std(moduli_s_krot,0,3);

krotShear = figure(97);
set(krotShear,'units','centimeters','Position',[1,1,15,12])
cmap = colormap(parula(krotParams.nKrot+1));
hold on
for jdx = 1:krotParams.nKrot
    errorbar(krotParams.nFilPerAsterList,meanMod_krot(:,jdx), ...
        (z/sqrt(krotParams.numRep))*stdevMod_krot(:,jdx),'Color', ...
        cmap(jdx,:),'LineWidth',0.75)
end
set(gca,'FontName',fname,'FontSize',fsz2)
lg = legend(string(krotParams.krot_list),'Location','northeast', ...
    'NumColumns',2);
lg.FontSize = insetFsz;
title(lg,{"Angular stiffness at"; "centers [${\rm pN}\cdot {\mu}{\rm m^{-1}}$]"})
xlabel('Astral number $a_n$')
xticks(krotParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
yl = ylim;
ylim([0,yl(2)])

% exportgraphics(krotShear,fullfile(figureFolder,['krotShear',filetype]), ...
%     'ContentType','vector')

%% SUPP: CROSSLINKER DENSITY %%

[stdevMod_crslnkDens,meanMod_crslnkDens] = std(moduli_s_crslnkDens,0,3);

crslnkDensShear = figure(96);
set(crslnkDensShear,'units','centimeters','Position',[1,1,15,15])
cmap = colormap(parula(crslnkDensParams.nCrslnkDens+1));
hold on
for jdx = 1:crslnkDensParams.nCrslnkDens
    errorbar(crslnkDensParams.nFilPerAsterList,meanMod_crslnkDens(:,jdx), ...
        (z/sqrt(crslnkDensParams.numRep))*stdevMod_crslnkDens(:,jdx),'Color', ...
        cmap(jdx,:),'LineWidth',0.75)
end
lg = legend(string(crslnkDensParams.crosslinker_dens_list), ...
    'Location','northeast');
title(lg,{"Crosslinker density"; "[${\rm particles}\cdot\mu {\rm m}^{-2}$]"})
xlabel('Astral number $a_n$')
xticks(crslnkDensParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
yl = ylim;
ylim([0,yl(2)])
set(gca,'FontName',fname,'FontSize',fsz2)
% exportgraphics(crslnkDensShear,fullfile(figureFolder, ...
%     ['crslnkDensShear',filetype]),'ContentType','vector')

%% SUPP: CROSSLINKER STIFFNESS %%

[stdevMod_crslnkStiff,meanMod_crslnkStiff] = std(moduli_s_crslnkStiff,0,3);

crslnkStiffShear = figure(95);
set(crslnkStiffShear,'units','centimeters','Position',[1,1,15,15])
cmap = colormap(parula(crslnkStiffParams.nCrslnkStiff+1));
hold on
for jdx = 1:crslnkStiffParams.nCrslnkStiff
    errorbar(crslnkStiffParams.nFilPerAsterList,meanMod_crslnkStiff(:,jdx), ...
        (z/sqrt(crslnkStiffParams.numRep))*stdevMod_crslnkStiff(:,jdx),'Color', ...
        cmap(jdx,:),'LineWidth',0.75)
end
lg = legend(string(crslnkStiffParams.crosslinker_stiff_list), ...
    'Location','northeast');
title(lg,{"Crosslinker stiffness"; "[${\rm pN}\cdot \mu {\rm m}^{-1}$]"})
xlabel('Astral number $a_n$')
xticks(crslnkStiffParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot\mu{\rm m}^{-1}$]')
yl = ylim;
ylim([0,yl(2)])
set(gca,'FontName',fname,'FontSize',fsz2)
% exportgraphics(crslnkStiffShear,fullfile(figureFolder, ...
%     ['crslnkStiffShear',filetype]),'ContentType','vector')

%% SUPP: MIKADO RIGIDITY (Wilhelm & Frey 2003) %%

[stdevMod_mikado,meanMod_mikado] = std(moduli_s_mikado,0,3);

nu = 4/3;
muOverNu = 2.3;
rho_c = 6.68; % W&F's reported rigidity percolation thresh for free hinges
deltaRho = mikadoParams.dens_list - rho_c;

xCoords = log10(repmat(mikadoParams.D_list,[1,mikadoParams.numDens]) .* ...
    (repmat(deltaRho,[mikadoParams.numD,1]).^nu));
yCoords = log10(meanMod_mikado .* ...
    repmat(mikadoParams.D_list,[1,mikadoParams.numDens]).^muOverNu);

linespecs = {'-ok','-^m','-vb','-<r','->c','-+g','-diamondm','-xr','-squareb'};

mikadoRig = figure(94);
set(mikadoRig,'units','centimeters','Position',[1,1,15,15])
hold on;
for jdx = 1:mikadoParams.numDens
    plot(xCoords(:,jdx),yCoords(:,jdx),linespecs{jdx},'LineWidth',0.7)
end
hold off;
legend(string(mikadoParams.dens_list),'Location','eastoutside')
xlabel('$\log_{10} (D(\rho-\rho_c)^\nu)$','Interpreter','latex')
ylabel('$\log_{10} (D^{\mu/\nu}\cdot G)$','Interpreter','latex')

% exportgraphics(mikadoRig,fullfile(figureFolder,['mikadoRig',filetype]), ...
%     "ContentType",'vector');

%% SUPP: TIMESERIES %%

% shear
timeserShear = figure(93); clf;
set(timeserShear,'units','centimeters','Position',[1,1,20,8], ...
    'defaultLineLineWidth',1.5)
timeserShear = tiledlayout(timeserShear,1,2);
nexttile(1)
hold on
xline(5,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams.shearForce/timeserParams.D))
xline(45,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,1), ...
        'Color',colors(netIdx,:))
end
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]')
set(gca,'FontSize',fsz2,'FontName',fname)

leglabels = {'',''};
leglabels = cat(2,leglabels,cellstr(string(timeserParams.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',insetFsz);
title(lg,'$a_n$','Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(5,'--k','LineWidth',1.5)
xline(45,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,2), ...
        'Color',colors(netIdx,:))
end
ylim([-1,1])
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]')
set(gca,'FontSize',fsz2,'FontName',fname)

% exportgraphics(timeserShear,fullfile(figureFolder, ...
%     ['timeserShear',filetype]),'ContentType','vector')

% tensile
timeserTens = figure(92); clf;
set(timeserTens,'units','centimeters','Position',[1,1,20,8], ...
    'defaultLineLineWidth',1.5)
timeserTens = tiledlayout(timeserTens,1,2);
nexttile(1)
hold on
xline(5,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams.shearForce/timeserParams.D))
xline(45,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,1), ...
        'Color',colors(netIdx,:))
end
hold off
xlim('tight')
ylim([-0.5,0.5])
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]')
set(gca,'FontSize',fsz2,'FontName',fname)

leglabels = {'',''};
leglabels = cat(2,leglabels,cellstr(string(timeserParams.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',insetFsz);
title(lg,'$a_n$','Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(5,'--k','LineWidth',1.5)
xline(45,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,2), ...
        'Color',colors(netIdx,:))
end
ylim([-0.5,0.5])
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]')
set(gca,'FontSize',fsz2,'FontName',fname)

% exportgraphics(timeserTens,fullfile(figureFolder, ...
%     ['timeserTens',filetype]),'ContentType','vector')

%% SUPP: DISPLACEMENT VS FORCE CURVES %%

selected_an = [1,4,8,16];
nNetTypes = length(selected_an);
networkLabel = @(an) sprintf('an%02i',an);

% shear
forceDispShear = figure(91); clf;
set(forceDispShear,'units','centimeters','Position',[1,1,20,20], ...
    'defaultLineLineWidth',1.5)
forceDispShear = tiledlayout(forceDispShear,2,2,'TileSpacing','compact');
for netIdx = 1:nNetTypes
    theseDisps = disps_s.(networkLabel(selected_an(netIdx)));
    nexttile(netIdx)
    hold on
    for repIdx = 1:sweep_s.numRep
        plot(sweep_s.forceVals,theseDisps(:,1,repIdx),'-xk','LineWidth', ...
            0.5,'Color',colors(netIdx,:))
    end
    hold off
    xlim('tight')
    ylim('tight')
    title(sprintf("$a_n = %i$",selected_an(netIdx)),'Interpreter','latex')
    set(gca,'FontSize',fsz2,'FontName',fname)
    xticks(sweep_s.forceVals)
end
xlabel(forceDispShear,'Shear force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispShear,'$x$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

% exportgraphics(forceDispShear,fullfile(figureFolder, ...
%     ['forceDispShear',filetype]),'ContentType','vector')

% tensile
forceDispTens = figure(90); clf;
set(forceDispTens,'units','centimeters','Position',[1,1,20,20], ...
    'defaultLineLineWidth',1.5)
forceDispTens = tiledlayout(forceDispTens,2,2);
for netIdx = 1:nNetTypes
    theseDisps = disps_t.(networkLabel(selected_an(netIdx)));
    nexttile(netIdx)
    hold on
    for repIdx = 1:sweep_t.numRep
        plot(sweep_t.forceVals,theseDisps(:,2,repIdx),'-xk','LineWidth', ...
            0.5,'Color',colors(netIdx,:))
    end
    hold off
    xlim('tight')
    ylim('tight')
    title(sprintf("$a_n = %i$",selected_an(netIdx)),'Interpreter','latex')
    set(gca,'FontSize',fsz2,'FontName',fname)
    xticks(sweep_t.forceVals)
end
xlabel(forceDispTens,'Tensile force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispTens,'$y$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

% exportgraphics(forceDispTens,fullfile(figureFolder, ...
%     ['forceDispTens',filetype]),'ContentType','vector')

%% SUPP: Percolation probability curves + fits

an_to_show = [1,2,4,8];
nNetTypes = length(an_to_show);
nSizes = length(Dlist);

spanDataFits = figure(89); clf;
set(spanDataFits,'units','centimeters','Position',[1,1,25,10])
spanDataFits = tiledlayout(spanDataFits,1,4);
legendLabels = cell(2*nSizes,1);
for idx = 1:nNetTypes
    nexttile(idx)
    hold on
    for kdx = 1:nSizes
        thisData = percDataUsed{an_to_show(idx),1,kdx};
        l = plot(smsplFits{an_to_show(idx),1,kdx}, ...
            thisData(:,1),thisData(:,2),'*');
        [l.Color] = deal(colors(kdx,:));
        [l.LineWidth] = deal(0.75);
        legendLabels{2*kdx - 1} = ''; % suppress data entry in legend
        legendLabels{2*kdx} = sprintf('%2.1f',Dlist(kdx));
    end
    yline(0.5,'--','LineWidth',1)
    hold off
    title(sprintf('$a_n = %i$',an_to_show(idx)))
    xlabel('')
    ylabel('')
    xscale('log')
    set(gca,'FontSize',fsz2,'FontName',fname)
    if idx == nNetTypes
        lg = legend(legendLabels,'FontSize',insetFsz);
        lg.Layout.Tile = 'east';
        title(lg,"$s \; [{\rm \mu m}]$",'FontSize',insetFsz)
        xlabel(spanDataFits,'Density $\rho$ $[{\rm \mu m^{-1}}]$', ...
            'Interpreter','latex','FontSize',fsz2)
        ylabel(spanDataFits,'$\Pr\{$Spanning comp.$\}$', ...
            'Interpreter','latex','FontSize',fsz2)
    else
        legend('off')
    end
end
% exportgraphics(spanDataFits,fullfile(figureFolder,['spanDataFits', ...
%     filetype]),'ContentType','vector')

connDataFits = figure(88); clf;
set(connDataFits,'units','centimeters','Position',[1,1,25,10])
connDataFits = tiledlayout(connDataFits,1,4);
legendLabels = cell(2*nSizes,1);
for idx = 1:nNetTypes
    nexttile(idx)
    hold on
    for kdx = 1:nSizes
        thisData = percDataUsed{an_to_show(idx),5,kdx};
        l = plot(smsplFits{an_to_show(idx),5,kdx}, ...
            thisData(:,1),thisData(:,2),'*');
        [l.Color] = deal(colors(kdx,:));
        [l.LineWidth] = deal(0.75);
        legendLabels{2*kdx - 1} = ''; % suppress data entry in legend
        legendLabels{2*kdx} = sprintf('%2.1f',Dlist(kdx));
    end
    yline(0.5,'--','LineWidth',1)
    hold off
    title(sprintf('$a_n = %i$',an_to_show(idx)))
    xlabel('')
    ylabel('')
    xscale('log')
    set(gca,'FontSize',fsz2,'FontName',fname)
    if idx == nNetTypes
        lg = legend(legendLabels,'FontSize',insetFsz);
        lg.Layout.Tile = 'east';
        title(lg,"$s \; [{\rm \mu m}]$",'FontSize',insetFsz)
        xlabel(connDataFits,'Density $\rho$ $[{\rm \mu m^{-1}}]$', ...
            'Interpreter','latex','FontSize',fsz2)
        ylabel(connDataFits,'$\Pr\{$Unique conn. comp.$\}$', ...
            'Interpreter','latex','FontSize',fsz2)
    else
        legend('off')
    end
end
% exportgraphics(connDataFits,fullfile(figureFolder,['connDataFits', ...
%     filetype]),'ContentType','vector')

%% SUPP: PRODUCTIVE SEGMENTS PER NODE

numAugNodes_byNet = cellfun(@length,springPerNodeCounts);
[augNodes_stdev,augNodes_means] = std(numAugNodes_byNet,0,2);
augNodes_ci_95 = augNodes_stdev * z / sqrt(springCountParams.N);

springPerNodeMeansPlot = figure(87); clf;
set(springPerNodeMeansPlot,'units','centimeters','Position',[1,1,10,16])
springPerNodeMeansPlot = tiledlayout(2,1);
nexttile(1)
errorbar(springCountParams.an_list,springPerNode_means_by_net_type, ...
    springPerNode_ci_95,"vertical",'LineWidth',1.5)
xticks([1,4:4:24])
ylim([3,4])
ylabel({"Productive"; "segments per node"})
set(gca,"FontName",fname,"FontSize",fsz2)
nexttile(2)
errorbar(springCountParams.an_list,augNodes_means,augNodes_ci_95, ...
    "vertical",'LineWidth',1.5)
xticks([1,4:4:24])
ylabel("Number of nodes")
set(gca,"FontName",fname,"FontSize",fsz2)
xlabel(springPerNodeMeansPlot,"Astral number $a_n$",'FontSize',fsz2, ...
    'Interpreter','latex')
% exportgraphics(springPerNodeMeansPlot,fullfile(figureFolder, ...
%     ['springsPerNodeMeans',filetype]),'ContentType','vector')

springPerNodeHists = figure(86); clf;
set(springPerNodeHists,'units','centimeters','Position',[1,1,24,16])
springPerNodeHists = tiledlayout(2,3);
for idx = 1:springCountParams.num_selected
    nexttile(idx)
    histogram(categorical(aggregatedCounts{idx}),'Normalization','count')
    title(sprintf("$a_n = %i$",springCountParams.selected_an(idx)))
    yscale('log')
    ylim([10^2,10^6+5e5])
    yticks(10.^(2:2:6))
    set(gca,'FontSize',insetFsz)
    if springCountParams.selected_an(idx) == 16
        labels = cellstr(string(0:16));
        [labels{2:2:16}] = deal('');
        ax = gca;
        ax.XAxis.TickLabels = labels;
    elseif springCountParams.selected_an(idx) == 24
        labels = cellstr(string(0:24));
        [labels{2:2:24}] = deal('');
        [labels{3:4:24}] = deal('');
        ax = gca;
        ax.XAxis.TickLabels = labels;
    end
end
xlabel(springPerNodeHists,'Productive segments per node','Interpreter', ...
    'latex','FontSize',fsz2)
ylabel(springPerNodeHists,sprintf("Aggregate counts ($N = %i$)", ...
    springCountParams.N),'Interpreter','latex','FontSize',fsz2)
% exportgraphics(springPerNodeHists,fullfile(figureFolder, ...
%     ['springPerNodeHists',filetype]),'ContentType','vector')

%% FIG1: SAMPLE DISPLACEMENT VS FORCE CURVE

sampleDispForceShear = figure(80); clf;
set(sampleDispForceShear,'defaultLineLineWidth',1.5,'defaultLineMarkerSize',20)
repIdx = 26;
hold on
plot(sweep_s.forceVals,disps_s.an01(:,1,repIdx),'-*b')
plot(sweep_s.forceVals,sweep_s.forceVals / moduli_s(1,repIdx),'-k')
plot([1.5,2.5,2.5],[1.5,1.5,2.5]./moduli_s(1,repIdx),':k','LineWidth',3)
text(2.75,0.133,'$G \propto 1/{\rm slope}$','FontSize',fsz)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xticks(sweep_s.forceVals)
xlim('tight')
xlabel('Shear force [pN]')
ylabel('$x$-Displacement [${\rm \mu m}$]')

% exportgraphics(sampleDispForceShear,fullfile(figureFolder, ...
%     ['sampleDispForceShear',filetype]),'ContentType','vector')

sampleDispForceTens = figure(79); clf;
set(sampleDispForceTens,'defaultLineLineWidth',1.5,'defaultLineMarkerSize',20)
repIdx = 26;
hold on
plot(sweep_t.forceVals,disps_t.an01(:,2,repIdx),'-*b')
plot(sweep_t.forceVals,sweep_t.forceVals / moduli_t(1,repIdx),'-k')
plot([6,10,10],[6,6,10]./moduli_t(1,repIdx),':k','LineWidth',3)
text(11,0.1,'$Y \propto 1/{\rm slope}$','FontSize',fsz)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xticks(sweep_t.forceVals)
xlim('tight')
xlabel('Tensile force [pN]')
ylabel('$y$-Displacement [${\rm \mu m}$]')

% exportgraphics(sampleDispForceTens,fullfile(figureFolder, ...
%     ['sampleDispForceTens',filetype]),'ContentType','vector')