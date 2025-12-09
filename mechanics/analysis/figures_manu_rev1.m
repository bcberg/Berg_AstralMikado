% figures_manu_rev1.m
% Brady Berg, 10/29/2025
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Importing data

% Ubuntu paths
saveDir = '~/Documents/astral-mikado-data';
% astralFuncFolder = '~/Documents/astral-mikado/';
% Windows paths
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';
% astralFuncFolder = 'C:\Users\bcber\Documents\astral-mikado\';
dataSubdir = 'mat_files';
figSubdir = 'subfigures';

filetype = '.pdf';

% shear modulus sweep - random angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'rep_shear_actin251015_l0.1_D1.mat'), ...
    'nFilPerAsterList','shearModuli','numRep','nNetTypes', ...
    'disps','forceVals','dens')
sweep_s_rnd = struct('nFilPerAsterList',nFilPerAsterList,'numRep',numRep, ...
    'nNetTypes',nNetTypes,'forceVals',forceVals,'dens',dens,'len_fil', ...
    0.1,'D',1);
disps_s_rnd = disps;
moduli_s_rnd = shearModuli;
clear('shearModuli','disps','nFilPerAsterList','numRep','nNetTypes', ...
    'forceVals','dens')

% Young's (tensile) modulus sweep - random angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'rep_tens_actin251016_l0.1_D1.mat'), ...
    'nFilPerAsterList','tensModuli','numRep','nNetTypes', ...
    'disps','forceVals','dens')
sweep_t = struct('nFilPerAsterList',nFilPerAsterList,'numRep',numRep, ...
    'nNetTypes',nNetTypes,'forceVals',forceVals,'dens',dens,'len_fil', ...
    0.1,'D',1);
disps_t = disps;
moduli_t = tensModuli;
clear('tensModuli','disps','nFilPerAsterList','numRep','nNetTypes', ...
    'forceVals','dens')

% modulus distribution sampling - random angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'mod_distr251020.mat'),'D','len_fil', ...
    'nFilPerAsterList','dens','nNetTypes','numRep','shearForces', ...
    'shearModuli','shearLineFits','percStatus_s','tensForces', ...
    'tensModuli','tensLineFits','percStatus_t')
moduli_s_hiN = shearModuli;
shearLineFits_hiN = shearLineFits;
percStatus_s_hiN = percStatus_s;
moduli_t_hiN = tensModuli;
tensLineFits_hiN = tensLineFits;
percStatus_t_hiN = percStatus_t;
modSampling = struct('D',D,'len_fil',len_fil,'nFilPerAsterList', ...
    nFilPerAsterList,'dens',dens,'nNetTypes',nNetTypes,'numRep',numRep, ...
    'shearForces',shearForces,'tensForces',tensForces);
clear('shearModuli','shearLineFits','percStatus_s','tensModuli', ...
    'tensLineFits','percStatus_t','D','len_fil','nFilPerAsterList', ...
    'dens','nNetTypes','numRep','shearForces','tensForces')

% dangling ends data - random angle asters
% note that l = 1 and D = 10 here
load(fullfile(saveDir,dataSubdir,'danglingEnds_250319.mat'),'numNets', ...
        'nFilPerAsterList','numNetTypes','endsLengths','endsMeans', ...
        'stdOfMeans','endsMeansTrunc','stdOfMeansTrunc','rho','l','D', ...
        'targetFilNum','labelPat','useFraction')
danglingParams = struct('nFilPerAsterList',nFilPerAsterList,'l',l, ...
    'D',D,'rho',rho,'numNets',numNets,'numNetTypes',numNetTypes, ...
    'targetFilNum',targetFilNum,'labelPat',labelPat);
clear('nFilPerAsterList','l','D','rho','numNets','numNetTypes', ...
    'targetFilNum','labelPat')

% example astral networks in MATLAB - random angle asters
% (unchanged from original submission)

% percolation data - random angle asters
load(fullfile(saveDir,dataSubdir,sprintf("percProbs_l%02i_D%02i",1,10)), ...
    'curves','densityRange')
curves_rnd = curves;
clear('curves')

% smoothing spline fits to percolation data - random angle asters
% includes estimates of critical densities (at l = 1, D = 10)
load(fullfile(saveDir,dataSubdir,'percFits_nonparam_tzOnly.mat'), ...
    'smsplFits','l','D','percDataUsed','percTypeChoice', ...
    'numPercTypesChosen','numNetTypes','astralNumList','to50')
smsplFits_rnd = smsplFits;
percDataUsed_rnd = percDataUsed;
to50_rnd = to50;
percParams_rnd = struct('densityRange',densityRange,'l',l,'D',D, ...
    'percTypeChoice',percTypeChoice,'numPercTypesChosen',numPercTypesChosen, ...
    'numNetTypes',numNetTypes,'astralNumList',astralNumList);
clear('curves','smsplFits','percDataUsed','to50','densityRange','l','D', ...
    'percTypeChoice','numPercTypesChosen','numNetTypes','astralNumList')

% smoothing spline fits to percolation data - random angle asters
% family of system sizes for l = 1
load(fullfile(saveDir,dataSubdir,'percMCRG_nonparamFits_tzOnly.mat'), ...
    'Dlist','percDataUsed','numNetTypes','numPercTypes','numD', ...
    'allAstralNums','allCurves')
% smsplFits_rnd_famOfD = smsplFits;     % need to redo fits
percDataUsed_rnd_famOfD = percDataUsed;
percParams_rnd_famOfD = struct('Dlist',Dlist,'l',1,'numD',numD, ...
    'numNetTypes',numNetTypes,'numPercTypes',numPercTypes,'allAstralNums', ...
    allAstralNums);
clear('Dlist','percDataUsed','numNetTypes','numPercTypes','numD', ...
    'allAstralNums')

% example astral networks in MATLAB - equal angle asters
% spanning example
load(fullfile(saveDir,dataSubdir,'spanEx_eq.mat'),'astralNum','rho', ...
    'l','D','targetFilNum','numAsters','net_eq_span', ...
    'crossings_eq_span','asters_eq_span')
tbSpanParams_eq = struct('astralNum',astralNum,'rho',rho,'l',l,'D',D, ...
    'targetFilNum',targetFilNum,'numAsters',numAsters);
clear('astralNum','rho','l','D','targetFilNum','numAsters')
% single connected component example
load(fullfile(saveDir,dataSubdir,'connEx_eq.mat'),'astralNum','rho', ...
    'l','D','targetFilNum','numAsters','net_eq_conn', ...
    'crossings_eq_conn','asters_eq_conn')
connParams_eq = struct('astralNum',astralNum,'rho',rho,'l',l,'D',D, ...
    'targetFilNum',targetFilNum,'numAsters',numAsters);
clear('astralNum','rho','l','D','targetFilNum','numAsters')

% percolation data - equal angle asters
load(fullfile(saveDir,dataSubdir, ...
    sprintf("percProbs_eqSpaced_l%02i_D%02i",1,10)),'curves','densityRange')
curves_eq = curves;
clear('curves')

% smoothing spline fits to percolation data - equal angle asters
% includes estimates of critical densities (at l = 1, D = 10)
load(fullfile(saveDir,dataSubdir,'percFits_nonparam_tzOnly_eqSpaced.mat'), ...
    'smsplFits','l','D','percDataUsed','percTypeChoice', ...
    'numPercTypesChosen','numNetTypes','astralNumList','to50')
smsplFits_eq = smsplFits;
percDataUsed_eq = percDataUsed;
to50_eq = to50;
percParams_eq = struct('densityRange',densityRange,'l',l,'D',D, ...
    'percTypeChoice',percTypeChoice,'numPercTypesChosen',numPercTypesChosen, ...
    'numNetTypes',numNetTypes,'astralNumList',astralNumList);
clear('curves','smsplFits','percDataUsed','to50','densityRange','l','D', ...
    'percTypeChoice','numPercTypesChosen','numNetTypes','astralNumList')

% shear modulus vs. astral number at additional filament densities
% random angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'rep_shear_otherRhos251021.mat'), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','forceVals', ...
    'numRep','nNetTypes','numDens','shearModuli')
moduli_s_otherRhos_rnd = shearModuli;
otherRhoParams_rnd = struct('dens',dens,'len_fil',len_fil,'D',D,'nFil', ...
    nFil,'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals, ...
    'numRep',numRep,'nNetTypes',nNetTypes,'numDens',numDens);
clear('dens','len_fil','D','nFil','nFilPerAsterList','forceVals', ...
    'numRep','nNetTypes','numDens','shearModuli')

% shear modulus vs. astral number - equal angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'rf_shear_eq_actin251111.mat'), ...
    'D','len_fil','nFilPerAsterList','shearModuli','numRep','nNetTypes',...
    'forceVals','dens')
sweep_s_eq = struct('nFilPerAsterList',nFilPerAsterList,'numRep',numRep, ...
    'nNetTypes',nNetTypes,'forceVals',forceVals,'dens',dens,'len_fil', ...
    len_fil,'D',D);
moduli_s_eq = shearModuli;
clear('shearModuli','nFilPerAsterList','numRep','nNetTypes', ...
    'forceVals','dens')

% shear modulus vs. astral number at additional filament densities
% equal angle asters
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,'rf_shear_aster_eq_otherRhos251115.mat'), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','forceVals', ...
    'numRep','nNetTypes','numDens','shearModuli')
moduli_s_otherRhos_eq = shearModuli;
otherRhoParams_eq = struct('dens',dens,'len_fil',len_fil,'D',D,'nFil', ...
    nFil,'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals, ...
    'numRep',numRep,'nNetTypes',nNetTypes,'numDens',numDens);
clear('dens','len_fil','D','nFil','nFilPerAsterList','forceVals', ...
    'numRep','nNetTypes','numDens','shearModuli')

% shear modulus vs. astral number - k_B*T at physiologic value
load(fullfile(saveDir,dataSubdir,'rf_shear_aster_actin_physkT251113.mat'), ...
    'D','len_fil','nFilPerAsterList','shearModuli','numRep',...
    'nNetTypes','forceVals','disps','dens','filtering')
moduli_s_physkT = shearModuli;
disps_s_physkT = disps;
sweep_s_physkT = struct('D',D,'len_fil',len_fil,'nFilPerAsterList', ...
    nFilPerAsterList,'numRep',numRep,'nNetTypes',nNetTypes,'forceVals', ...
    forceVals,'dens',dens,'filtering',filtering);
clear('D','len_fil','nFilPerAsterList','shearModuli','numRep',...
    'nNetTypes','forceVals','disps','dens','filtering')

% Young's modulus vs. astral number - k_B*T at physiologic value
load(fullfile(saveDir,dataSubdir,'rf_tens_aster_actin_physkT251114.mat'), ...
    'D','len_fil','nFilPerAsterList','tensModuli','numRep',...
    'nNetTypes','forceVals','disps','dens','filtering')
moduli_t_physkT = tensModuli;
disps_t_physkT = disps;
sweep_t_physkT = struct('D',D,'len_fil',len_fil,'nFilPerAsterList', ...
    nFilPerAsterList,'numRep',numRep,'nNetTypes',nNetTypes,'forceVals', ...
    forceVals,'dens',dens,'filtering',filtering);
clear('D','len_fil','nFilPerAsterList','tensModuli','numRep',...
    'nNetTypes','forceVals','disps','dens','filtering')

% crosslinker density sweep
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,"crosslinker_dens_shear251126.mat"), ...
    'dens','len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nNetTypes','crosslinker_per_fil_list','crosslinker_dens_list', ...
    'nCrslnkDens','shearModuli','filtering')
crslnkDensParams = struct('dens',dens,'len_fil',len_fil,'D',D, ...
    'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals, ...
    'numRep',numRep,'nNetTypes',nNetTypes,'crosslinker_per_fil_list', ...
    crosslinker_per_fil_list,'crosslinker_dens_list',crosslinker_dens_list, ...
    'nCrslnkDens',nCrslnkDens,'filtering',filtering);
moduli_s_crslnkDens = shearModuli;
clear('dens','len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nNetTypes','crosslinker_per_fil_list','crosslinker_dens_list', ...
    'nCrslnkDens','shearModuli','filtering')

% crosslinker stiffness sweep
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,"crosslinker_stiff_shear251127.mat"), ...
    'dens','len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nNetTypes','crosslinker_stiff_list','nCrslnkStiff','shearModuli', ...
    'filtering')
crslnkStiffParams = struct('dens',dens,'len_fil',len_fil,'D',D, ...
    'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals, ...
    'numRep',numRep,'nNetTypes',nNetTypes,'crosslinker_stiff_list', ...
    crosslinker_stiff_list,'nCrslnkStiff',nCrslnkStiff,'filtering', ...
    filtering);
moduli_s_crslnkStiff = shearModuli;
clear('dens','len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nNetTypes','crosslinker_stiff_list','nCrslnkStiff','shearModuli', ...
    'filtering')

% filament bending rigidity sweep
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,"kbend_shear251129.mat"),'dens', ...
    'len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nKbend','nNetTypes','kbend_list','shearModuli','filtering')
kbendParams = struct('dens',dens,'len_fil',len_fil,'D',D, ...
    'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals, ...
    'numRep',numRep,'nKbend',nKbend,'nNetTypes',nNetTypes, ...
    'kbend_list',kbend_list,'filtering',filtering);
moduli_s_kbend = shearModuli;
clear('dens','len_fil','D','nFilPerAsterList','forceVals','numRep', ...
    'nKbend','nNetTypes','kbend_list','shearModuli','filtering')

% angular stiffness at astral centers sweep
% uses actin-like params with k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,"krot_shear251202.mat"),'dens', ...
    'len_fil','D','nFilPerAsterList','forceVals','numRep','nKrot', ...
    'nNetTypes','krot_list','shearModuli','filtering')
krotParams = struct('dens',dens,'len_fil',len_fil,'D',D, ...
    'nFilPerAsterList',nFilPerAsterList,'forceVals',forceVals,'numRep', ...
    numRep,'nKrot',nKrot,'nNetTypes',nNetTypes,'krot_list',krot_list, ...
    'filtering',filtering);
moduli_s_krot = shearModuli;
clear('dens','len_fil','D','nFilPerAsterList','forceVals','numRep','nKrot', ...
    'nNetTypes','krot_list','shearModuli','filtering')

% timeseries at various astral numbers - k_B*T 1/200 physiologic
load(fullfile(saveDir,dataSubdir,"timeseries_astralNum251014.mat"), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce', ...
    'shearCoMs','tensCoMs','netLabelPat')
shearCoMs_lokT = shearCoMs;
tensCoMs_lokT = tensCoMs;
timeserParams = struct('dens',dens,'len_fil',len_fil,'D',D,'nFil',nFil, ...
    'nFilPerAsterList',nFilPerAsterList,'nNetTypes',nNetTypes, ...
    'frameVals',frameVals,'secsPerFrame',secsPerFrame,'shearForce', ...
    shearForce,'tensForce',tensForce,'netLabelPat',netLabelPat);
clear('dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce','netLabelPat', ...
    'shearCoMs','tensCoMs')

% timeseries at various astral numbers - k_B*T at physiologic value
load(fullfile(saveDir,dataSubdir,"timeseries_astralNum_physkT251113.mat"), ...
    'dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce', ...
    'shearCoMs','tensCoMs','netLabelPat')
shearCoMs_physkT = shearCoMs;
tensCoMs_physkT = tensCoMs;
timeserParams_physkT = struct('dens',dens,'len_fil',len_fil,'D',D,'nFil', ...
    nFil,'nFilPerAsterList',nFilPerAsterList,'nNetTypes',nNetTypes, ...
    'frameVals',frameVals,'secsPerFrame',secsPerFrame,'shearForce', ...
    shearForce,'tensForce',tensForce,'netLabelPat',netLabelPat);
clear('dens','len_fil','D','nFil','nFilPerAsterList','nNetTypes', ...
    'frameVals','secsPerFrame','shearForce','tensForce','netLabelPat', ...
    'shearCoMs','tensCoMs')

% counting productive segments per node
% (unchanged from original submission)

% percolation in Cytosim
% for later

%% General figure settings

%%%%% font specifications %%%%%
fsz = 22;
fsz2 = 18;
fszInset = 16;
fname = 'CMU Serif';

%%%%% standard colors %%%%%
colors = colororder("gem12");

%%%%% for 95% CIs %%%%%
z = norminv(0.975);

%%%%% common struct field %%%%%
networkLabel = @(an) sprintf('an%02i',an);

%%%%% CUSTOM COLORS %%%%%
lgray = [197,190,181]/255;
gray = [85,87,89]/255;
uciBlue = [0,100,164]/255;
uciGold = [255,210,0]/255;
uciOrange = [247,141,45]/255;

%% SAMPLE DISPLACEMENT VS. FORCE CURVE %%

sampleDispForceShear = figure(1); clf;
set(sampleDispForceShear,'defaultLineLineWidth',1.5, ...
    'defaultLineMarkerSize',20,'units','centimeters','Position',[1,1,20,16])
repIdx = 7;
hold on
plot(sweep_s_rnd.forceVals,disps_s_rnd.an01(:,1,repIdx),'-*b')
plot(sweep_s_rnd.forceVals,sweep_s_rnd.forceVals / moduli_s_rnd(1,repIdx),'-k')
plot([0.1,0.15,0.15],[0.1,0.1,0.15]./moduli_s_rnd(1,repIdx),':k','LineWidth',3)
text(0.16,0.02,'$G \propto 1/{\rm slope}$','FontSize',fsz)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xticks(sweep_s_rnd.forceVals)
xlim('tight')
ylim('tight')
xlabel('Shear force [pN]')
ylabel('$x$-Displacement [${\rm \mu m}$]')

exportgraphics(sampleDispForceShear,fullfile(saveDir,figSubdir, ...
    ['sampleDispForceShear',filetype]),'ContentType','vector')

sampleDispForceTens = figure(2); clf;
set(sampleDispForceTens,'defaultLineLineWidth',1.5, ...
    'defaultLineMarkerSize',20,'units','centimeters','Position',[1,1,20,16])
repIdx = 19;
hold on
plot(sweep_t.forceVals,disps_t.an01(:,2,repIdx),'-*b')
plot(sweep_t.forceVals,sweep_t.forceVals / moduli_t(1,repIdx),'-k')
plot([0.4,0.6,0.6],[0.4,0.4,0.6]./moduli_t(1,repIdx),':k','LineWidth',3)
text(0.635,0.0145,'$Y \propto 1/{\rm slope}$','FontSize',fsz)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xticks(sweep_t.forceVals)
xlim('tight')
xlabel('Tensile force [pN]')
ylabel('$y$-Displacement [${\rm \mu m}$]')

exportgraphics(sampleDispForceTens,fullfile(saveDir,figSubdir, ...
    ['sampleDispForceTens',filetype]),'ContentType','vector')

%% MODULI VS. ASTRAL NUMBER (dens=75 um^-1, 1/200 physio kT) %%

[stdevs_s,means_s] = std(moduli_s_rnd,0,2); % compute summary stats across cols
[stdevs_t,means_t] = std(moduli_t,0,2); % compute summary stats across cols

modPlot = figure(3); clf;
set(modPlot,'units','centimeters','Position',[1,1,20,22])
modPlot = tiledlayout(modPlot,2,1,'TileSpacing','loose');
nexttile(1)
hold on
for idx = 1:sweep_s_rnd.nNetTypes
    plot(sweep_s_rnd.nFilPerAsterList(idx),moduli_s_rnd(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
errorbar(sweep_s_rnd.nFilPerAsterList, means_s, ...
    z*stdevs_s/sqrt(sweep_s_rnd.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,sweep_s_rnd.nFilPerAsterList(end)])
xticks([1,4:4:sweep_s_rnd.nFilPerAsterList(end)])
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
xticks([1,4:4:sweep_t.nFilPerAsterList(end)])
ylabel({"Young's modulus"; "$Y$ [${\rm pN}\cdot\mu{\rm m^{-1}}$]"})
xlabel(modPlot,'Astral number $a_n$','FontSize',fsz,'Interpreter','latex')
exportgraphics(modPlot,fullfile(saveDir,figSubdir,['moduliVsAstralNum', ...
    filetype]),'ContentType','vector')

%% MODULUS DISTRIBUTION PLOTS %%

modDistPlot = figure(4); clf;
set(modDistPlot,'units','centimeters','Position',[1,1,18.5,20])
modDistPlot = tiledlayout(modDistPlot,2,1,'TileSpacing','loose');
% legLabels = arrayfun(@(x) sprintf('$a_n=%i$',x), ...
%     modSampling.nFilPerAsterList,'UniformOutput',false);
% modSampling.nFilPerAsterList is [1;2;3;4;8;12;16]
% chosenTypes = [1,4,5,7];    % indices of astral nums 1,4,8,16
chosenTypes = [1,4,5,6,7];  % indices of astral nums 1,4,8,12,16
numCurves = length(chosenTypes);
distColors = colormap(parula(numCurves+1));
nexttile(1)
hold on
for idx = 1:numCurves
    l = cdfplot(moduli_s_hiN(chosenTypes(idx),:));
    l.LineWidth = 2;
    l.Color = distColors(idx,:);
end
yline(0.5,'--k','LineWidth',1.5)
hold off
xlabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]', ...
    'FontSize',fsz2)
xlim([0,inf])
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
yline(0.5,'--k','LineWidth',1.5)
hold off
xlabel("Young's modulus $Y$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]", ...
    'FontSize',fsz2)
xlim([0,inf])
ylabel('')
title('')
set(gca,'FontName',fname,'FontSize',fsz)
ylabel(modDistPlot,'Empirical CDF','Interpreter','latex','FontSize', ...
    fsz,'FontName',fname)
lg = legend(string(modSampling.nFilPerAsterList(chosenTypes)),'Location','east', ...
    'FontSize',fsz2);
title(lg,'$a_n$')
lg.Layout.Tile = 'east';
exportgraphics(modDistPlot,fullfile(saveDir,figSubdir, ...
    ['modECDFs',filetype]),'ContentType','vector')

%% EXAMPLE HISTOGRAMS OF MODULI %%

% modSampling.nFilPerAsterList is [1;2;3;4;8;12;16]

shear1hist = figure(5); clf;
set(shear1hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_s_hiN(1,:),'Normalization','probability');
h.FaceColor = distColors(1,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',fszInset)
% ylabel('Prob. density')
exportgraphics(shear1hist,fullfile(saveDir,figSubdir, ...
    ['shear1hist',filetype]),'ContentType','vector')

shear16hist = figure(6); clf;
set(shear16hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_s_hiN(7,:),'Normalization','probability');
h.FaceColor = distColors(5,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',fszInset)
% ylabel('Prob. density')
exportgraphics(shear16hist,fullfile(saveDir,figSubdir, ...
    ['shear16hist',filetype]),'ContentType','vector')

tens1hist = figure(7); clf;
set(tens1hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_t_hiN(1,:),'Normalization','probability');
h.FaceColor = distColors(1,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',fszInset)
% ylabel('Prob. density')
exportgraphics(tens1hist,fullfile(saveDir,figSubdir, ...
    ['tens1hist',filetype]),'ContentType','vector')

tens16hist = figure(8); clf;
set(tens16hist,'units','centimeters','Position',[1,1,4,4])
h = histogram(moduli_t_hiN(7,:),'Normalization','probability');
h.FaceColor = distColors(5,:);
h.FaceAlpha = 1;
set(gca,'FontName',fname,'FontSize',fszInset)
% ylabel('Prob. density')
exportgraphics(tens16hist,fullfile(saveDir,figSubdir, ...
    ['tens16hist',filetype]),'ContentType','vector')

%% PROBABILITY OF RELATIVE WEAKNESS %%

beta = 0.5;
med_s = median(moduli_s_hiN,2);
med_t = median(moduli_t_hiN,2);
probPartition_s = zeros(3,modSampling.nNetTypes);
probPartition_t = zeros(3,modSampling.nNetTypes);
for idx = 1:modSampling.nNetTypes
    probPartition_s(1,idx) = sum(moduli_s_hiN(idx,:) < beta*med_s(idx)) / ...
        modSampling.numRep;
    probPartition_s(3,idx) = sum(moduli_s_hiN(idx,:) > (1/beta)*med_s(idx)) /...
        modSampling.numRep;
    probPartition_t(1,idx) = sum(moduli_t_hiN(idx,:) < beta*med_t(idx)) / ...
        modSampling.numRep;
    probPartition_t(3,idx) = sum(moduli_t_hiN(idx,:) > (1/beta)*med_t(idx))/...
        modSampling.numRep;
end
probPartition_s(2,:) = ones(1,modSampling.nNetTypes) - ...
    probPartition_s(1,:) - probPartition_s(3,:);
probPartition_t(2,:) = ones(1,modSampling.nNetTypes) - ...
    probPartition_t(1,:) - probPartition_t(3,:);

relativelyWeak = figure(9); clf;
set(relativelyWeak,'units','centimeters','Position',[1,1,15,20], ...
    'defaultLineLineWidth',2)
relativelyWeak = tiledlayout(relativelyWeak,2,1,'TileSpacing','loose');
nexttile(1)
set(gca,'FontName',fname,'FontSize',fsz)
hold on
plot(modSampling.nFilPerAsterList(1:7),probPartition_s(1,1:7),'o-', ...
    'DisplayName','Weak','Color',colors(1,:),'MarkerFaceColor', ...
    colors(1,:),'MarkerEdgeColor',colors(1,:))
text(1.5,0.18,sprintf('$\\Pr\\{G < %1.1f\\cdot G_{\\rm med}\\}$',beta), ...
    'Color',colors(1,:),'FontSize',fszInset)
hold off
xticks([1,4:4:16])
xlim([1,16])
% ylim([0,0.4])
ylabel('Shear','FontSize',fsz)

nexttile(2)
set(gca,'FontName',fname,'FontSize',fsz)
hold on
plot(modSampling.nFilPerAsterList(1:7),probPartition_t(1,1:7),'o-', ...
    'DisplayName','Weak','Color',colors(1,:),'MarkerFaceColor', ...
    colors(1,:),'MarkerEdgeColor',colors(1,:))
text(1.5,0.36,sprintf('$\\Pr\\{Y < %1.1f\\cdot Y_{\\rm med}\\}$',beta), ...
    'Color',colors(1,:),'FontSize',fszInset)
hold off
xticks([1,4:4:16])
xlim([1,16])
% ylim([0,0.4])
ylabel("Young's",'FontSize',fsz)
xlabel(relativelyWeak,'Astral number $a_n$', ...
    'FontName',fname,'FontSize',fsz,'Interpreter','latex')
ylabel(relativelyWeak,'Proportion weaker than half of median', ...
    'FontName',fname,'FontSize',fsz,'Interpreter','latex')
exportgraphics(relativelyWeak,fullfile(saveDir,figSubdir, ...
    ['relativelyWeak',filetype]),'ContentType','vector')

%% DANGLING ENDS DISTRIBUTIONS %%

endsHists = figure(10); clf;
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
    % dangling end data computed for l = 1, D = 10 but Cytosim simulations
    % are scaled down by a length factor of 10
    histogram((sweep_s_rnd.len_fil/danglingParams.l) * theseEnds, ...
        'Normalization','probability','BinWidth',0.01);
    ylim([0,0.4])
    title(sprintf('$a_n = %i$',astralNum))
    set(gca,'FontName',fname,'FontSize',fszInset)
end
xlabel(endsHists,'Dangling end length [$\mu\rm{m}$]','Interpreter', ...
    'latex','FontSize',fsz2,'FontName',fname)
ylabel(endsHists,'Probability','FontName',fname,'FontSize',fsz2)
exportgraphics(endsHists,fullfile(saveDir,figSubdir, ...
    ['endsHists',filetype]),'ContentType','vector')

%% MEAN, TOTAL DANGLING ENDS %%

linewidth = 1.5;
% useFraction statistics
[useSTD,useMean] = std(useFraction,0,2);

% dangling end data computed for l = 1, D = 10 but Cytosim simulations
% are scaled down by a length factor of 10
lenScale = (sweep_s_rnd.len_fil/danglingParams.l);

endsSummary = figure(11); clf;
set(endsSummary,'units','centimeters','Position',[1,1,15,12])
colororder([0 0 0; 0.8500 0.3250 0.0980])
yyaxis left
hold on
l1 = errorbar(danglingParams.nFilPerAsterList,lenScale*mean(endsMeans,2), ...
    lenScale*z*stdOfMeans/sqrt(danglingParams.numNets),'LineStyle','-', ...
    'Color','k','LineWidth',linewidth);
l2 = errorbar(danglingParams.nFilPerAsterList,lenScale*mean(endsMeansTrunc,2), ...
    lenScale*z*stdOfMeansTrunc/sqrt(danglingParams.numNets),'LineStyle','-', ...
    'Color',colors(3,:),'LineWidth',linewidth);
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
set(gca,'FontName',fname,'FontSize',fsz2)
exportgraphics(endsSummary,fullfile(saveDir,figSubdir, ...
    ['endsSummary',filetype]),'ContentType','vector')

%% PERCOLATION PROBABILITY HEATMAPS - RANDOM ANGLE ASTERS %%

percFsz = 16;
percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
chosenTypes = [1,5];
numDensVals = length(percParams_rnd.densityRange);
heatmaps = zeros(24,numDensVals,2);
for idx = 1:24
    theseCurves = curves_rnd.(networkLabel(idx));
    for pdx = 1:2
        heatmaps(idx,:,pdx) = clip(spline(log10(theseCurves(1,:)), ...
            theseCurves(1+chosenTypes(pdx),:), ...
            log10(percParams_rnd.densityRange)),0,1);
    end
end
% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
densityScaleFactor = percParams_rnd.l / sweep_s_rnd.len_fil;

percMaps_rnd = figure(12); clf;
set(percMaps_rnd,'units','centimeters','Position',[1,1,16,24])
percMaps_rnd = tiledlayout(percMaps_rnd,2,1,'TileSpacing','loose');
nexttile(1)
pcolor(densityScaleFactor*percParams_rnd.densityRange,1:24,heatmaps(:,:,1))
xscale('log')
yticks([1,6,12,18,24])
title('Spanning percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
nexttile(2)
pcolor(densityScaleFactor*percParams_rnd.densityRange,1:24,heatmaps(:,:,2))
xscale('log')
yticks([1,6,12,18,24])
title('Connectivity percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
xlabel(percMaps_rnd,'Density $\rho$ [$\mu{\rm m^{-1}}$]','Interpreter', ...
    'latex','FontSize',percFsz,'FontName',fname)
ylabel(percMaps_rnd,'Astral number $a_n$','Interpreter','latex', ...
    'FontSize',percFsz,'FontName',fname)
bar = colorbar;
title(bar,'Probability')
bar.Layout.Tile = 'east';
exportgraphics(percMaps_rnd,fullfile(saveDir,figSubdir,['percMaps_rnd', ...
    filetype]),'ContentType','vector')

%% PERCOLATION THRESHOLD ESTIMATES - RANDOM ANGLE ASTERS %%

threshPlot_rnd = figure(13); clf;
set(threshPlot_rnd,'units','centimeters','Position',[1,1,20,15], ...
    'defaultLineLineWidth',2.5)
hold on
for idx = 1:percParams_rnd.numPercTypesChosen
    plot(1:24,densityScaleFactor*to50_rnd(idx,1:24), ...
        'Color',colors(idx,:),'LineStyle','-')
end
yline(60,':','60','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','bottom')
yline(75,'-','75','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','middle')
yline(90,'--','90','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','top')
hold off
xticks([1,4:4:24])
xlim([1,24])
xlabel('Astral number $a_n$')
ylim([30,1e3])
yscale('log')
ylabel('Critical density $\rho_c$ [$\mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz,'FontName',fname)
lg = legend({'Spanning','Connectivity'},'Location', ...
    'northeast','FontSize',16,'FontName',fname);
title(lg,'Random angle')
exportgraphics(threshPlot_rnd,fullfile(saveDir,figSubdir, ...
    ['threshPlot_rnd',filetype]),'ContentType','Vector')

%% CORRELATION BETWEEN deltaRHO AND MODULI - RANDOM ANGLE ASTERS %%

% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
densityScaleFactor = percParams_rnd.l / sweep_s_rnd.len_fil;

connRhoc_rnd = densityScaleFactor * transpose(to50_rnd(2,1:24));
simDensities = repmat([60,75,90],[24,1]);
deltaRho_conn_rnd = simDensities - repmat(connRhoc_rnd,[1,3]);

moduli_loDens_rnd = squeeze(moduli_s_otherRhos_rnd(:,1,:));
moduli_hiDens_rnd = squeeze(moduli_s_otherRhos_rnd(:,2,:));

% correlation coefficient over all a_n
% (except where deltaRho > 0)
allModuli = [reshape(moduli_loDens_rnd,[],1); reshape(moduli_s_rnd,[],1);
    reshape(moduli_hiDens_rnd,[],1)];
deltaRho_vector = [repmat(deltaRho_conn_rnd(:,1),otherRhoParams_rnd.numRep,1); 
    repmat(deltaRho_conn_rnd(:,2),sweep_s_rnd.numRep,1); 
    repmat(deltaRho_conn_rnd(:,3),otherRhoParams_rnd.numRep,1)];
filter = logical((allModuli>0) .* (deltaRho_vector<0));
correlationMatrix = -corrcoef(log10(allModuli(filter)), ...
    log10(-deltaRho_vector(filter)));
corrCoeff_rnd = correlationMatrix(1,2);

% correlation coefficient excluding a_n = 1 (since its threshold estimate
% is ~10x estimates from other literature, e.g. Wilhelm & Frey 2003)
allModuli_no1 = [reshape(moduli_loDens_rnd(2:24,:),[],1); 
    reshape(moduli_s_rnd(2:24,:),[],1);
    reshape(moduli_hiDens_rnd(2:24,:),[],1)];
deltaRho_vector_no1 = [repmat(deltaRho_conn_rnd(2:24,1),otherRhoParams_rnd.numRep,1); 
    repmat(deltaRho_conn_rnd(2:24,2),sweep_s_rnd.numRep,1); 
    repmat(deltaRho_conn_rnd(2:24,3),otherRhoParams_rnd.numRep,1)];
filter = logical((allModuli_no1>0) .* (deltaRho_vector_no1<0));
correlationMatrix_no1 = -corrcoef(log10(allModuli_no1(filter)), ...
    log10(-deltaRho_vector_no1(filter)));
corrCoeff_no1_rnd = correlationMatrix_no1(1,2);

corrPlot_conn_rnd = figure(14); clf;
set(corrPlot_conn_rnd,'units','centimeters','Position',[1,1,20,14])
hold on
% first 3 plot calls facilitate legend format
plot(deltaRho_conn_rnd(1,1),moduli_loDens_rnd(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_rnd(1,2),moduli_s_rnd(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_rnd(1,3),moduli_hiDens_rnd(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_conn_rnd(:,1),[1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd, ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_conn_rnd(:,2),[1,sweep_s_rnd.numRep]),moduli_s_rnd,'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_conn_rnd(:,3),[1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd, ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-150,-6.5])
xscale('log')
% ylim([-Inf,60])
yscale('log')
xlabel('Density rel. to critical density $\rho - \rho^{\rm conn}_c(a_n)$')
ylabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz2,'FontName',fname)
leglabels = {sprintf('%1.0f',otherRhoParams_rnd.dens(1)),sprintf('%1.0f', ...
    sweep_s_rnd.dens),sprintf('%1.0f',otherRhoParams_rnd.dens(2))};
lg = legend(leglabels,'Location','southeast','FontSize',fszInset);
title(lg,{'Density $\rho$'; '[${\rm \mu m^{-1}}$]'})
exportgraphics(corrPlot_conn_rnd,fullfile(saveDir,figSubdir, ...
    ['corrPlot_conn_rnd',filetype]),'ContentType','vector')

%% PERCOLATION DEFINITION VISUALS - EQUAL ANGLE ASTERS %%

tbSpanEx_eq = figure(15); clf;
set(tbSpanEx_eq,'units','centimeters','Position',[1,1,16,12])
hold on
plot([0,tbSpanParams_eq.D],[0,0],'-','LineWidth',4,'Color',uciBlue)
plot([0,tbSpanParams_eq.D],[tbSpanParams_eq.D,tbSpanParams_eq.D],'-', ...
    'LineWidth',4,'Color',uciBlue)
for idx = 1:tbSpanParams_eq.numAsters
    for jdx = 1:tbSpanParams_eq.astralNum
        plot(asters_eq_span.centers(idx,1) + tbSpanParams_eq.l * ...
            cos(asters_eq_span.orients(idx,jdx)) * (0:1), ...
            asters_eq_span.centers(idx,2) + tbSpanParams_eq.l * ...
            sin(asters_eq_span.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',3,'Color','k')
    end
end
for idx = 1:size(net_eq_span.springs,1)
    augNodeA = net_eq_span.springs(idx,1);
    augNodeB = net_eq_span.springs(idx,2);
    coords = [net_eq_span.augNodes(augNodeA,1:2);
        net_eq_span.augNodes(augNodeB,1:2)];
    plot(coords(:,1), coords(:,2),'.-','MarkerSize',18,'LineWidth',3,...
        'Color',uciOrange)
    % if augNodeA <= tbSpanParams_eq.numAsters && tbSpanParams_eq.astralNum >= 2
    %     % i.e., if augNodeA is an astral center
    %     plot(net_eq_span.augNodes(augNodeA,1),net_eq_span.augNodes(augNodeA,2), ...
    %         '.r','MarkerSize',24)
    % elseif augNodeB <= tbSpanParams_eq.numAsters && tbSpanParams_eq.astralNum >= 2
    %     % likewise for augNodeB
    %     plot(net_eq_span.augNodes(augNodeB,1),net_eq_span.augNodes(augNodeB,2), ...
    %         '.r','MarkerSize',24)
    % end
end
if tbSpanParams_eq.astralNum >= 2
    plot(asters_eq_span.centers(:,1),...
        asters_eq_span.centers(:,2),'.r','MarkerSize',24)
    % could use "scatter" to enable transparency for center marks
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - tbSpanParams_eq.D/2)) * [-1,1] + tbSpanParams_eq.D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - tbSpanParams_eq.D/2)) * [-1,1] + tbSpanParams_eq.D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
exportgraphics(tbSpanEx_eq,fullfile(saveDir,figSubdir,['tbSpanEx_eq',filetype]), ...
    'ContentType','vector')

connEx_eq = figure(16); clf;
set(connEx_eq,'units','centimeters','Position',[1,1,16,12])
hold on
plot([0,connParams_eq.D],[0,0],'-','LineWidth',4,'Color',uciBlue)
plot([0,connParams_eq.D],[connParams_eq.D,connParams_eq.D],'-','LineWidth', ...
    4,'Color',uciBlue)
for idx = 1:connParams_eq.numAsters
    for jdx = 1:connParams_eq.astralNum
        plot(asters_eq_conn.centers(idx,1) + connParams_eq.l * ...
            cos(asters_eq_conn.orients(idx,jdx)) * (0:1), ...
            asters_eq_conn.centers(idx,2) + connParams_eq.l * ...
            sin(asters_eq_conn.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',3,'Color','k')
    end
end
for idx = 1:size(net_eq_conn.springs,1)
    augNodeA = net_eq_conn.springs(idx,1);
    augNodeB = net_eq_conn.springs(idx,2);
    coords = [net_eq_conn.augNodes(augNodeA,1:2);
        net_eq_conn.augNodes(augNodeB,1:2)];
    plot(coords(:,1), coords(:,2),'.-','MarkerSize',18,'LineWidth',3, ...
        'Color',uciOrange)
    % if augNodeA <= connParams_eq.numAsters && connParams_eq.astralNum >= 2
    %     % i.e., if augNodeA is an astral center
    %     plot(net_eq_conn.augNodes(augNodeA,1),net_eq_conn.augNodes(augNodeA,2), ...
    %         '.r','MarkerSize',24)
    % elseif augNodeB <= connParams_eq.numAsters && connParams_eq.astralNum >= 2
    %     % likewise for augNodeB
    %     plot(net_eq_conn.augNodes(augNodeB,1),net_eq_conn.augNodes(augNodeB,2), ...
    %         '.r','MarkerSize',24)
    % end
end
if connParams_eq.astralNum >= 2
    plot(asters_eq_conn.centers(:,1),...
        asters_eq_conn.centers(:,2),'.r','MarkerSize',24)
    % could use "scatter" to enable transparency for center marks
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - connParams_eq.D/2)) * [-1,1] + connParams_eq.D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - connParams_eq.D/2)) * [-1,1] + connParams_eq.D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
exportgraphics(connEx_eq,fullfile(saveDir,figSubdir,['connEx_eq',filetype]), ...
    'ContentType','vector')


%% PERCOLATION PROBABILITY HEATMAPS - EQUAL ANGLE ASTERS %%

percFsz = 16;
% percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
%     "TB AND LR", "Single connected component"};
chosenTypes = [1,5];
numDensVals = length(percParams_eq.densityRange);
heatmaps = zeros(24,numDensVals,2);
for idx = 1:24
    theseCurves = curves_eq.(networkLabel(idx));
    for pdx = 1:2
        heatmaps(idx,:,pdx) = clip(spline(log10(theseCurves(1,:)), ...
            theseCurves(1+chosenTypes(pdx),:), ...
            log10(percParams_eq.densityRange)),0,1);
    end
end
% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
densityScaleFactor = percParams_eq.l / sweep_s_eq.len_fil;

percMaps_eq = figure(17); clf;
set(percMaps_eq,'units','centimeters','Position',[1,1,16,24])
percMaps_eq = tiledlayout(percMaps_eq,2,1,'TileSpacing','loose');
nexttile(1)
pcolor(densityScaleFactor*percParams_eq.densityRange,1:24,heatmaps(:,:,1))
xscale('log')
yticks([1,6,12,18,24])
title('Spanning percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
nexttile(2)
pcolor(densityScaleFactor*percParams_eq.densityRange,1:24,heatmaps(:,:,2))
xscale('log')
yticks([1,6,12,18,24])
title('Connectivity percolation')
set(gca,'FontName',fname,'FontSize',percFsz)
xlabel(percMaps_eq,'Density $\rho$ [$\mu{\rm m^{-1}}$]','Interpreter', ...
    'latex','FontSize',percFsz,'FontName',fname)
ylabel(percMaps_eq,'Astral number $a_n$','Interpreter','latex', ...
    'FontSize',percFsz,'FontName',fname)
bar = colorbar;
title(bar,'Probability')
bar.Layout.Tile = 'east';
exportgraphics(percMaps_eq,fullfile(saveDir,figSubdir,['percMaps_eq', ...
    filetype]),'ContentType','vector')

%% PERCOLATION THRESHOLD ESTIMATES - EQUAL ANGLE ASTERS %%

threshPlot_eq = figure(18); clf;
set(threshPlot_eq,'units','centimeters','Position',[1,1,20,15], ...
    'defaultLineLineWidth',2.5)
hold on
for idx = 1:percParams_eq.numPercTypesChosen
    plot(1:24,densityScaleFactor*to50_eq(idx,1:24), ...
        'Color',colors(idx,:),'LineStyle','-')
end
yline(60,':','60','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','bottom')
yline(75,'-','75','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','middle')
yline(90,'--','90','Interpreter','latex','LineWidth',2, ...
    'FontSize',fszInset,'LabelVerticalAlignment','top')
hold off
xticks([1,4:4:24])
xlim([1,24])
xlabel('Astral number $a_n$')
ylim([30,1e3])
yscale('log')
ylabel('Critical density $\rho_c$ [$\mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz,'FontName',fname)
lg = legend({'Spanning','Connectivity'},'Location', ...
    'northeast','FontSize',16,'FontName',fname);
title(lg,'Equal angle')
exportgraphics(threshPlot_eq,fullfile(saveDir,figSubdir,['threshPlot_eq', ...
    filetype]),'ContentType','Vector')

%% CORRELATION BETWEEN deltaRHO AND MODULI - EQUAL ANGLE ASTERS %%

% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
densityScaleFactor = percParams_eq.l / sweep_s_eq.len_fil;

connRhoc_eq = densityScaleFactor * transpose(to50_eq(2,1:24));
simDensities = repmat([60,75,90],[24,1]);
deltaRho_conn_eq = simDensities - repmat(connRhoc_eq,[1,3]);

moduli_loDens_eq = squeeze(moduli_s_otherRhos_eq(:,1,:));
moduli_hiDens_eq = squeeze(moduli_s_otherRhos_eq(:,2,:));

% correlation coefficient over all a_n
% (except where deltaRho > 0)
allModuli = [reshape(moduli_loDens_eq,[],1); reshape(moduli_s_eq,[],1);
    reshape(moduli_hiDens_eq,[],1)];
deltaRho_vector = [repmat(deltaRho_conn_eq(:,1),otherRhoParams_eq.numRep,1); 
    repmat(deltaRho_conn_eq(:,2),sweep_s_eq.numRep,1); 
    repmat(deltaRho_conn_eq(:,3),otherRhoParams_eq.numRep,1)];
filter = logical((allModuli>0) .* (deltaRho_vector<0));
correlationMatrix = -corrcoef(log10(allModuli(filter)), ...
    log10(-deltaRho_vector(filter)));
corrCoeff_eq = correlationMatrix(1,2);

% correlation coefficient excluding a_n = 1 (since its threshold estimate
% is ~10x estimates from other literature, e.g. Wilhelm & Frey 2003)
allModuli_no1 = [reshape(moduli_loDens_eq(2:24,:),[],1); 
    reshape(moduli_s_eq(2:24,:),[],1);
    reshape(moduli_hiDens_eq(2:24,:),[],1)];
deltaRho_vector_no1 = [repmat(deltaRho_conn_eq(2:24,1),otherRhoParams_eq.numRep,1); 
    repmat(deltaRho_conn_eq(2:24,2),sweep_s_eq.numRep,1); 
    repmat(deltaRho_conn_eq(2:24,3),otherRhoParams_eq.numRep,1)];
filter = logical((allModuli_no1>0) .* (deltaRho_vector_no1<0));
correlationMatrix_no1 = -corrcoef(log10(allModuli_no1(filter)), ...
    log10(-deltaRho_vector_no1(filter)));
corrCoeff_no1_eq = correlationMatrix_no1(1,2);

corrPlot_conn_eq = figure(19); clf;
set(corrPlot_conn_eq,'units','centimeters','Position',[1,1,20,14])
hold on
% first 3 plot calls facilitate legend format
plot(deltaRho_conn_eq(1,1),moduli_loDens_eq(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_eq(1,2),moduli_s_eq(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_eq(1,3),moduli_hiDens_eq(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_conn_eq(:,1),[1,otherRhoParams_eq.numRep]),moduli_loDens_eq, ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_conn_eq(:,2),[1,sweep_s_eq.numRep]),moduli_s_eq,'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_conn_eq(:,3),[1,otherRhoParams_eq.numRep]),moduli_hiDens_eq, ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-150,-5])
xscale('log')
% ylim([-Inf,60])
yscale('log')
xlabel('Density rel. to critical density $\rho - \rho^{\rm conn}_c(a_n)$')
ylabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz2,'FontName',fname)
leglabels = {sprintf('%1.0f',otherRhoParams_eq.dens(1)),sprintf('%1.0f', ...
    sweep_s_eq.dens),sprintf('%1.0f',otherRhoParams_eq.dens(2))};
lg = legend(leglabels,'Location','southeast','FontSize',fszInset);
title(lg,{'Density $\rho$'; '[${\rm \mu m^{-1}}$]'})
exportgraphics(corrPlot_conn_eq,fullfile(saveDir,figSubdir, ...
    ['corrPlot_conn_eq',filetype]),'ContentType','vector')

%% SUPP: TIMESERIES (1/200 physio kBT) %%

% shear
timeserShear_lokT = figure(99); clf;
set(timeserShear_lokT,'units','centimeters','Position',[1,1,22,8], ...
    'defaultLineLineWidth',1.5)
timeserShear_lokT = tiledlayout(timeserShear_lokT,1,2);
nexttile(1)
hold on
xline(50,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams.shearForce/timeserParams.D))
xline(300,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs_lokT.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,1), ...
        'Color',colors(netIdx,:))
end
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

leglabels = {'',''};
leglabels = cat(2,leglabels,cellstr(string(timeserParams.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',fszInset);
title(lg,{'Astral'; 'num. $a_n$'},'Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(50,'--k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs_lokT.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,2), ...
        'Color',colors(netIdx,:))
end
ylim([-0.02,0.02])
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

exportgraphics(timeserShear_lokT,fullfile(saveDir,figSubdir, ...
    ['timeserShear_lokT',filetype]),'ContentType','vector')

% tensile
timeserTens_lokT = figure(98); clf;
set(timeserTens_lokT,'units','centimeters','Position',[1,1,22,8], ...
    'defaultLineLineWidth',1.5)
timeserTens_lokT = tiledlayout(timeserTens_lokT,1,2);
nexttile(1)
hold on
xline(50,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams.shearForce/timeserParams.D))
xline(300,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs_lokT.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,1), ...
        'Color',colors(netIdx,:))
end
hold off
xlim('tight')
ylim([-0.04,0.04])
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

leglabels = {'',''};
leglabels = cat(2,leglabels,cellstr(string(timeserParams.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',fszInset);
title(lg,{'Astral'; 'num. $a_n$'},'Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(50,'--k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams.nNetTypes
    netLabel = sprintf(timeserParams.netLabelPat, ...
        timeserParams.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs_lokT.(netLabel);
    plot(timeserParams.frameVals*timeserParams.secsPerFrame,theseCoMs(:,2), ...
        'Color',colors(netIdx,:))
end
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

exportgraphics(timeserTens_lokT,fullfile(saveDir,figSubdir, ...
    ['timeserTens_lokT',filetype]),'ContentType','vector')

%% SUPP: DISPLACEMENT VS FORCE CURVES (1/200 physio kBT) %%

selected_an = [1,4,8,16];
nNetTypes = length(selected_an);

% shear
forceDispShear_lokT = figure(97); clf;
set(forceDispShear_lokT,'units','centimeters','Position',[1,1,25,25], ...
    'defaultLineLineWidth',1.5)
forceDispShear_lokT = tiledlayout(forceDispShear_lokT,2,2, ...
    'TileSpacing','compact');
for netIdx = 1:nNetTypes
    theseDisps = disps_s_rnd.(networkLabel(selected_an(netIdx)));
    nexttile(netIdx)
    hold on
    for repIdx = 1:sweep_s_rnd.numRep
        plot(sweep_s_rnd.forceVals,theseDisps(:,1,repIdx),'-xk','LineWidth', ...
            0.5,'Color',colors(netIdx,:))
    end
    hold off
    xlim('tight')
    ylim('tight')
    title(sprintf("$a_n = %i$",selected_an(netIdx)),'Interpreter','latex')
    set(gca,'FontSize',fszInset,'FontName',fname)
    xticks(sweep_s_rnd.forceVals)
end
xlabel(forceDispShear_lokT,'Shear force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispShear_lokT,'$x$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

exportgraphics(forceDispShear_lokT,fullfile(saveDir,figSubdir, ...
    ['forceDispShear_lokT',filetype]),'ContentType','vector')

% tensile
forceDispTens_lokT = figure(96); clf;
set(forceDispTens_lokT,'units','centimeters','Position',[1,1,25,25], ...
    'defaultLineLineWidth',1.5)
forceDispTens_lokT = tiledlayout(forceDispTens_lokT,2,2, ...
    'TileSpacing','compact');
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
    set(gca,'FontSize',fszInset,'FontName',fname)
    xticks(sweep_t.forceVals)
end
xlabel(forceDispTens_lokT,'Tensile force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispTens_lokT,'$y$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

exportgraphics(forceDispTens_lokT,fullfile(saveDir,figSubdir, ...
    ['forceDispTens_lokT',filetype]),'ContentType','vector')

%% SUPP: TIMESERIES (phys kBT) %%

% shear
timeserShear_physkT = figure(95); clf;
set(timeserShear_physkT,'units','centimeters','Position',[1,1,22,8], ...
    'defaultLineLineWidth',1.5)
timeserShear_physkT = tiledlayout(timeserShear_physkT,1,2);
nexttile(1)
hold on
xline(41,'--k','LineWidth',1.5)
xline(50,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams_physkT.shearForce/timeserParams_physkT.D))
xline(291,'-k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams_physkT.nNetTypes
    netLabel = sprintf(timeserParams_physkT.netLabelPat, ...
        timeserParams_physkT.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs_physkT.(netLabel);
    plot(timeserParams_physkT.frameVals*timeserParams_physkT.secsPerFrame, ...
        theseCoMs(:,1),'Color',colors(netIdx,:))
end
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

leglabels = {'','','',''};
leglabels = cat(2,leglabels, ...
    cellstr(string(timeserParams_physkT.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',fszInset);
title(lg,{'Astral'; 'num. $a_n$'},'Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(41,'--k','LineWidth',1.5)
xline(50,'--k','LineWidth',1.5)
xline(291,'-k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams_physkT.nNetTypes
    netLabel = sprintf(timeserParams_physkT.netLabelPat, ...
        timeserParams_physkT.nFilPerAsterList(netIdx));
    theseCoMs = shearCoMs_physkT.(netLabel);
    plot(timeserParams_physkT.frameVals*timeserParams_physkT.secsPerFrame, ...
        theseCoMs(:,2),'Color',colors(netIdx,:))
end
ylim([-0.075,0.075])
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

exportgraphics(timeserShear_physkT,fullfile(saveDir,figSubdir, ...
    ['timeserShear_physkT',filetype]),'ContentType','vector')

% tensile
timeserTens_physkT = figure(94); clf;
set(timeserTens_physkT,'units','centimeters','Position',[1,1,22,8], ...
    'defaultLineLineWidth',1.5)
timeserTens_physkT = tiledlayout(timeserTens_physkT,1,2);
nexttile(1)
hold on
xline(41,'--k','LineWidth',1.5)
xline(50,'--k','LineWidth',1.5,'DisplayName', ...
    sprintf("%2.1f [${\\rm pN/\\mu m}$] applied", ...
    timeserParams_physkT.shearForce/timeserParams_physkT.D))
xline(291,'-k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5,'DisplayName','Final position')
for netIdx = 1:timeserParams_physkT.nNetTypes
    netLabel = sprintf(timeserParams_physkT.netLabelPat, ...
        timeserParams_physkT.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs_physkT.(netLabel);
    plot(timeserParams_physkT.frameVals*timeserParams_physkT.secsPerFrame, ...
        theseCoMs(:,1),'Color',colors(netIdx,:))
end
hold off
xlim('tight')
ylim([-0.075,0.075])
xlabel('Time [s]')
ylabel('$x$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

leglabels = {'','','',''};
leglabels = cat(2,leglabels, ...
    cellstr(string(timeserParams_physkT.nFilPerAsterList)));
lg = legend(leglabels,'FontName',fname,'FontSize',fszInset);
title(lg,{'Astral'; 'num. $a_n$'},'Interpreter','latex')
lg.Layout.Tile = 'east';

nexttile(2)
hold on
xline(41,'--k','LineWidth',1.5)
xline(50,'--k','LineWidth',1.5)
xline(291,'-k','LineWidth',1.5)
xline(300,'-k','LineWidth',1.5)
for netIdx = 1:timeserParams_physkT.nNetTypes
    netLabel = sprintf(timeserParams_physkT.netLabelPat, ...
        timeserParams_physkT.nFilPerAsterList(netIdx));
    theseCoMs = tensCoMs_physkT.(netLabel);
    plot(timeserParams_physkT.frameVals*timeserParams_physkT.secsPerFrame, ...
        theseCoMs(:,2),'Color',colors(netIdx,:))
end
% ylim([-0.5,0.5])
hold off
xlim('tight')
xlabel('Time [s]')
ylabel('$y$-center of mass [${\rm \mu m}$]','FontSize',fszInset)
set(gca,'FontSize',fszInset,'FontName',fname)

exportgraphics(timeserTens_physkT,fullfile(saveDir,figSubdir, ...
    ['timeserTens_physkT',filetype]),'ContentType','vector')

%% SUPP: DISPLACEMENT VS FORCE CURVES (phys kBT) %%

selected_an = [1,4,8,16];
nNetTypes = length(selected_an);

% shear
forceDispShear_physkT = figure(93); clf;
set(forceDispShear_physkT,'units','centimeters','Position',[1,1,25,25], ...
    'defaultLineLineWidth',1.5)
forceDispShear_physkT = tiledlayout(forceDispShear_physkT,2,2, ...
    'TileSpacing','compact');
for netIdx = 1:nNetTypes
    theseDisps = disps_s_physkT.(networkLabel(selected_an(netIdx)));
    nexttile(netIdx)
    hold on
    for repIdx = 1:sweep_s_physkT.numRep
        plot(sweep_s_physkT.forceVals,theseDisps(:,1,repIdx), ...
            '-xk','LineWidth',0.5,'Color',colors(netIdx,:))
    end
    hold off
    xlim('tight')
    ylim('tight')
    title(sprintf("$a_n = %i$",selected_an(netIdx)),'Interpreter','latex')
    set(gca,'FontSize',fszInset,'FontName',fname)
    xticks(sweep_s_physkT.forceVals)
end
xlabel(forceDispShear_physkT,'Shear force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispShear_physkT,'$x$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

exportgraphics(forceDispShear_physkT,fullfile(saveDir,figSubdir, ...
    ['forceDispShear_physkT',filetype]),'ContentType','vector')

% tensile
forceDispTens_physkT = figure(92); clf;
set(forceDispTens_physkT,'units','centimeters','Position',[1,1,25,25], ...
    'defaultLineLineWidth',1.5)
forceDispTens_physkT = tiledlayout(forceDispTens_physkT,2,2, ...
    'TileSpacing','compact');
for netIdx = 1:nNetTypes
    theseDisps = disps_t_physkT.(networkLabel(selected_an(netIdx)));
    nexttile(netIdx)
    hold on
    for repIdx = 1:sweep_t_physkT.numRep
        plot(sweep_t_physkT.forceVals,theseDisps(:,2,repIdx), ...
            '-xk','LineWidth',0.5,'Color',colors(netIdx,:))
    end
    hold off
    xlim('tight')
    ylim('tight')
    title(sprintf("$a_n = %i$",selected_an(netIdx)),'Interpreter','latex')
    set(gca,'FontSize',fszInset,'FontName',fname)
    xticks(sweep_t_physkT.forceVals)
end
xlabel(forceDispTens_physkT,'Tensile force [pN]', ...
    'FontSize',fsz,'FontName',fname)
ylabel(forceDispTens_physkT,'$y$-displacement [${\rm \mu m}$]', ...
    'FontSize',fsz,'FontName',fname,'Interpreter','latex')

exportgraphics(forceDispTens_physkT,fullfile(saveDir,figSubdir, ...
    ['forceDispTens_physkT',filetype]),'ContentType','vector')

%% SUPP: MODULI VS. ASTRAL NUMBER (dens=75 um^-1, physkT) %%

[stdevs_s_physkT,means_s_physkT] = std(moduli_s_physkT,0,2);
[stdevs_t_physkT,means_t_physkT] = std(moduli_t_physkT,0,2);

modPlot_physkT = figure(91); clf;
set(modPlot_physkT,'units','centimeters','Position',[1,1,20,22])
modPlot_physkT = tiledlayout(modPlot_physkT,2,1,'TileSpacing','loose');
nexttile(1)
hold on
for idx = 1:sweep_s_physkT.nNetTypes
    plot(sweep_s_physkT.nFilPerAsterList(idx),moduli_s_physkT(idx,:), ...
        'Marker','o','MarkerSize',3,'MarkerFaceColor',lgray, ...
        'MarkerEdgeColor',lgray)
end
errorbar(sweep_s_physkT.nFilPerAsterList, means_s_physkT, ...
    z*stdevs_s_physkT/sqrt(sweep_s_physkT.numRep),'-o','LineWidth',1.5, ...
    'Color',uciBlue)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,sweep_s_physkT.nFilPerAsterList(end)])
xticks([1,4:4:sweep_s_physkT.nFilPerAsterList(end)])
% xlabel('Astral number $a_n$')
ylabel({'Shear modulus'; '$G$ [${\rm pN}\cdot\mu{\rm m^{-1}}$]'})

nexttile(2)
hold on
for idx = 1:sweep_t_physkT.nNetTypes
    plot(sweep_t_physkT.nFilPerAsterList(idx),moduli_t_physkT(idx,:), ...
        'Marker','o','MarkerSize',3,'MarkerFaceColor',lgray, ...
        'MarkerEdgeColor', lgray)
end
errorbar(sweep_t_physkT.nFilPerAsterList, means_t_physkT, ...
    z*stdevs_t_physkT/sqrt(sweep_t_physkT.numRep),'-o','LineWidth',1.5, ...
    'Color',uciBlue)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,sweep_t_physkT.nFilPerAsterList(end)])
xticks([1,4:4:sweep_t_physkT.nFilPerAsterList(end)])
ylabel({"Young's modulus"; "$Y$ [${\rm pN}\cdot\mu{\rm m^{-1}}$]"})
xlabel(modPlot_physkT,'Astral number $a_n$','FontSize',fsz,'Interpreter','latex')
exportgraphics(modPlot_physkT,fullfile(saveDir,figSubdir, ...
    ['moduliVsAstralNum_physkT',filetype]),'ContentType','vector')

%% SUPP: SHEAR MODULI VS. ASTRAL NUMBER AT 2 ADD'L FILAMENT DENSITIES %%
% random angle asters

[stdevMod_otherRhos_rnd,meanMod_otherRhos_rnd] = std(moduli_s_otherRhos_rnd,0,3);
markerList = {'v','^'};
otherRhos_rnd = figure(90); clf;
set(otherRhos_rnd,'units','centimeters','Position',[1,1,20,14])
hold on
for densIdx = 1:otherRhoParams_rnd.numDens
    for netIdx = 1:otherRhoParams_rnd.nNetTypes
        plot(otherRhoParams_rnd.nFilPerAsterList(netIdx), ...
            squeeze(moduli_s_otherRhos_rnd(netIdx,densIdx,:)), ...
            'Marker',markerList{densIdx},'MarkerSize',3,'MarkerFaceColor', ...
            lgray,'MarkerEdgeColor',lgray)
    end
end
for idx = 1:sweep_s_rnd.nNetTypes
    plot(sweep_s_rnd.nFilPerAsterList(idx),moduli_s_rnd(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
l1 = errorbar(otherRhoParams_rnd.nFilPerAsterList,meanMod_otherRhos_rnd(:,1), ...
    z*stdevMod_otherRhos_rnd(:,1)/sqrt(otherRhoParams_rnd.numRep),'Color','k', ...
    'Marker',markerList{1},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams_rnd.dens(1)));
l2 = errorbar(otherRhoParams_rnd.nFilPerAsterList,meanMod_otherRhos_rnd(:,2), ...
    z*stdevMod_otherRhos_rnd(:,2)/sqrt(otherRhoParams_rnd.numRep),'Color', ...
    uciOrange,'Marker',markerList{2},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams_rnd.dens(2)));
l3 = errorbar(sweep_s_rnd.nFilPerAsterList, means_s, ...
    z*stdevs_s/sqrt(sweep_s_rnd.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue,'DisplayName',sprintf('%1.1f',sweep_s_rnd.dens));
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,otherRhoParams_rnd.nFilPerAsterList(end)])
xticks([1,4:4:otherRhoParams_rnd.nFilPerAsterList(end)])
xlabel('Astral number $a_n$')
yl = ylim;
ylim([0,yl(2)])
ylabel({'Shear modulus $G$'; '[${\rm pN} \cdot \mu{\rm m}^{-1}$]'})
lg = legend([l1,l3,l2],'Location','northeast','FontSize',fsz2);
title(lg,{'Density $\rho$';'[${\rm \mu m^{-1}}$]'},'FontSize',fsz2)
exportgraphics(otherRhos_rnd,fullfile(saveDir,figSubdir, ...
    ['otherRhos_rnd',filetype]),'ContentType','vector')

%% SUPP: SHEAR MODULI VS. ASTRAL NUMBER AT 2 ADD'L FILAMENT DENSITIES %%
% equal angle asters

[stdevs_s_eq,means_s_eq] = std(moduli_s_eq,0,2);
[stdevMod_otherRhos_eq,meanMod_otherRhos_eq] = std(moduli_s_otherRhos_eq,0,3);
markerList = {'v','^'};
otherRhos_eq = figure(89); clf;
set(otherRhos_eq,'units','centimeters','Position',[1,1,20,14])
hold on
for densIdx = 1:otherRhoParams_eq.numDens
    for netIdx = 1:otherRhoParams_eq.nNetTypes
        plot(otherRhoParams_eq.nFilPerAsterList(netIdx), ...
            squeeze(moduli_s_otherRhos_eq(netIdx,densIdx,:)), ...
            'Marker',markerList{densIdx},'MarkerSize',3,'MarkerFaceColor', ...
            lgray,'MarkerEdgeColor',lgray)
    end
end
for idx = 1:sweep_s_eq.nNetTypes
    plot(sweep_s_eq.nFilPerAsterList(idx),moduli_s_eq(idx,:),'Marker','o', ...
        'MarkerSize',3,'MarkerFaceColor',lgray,'MarkerEdgeColor', lgray)
end
l1 = errorbar(otherRhoParams_eq.nFilPerAsterList,meanMod_otherRhos_eq(:,1), ...
    z*stdevMod_otherRhos_eq(:,1)/sqrt(otherRhoParams_eq.numRep),'Color','k', ...
    'Marker',markerList{1},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams_eq.dens(1)));
l2 = errorbar(otherRhoParams_eq.nFilPerAsterList,meanMod_otherRhos_eq(:,2), ...
    z*stdevMod_otherRhos_eq(:,2)/sqrt(otherRhoParams_eq.numRep),'Color', ...
    uciOrange,'Marker',markerList{2},'LineWidth',1.5,'DisplayName', ...
    sprintf('%1.1f',otherRhoParams_eq.dens(2)));
l3 = errorbar(sweep_s_eq.nFilPerAsterList, means_s_eq, ...
    z*stdevs_s_eq/sqrt(sweep_s_eq.numRep), '-o', 'LineWidth', 1.5, 'Color', ...
    uciBlue,'DisplayName',sprintf('%1.1f',sweep_s_eq.dens));
hold off
set(gca,'FontName',fname,'FontSize',fsz)
xlim([0,otherRhoParams_eq.nFilPerAsterList(end)])
xticks([1,4:4:otherRhoParams_eq.nFilPerAsterList(end)])
xlabel('Astral number $a_n$')
yl = ylim;
ylim([0,yl(2)])
ylabel({'Shear modulus $G$'; '[${\rm pN} \cdot \mu{\rm m}^{-1}$]'})
lg = legend([l1,l3,l2],'Location','northeast','FontSize',fsz2);
title(lg,{'Density $\rho$';'[${\rm \mu m^{-1}}$]'},'FontSize',fsz2)
exportgraphics(otherRhos_eq,fullfile(saveDir,figSubdir, ...
    ['otherRhos_eq',filetype]),'ContentType','vector')

%% SUPP: COMBINED CONN CORRELATION PLOT %%

alpha = 0.4;
corrPlot_conn_rndeq = figure(88); clf;
set(corrPlot_conn_rndeq,'units','centimeters','Position',[1,1,20,14])
hold on
% first 6 scatter/plot calls facilitate legend format
scatter(deltaRho_conn_rnd(1,1),moduli_loDens_rnd(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineWidth',0.75, ...
    'MarkerFaceColor','k','MarkerFaceAlpha',alpha)
scatter(deltaRho_conn_rnd(1,2),moduli_s_rnd(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75, ...
    'MarkerFaceColor',uciBlue,'MarkerFaceAlpha',alpha)
scatter(deltaRho_conn_rnd(1,3),moduli_hiDens_rnd(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineWidth',0.75, ...
    'MarkerFaceColor',uciOrange,'MarkerFaceAlpha',alpha)
plot(deltaRho_conn_eq(1,1),moduli_loDens_eq(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_eq(1,2),moduli_s_eq(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_conn_eq(1,3),moduli_hiDens_eq(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_conn_rnd(:,1),[1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd, ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75,'MarkerFaceColor', ...
    'k','MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_conn_rnd(:,2),[1,sweep_s_rnd.numRep]),moduli_s_rnd,'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75,'MarkerFaceColor',uciBlue, ...
    'MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_conn_rnd(:,3),[1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd, ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75, ...
    'MarkerFaceColor',uciOrange,'MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_conn_eq(:,1),[1,otherRhoParams_eq.numRep]),moduli_loDens_eq, ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_conn_eq(:,2),[1,sweep_s_eq.numRep]),moduli_s_eq,'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_conn_eq(:,3),[1,otherRhoParams_eq.numRep]),moduli_hiDens_eq, ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-160,-5])
xscale('log')
% ylim([-Inf,60])
yscale('log')
xlabel('Density rel. to critical density $\rho - \rho^{\rm conn}_c(a_n)$')
ylabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
set(gca,'FontSize',fsz2,'FontName',fname)
leglabels = {sprintf('rnd. %1.1f',otherRhoParams_rnd.dens(1)), ...
    sprintf('rnd. %1.1f',sweep_s_rnd.dens), ...
    sprintf('rnd. %1.1f',otherRhoParams_rnd.dens(2)), ...
    sprintf('eq. %1.1f',otherRhoParams_rnd.dens(1)), ...
    sprintf('eq. %1.1f',sweep_s_rnd.dens), ...
    sprintf('eq. %1.1f',otherRhoParams_rnd.dens(2))};
lg = legend(leglabels,'Location','southeast','FontSize',fszInset, ...
    'NumColumns',2);
title(lg,{'Density $\rho$ [${\rm \mu m^{-1}}$]'})
exportgraphics(corrPlot_conn_rndeq,fullfile(saveDir,figSubdir, ...
    ['corrPlot_conn_rndeq',filetype]),'ContentType','vector')

%% SUPP: CORRELATION PLOTS W.R.T. SPANNING THRESHOLDS %%

% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
densityScaleFactor = percParams_rnd.l / sweep_s_rnd.len_fil;
% densityScaleFactor = percParams_eq.l / sweep_s_eq.len_fil; % redundant
simDensities = repmat([60,75,90],[24,1]);

spanRhoc_rnd = densityScaleFactor * transpose(to50_rnd(1,1:24));
deltaRho_span_rnd = simDensities - repmat(spanRhoc_rnd,[1,3]);
spanRhoc_eq = densityScaleFactor * transpose(to50_eq(1,1:24));
deltaRho_span_eq = simDensities - repmat(spanRhoc_eq,[1,3]);

corrPlot_span_rnd = figure(87); clf;
set(corrPlot_span_rnd,'units','centimeters','Position',[1,1,25,10])
corrPlot_span_rnd = tiledlayout(1,2,"TileSpacing","compact");
negFilter_rnd = deltaRho_span_rnd < 0;
nexttile(1) % deltaRho_span_rnd < 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
% first 3 plot calls facilitate legend format
plot(deltaRho_span_rnd(1,1),moduli_loDens_rnd(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_rnd(9,2),moduli_s_rnd(9,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_rnd(12,3),moduli_hiDens_rnd(12,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,1),1), ...
    [1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd(negFilter_rnd(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,2),2), ...
    [1,sweep_s_rnd.numRep]),moduli_s_rnd(negFilter_rnd(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,3),3), ...
    [1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd(negFilter_rnd(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
% xlim([-150,-6.5])
xscale('log')
ylim([0.6,41])
yscale('log')
nexttile(2) % deltaRho_span_rnd > 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
plot(deltaRho_span_rnd(2,1),moduli_loDens_rnd(2,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_rnd(1,2),moduli_s_rnd(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_rnd(1,3),moduli_hiDens_rnd(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,1),1), ...
    [1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd(~negFilter_rnd(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,2),2), ...
    [1,sweep_s_rnd.numRep]),moduli_s_rnd(~negFilter_rnd(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,3),3), ...
    [1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd(~negFilter_rnd(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
% xlim([-150,-6.5])
xscale('log')
ylim([0.6,41])
yscale('log')

xlabel(corrPlot_span_rnd, ...
    'Density rel. to critical density $\rho - \rho^{\rm span}_c(a_n)$', ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
ylabel(corrPlot_span_rnd, ...
    {'Shear modulus $G$'; '[${\rm pN}\cdot \mu{\rm m^{-1}}$]'}, ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
leglabels = {sprintf('%1.0f',otherRhoParams_rnd.dens(1)),sprintf('%1.0f', ...
    sweep_s_rnd.dens),sprintf('%1.0f',otherRhoParams_rnd.dens(2))};
lg = legend(leglabels,'FontSize',fszInset);
title(lg,{'Density $\rho$'; '[${\rm \mu m^{-1}}$]'})
lg.Layout.Tile = 'east';
exportgraphics(corrPlot_span_rnd,fullfile(saveDir,figSubdir, ...
    ['corrPlot_span_rnd',filetype]),'ContentType','vector')

corrPlot_span_eq = figure(86); clf;
set(corrPlot_span_eq,'units','centimeters','Position',[1,1,25,10])
corrPlot_span_eq = tiledlayout(1,2,"TileSpacing","compact");
negFilter_eq = deltaRho_span_eq < 0;
nexttile(1) % deltaRho_span_eq < 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
% first 3 plot calls facilitate legend format
plot(deltaRho_span_eq(1,1),moduli_loDens_eq(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(10,2),moduli_s_eq(10,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(13,3),moduli_hiDens_eq(13,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,1),1), ...
    [1,otherRhoParams_eq.numRep]),moduli_loDens_eq(negFilter_eq(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,2),2), ...
    [1,sweep_s_eq.numRep]),moduli_s_eq(negFilter_eq(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,3),3), ...
    [1,otherRhoParams_eq.numRep]),moduli_hiDens_eq(negFilter_eq(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-100,-0.12])
xscale('log')
ylim([0.6,65])
yscale('log')
nexttile(2) % deltaRho_span_eq > 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
plot(deltaRho_span_eq(2,1),moduli_loDens_eq(2,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(1,2),moduli_s_eq(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(1,3),moduli_hiDens_eq(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,1),1), ...
    [1,otherRhoParams_eq.numRep]),moduli_loDens_eq(~negFilter_eq(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,2),2), ...
    [1,sweep_s_eq.numRep]),moduli_s_eq(~negFilter_eq(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,3),3), ...
    [1,otherRhoParams_eq.numRep]),moduli_hiDens_eq(~negFilter_eq(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
% xlim([-150,-6.5])
xscale('log')
ylim([0.6,65])
yscale('log')

xlabel(corrPlot_span_eq, ...
    'Density rel. to critical density $\rho - \rho^{\rm span}_c(a_n)$', ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
ylabel(corrPlot_span_eq, ...
    {'Shear modulus $G$'; '[${\rm pN}\cdot \mu{\rm m^{-1}}$]'}, ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
leglabels = {sprintf('%1.0f',otherRhoParams_eq.dens(1)),sprintf('%1.0f', ...
    sweep_s_eq.dens),sprintf('%1.0f',otherRhoParams_eq.dens(2))};
lg = legend(leglabels,'FontSize',fszInset);
title(lg,{'Density $\rho$'; '[${\rm \mu m^{-1}}$]'})
lg.Layout.Tile = 'east';
exportgraphics(corrPlot_span_eq,fullfile(saveDir,figSubdir, ...
    ['corrPlot_span_eq',filetype]),'ContentType','vector')

alpha = 0.4;
corrPlot_span_rndeq = figure(85); clf;
set(corrPlot_span_rndeq,'units','centimeters','Position',[1,1,28.5,10])
corrPlot_span_rndeq = tiledlayout(1,2,"TileSpacing","compact");
nexttile(1)     % deltaRho_span_* < 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
% first 6 scatter/plot calls facilitate legend format
scatter(deltaRho_span_rnd(1,1),moduli_loDens_rnd(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineWidth',0.75, ...
    'MarkerFaceColor','k','MarkerFaceAlpha',alpha)
scatter(deltaRho_span_rnd(9,2),moduli_s_rnd(9,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75, ...
    'MarkerFaceColor',uciBlue,'MarkerFaceAlpha',alpha)
scatter(deltaRho_span_rnd(12,3),moduli_hiDens_rnd(12,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineWidth',0.75, ...
    'MarkerFaceColor',uciOrange,'MarkerFaceAlpha',alpha)
plot(deltaRho_span_eq(1,1),moduli_loDens_eq(1,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(10,2),moduli_s_eq(10,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(13,3),moduli_hiDens_eq(13,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,1),1), ...
    [1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd(negFilter_rnd(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75, ...
    'MarkerFaceColor','k','MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,2),2), ...
    [1,sweep_s_rnd.numRep]),moduli_s_rnd(negFilter_rnd(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75,'MarkerFaceColor',uciBlue, ...
    'MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_span_rnd(negFilter_rnd(:,3),3), ...
    [1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd(negFilter_rnd(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75, ...
    'MarkerFaceColor',uciOrange,'MarkerFaceAlpha',alpha)
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,1),1), ...
    [1,otherRhoParams_eq.numRep]),moduli_loDens_eq(negFilter_eq(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,2),2), ...
    [1,sweep_s_eq.numRep]),moduli_s_eq(negFilter_eq(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(negFilter_eq(:,3),3), ...
    [1,otherRhoParams_eq.numRep]),moduli_hiDens_eq(negFilter_eq(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([-100,-0.13])
xscale('log')
ylim([0.6,65])
yscale('log')
nexttile(2) % deltaRho_span_eq > 0
set(gca,'FontSize',fsz2,'FontName',fname)
hold on
scatter(deltaRho_span_rnd(2,1),moduli_loDens_rnd(2,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineWidth',0.75, ...
    'MarkerFaceColor','k','MarkerFaceAlpha',alpha)
scatter(deltaRho_span_rnd(1,2),moduli_s_rnd(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75, ...
    'MarkerFaceColor',uciBlue,'MarkerFaceAlpha',alpha)
scatter(deltaRho_span_rnd(1,3),moduli_hiDens_rnd(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineWidth',0.75, ...
    'MarkerFaceColor',uciOrange,'MarkerFaceAlpha',alpha)
plot(deltaRho_span_eq(2,1),moduli_loDens_eq(2,1),'Marker','v', ...
    'MarkerEdgeColor','k','LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(1,2),moduli_s_eq(1,1),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineStyle','none','LineWidth',0.75)
plot(deltaRho_span_eq(1,3),moduli_hiDens_eq(1,1),'Marker','^', ...
    'MarkerEdgeColor',uciOrange,'LineStyle','none','LineWidth',0.75)
% these are the full scatter plots
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,1),1), ...
    [1,otherRhoParams_rnd.numRep]),moduli_loDens_rnd(~negFilter_rnd(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,2),2), ...
    [1,sweep_s_rnd.numRep]),moduli_s_rnd(~negFilter_rnd(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_rnd(~negFilter_rnd(:,3),3), ...
    [1,otherRhoParams_rnd.numRep]),moduli_hiDens_rnd(~negFilter_rnd(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,1),1), ...
    [1,otherRhoParams_eq.numRep]),moduli_loDens_eq(~negFilter_eq(:,1),:), ...
    'Marker','v','MarkerEdgeColor','k','LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,2),2), ...
    [1,sweep_s_eq.numRep]),moduli_s_eq(~negFilter_eq(:,2),:),'Marker','o', ...
    'MarkerEdgeColor',uciBlue,'LineWidth',0.75)
scatter(repmat(deltaRho_span_eq(~negFilter_eq(:,3),3), ...
    [1,otherRhoParams_eq.numRep]),moduli_hiDens_eq(~negFilter_eq(:,3),:), ...
    'Marker','^','MarkerEdgeColor',uciOrange,'LineWidth',0.75)
hold off
xlim([0.5,50])
xscale('log')
ylim([0.6,65])
yscale('log')

xlabel(corrPlot_span_rndeq, ...
    'Density rel. to critical density $\rho - \rho^{\rm span}_c(a_n)$', ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
ylabel(corrPlot_span_rndeq, ...
    {'Shear modulus $G$'; '[${\rm pN}\cdot \mu{\rm m^{-1}}$]'}, ...
    'Interpreter','latex','FontSize',fsz,'FontName',fname)
leglabels = {sprintf('rnd. %1.1f',otherRhoParams_rnd.dens(1)), ...
    sprintf('rnd. %1.1f',sweep_s_rnd.dens), ...
    sprintf('rnd. %1.1f',otherRhoParams_rnd.dens(2)), ...
    sprintf('eq. %1.1f',otherRhoParams_rnd.dens(1)), ...
    sprintf('eq. %1.1f',sweep_s_rnd.dens), ...
    sprintf('eq. %1.1f',otherRhoParams_rnd.dens(2))};
lg = legend(leglabels,'FontSize',fszInset,'NumColumns',2);
title(lg,'Density $\rho$ [${\rm \mu m^{-1}}$]')
lg.Layout.Tile = 'east';
exportgraphics(corrPlot_span_rndeq,fullfile(saveDir,figSubdir, ...
    ['corrPlot_span_rndeq',filetype]),'ContentType','vector')

%% SUPP: CROSSLINKER DENSITY %%

[stdevMod_crslnkDens,meanMod_crslnkDens] = std(moduli_s_crslnkDens,0,3);

crslnkDensShear_famOfDens = figure(84); clf;
set(crslnkDensShear_famOfDens,'units','centimeters','Position',[1,1,20,16])
cmap = colormap(parula(crslnkDensParams.nCrslnkDens+1));
hold on
for jdx = 1:crslnkDensParams.nCrslnkDens
    errorbar(crslnkDensParams.nFilPerAsterList,meanMod_crslnkDens(:,jdx), ...
        (z/sqrt(crslnkDensParams.numRep))*stdevMod_crslnkDens(:,jdx), ...
        'Color',cmap(jdx,:),'LineWidth',1.5)
end
hold off
lg = legend(string(crslnkDensParams.crosslinker_dens_list), ...
    'Location','northeast');
title(lg,{"Crosslinker density"; "[${\rm particles}\cdot\mu {\rm m}^{-2}$]"})
xlabel('Astral number $a_n$')
xticks(crslnkDensParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
set(gca,'FontName',fname,'FontSize',fsz)
exportgraphics(crslnkDensShear_famOfDens,fullfile(saveDir,figSubdir, ...
    ['crslnkDensShear_famOfDens',filetype]),'ContentType','vector')

crslnkDensShear_famOfan = figure(83); clf;
set(crslnkDensShear_famOfan,'units','centimeters','Position',[1,1,20,16])
an_to_show = [1,2,4,5,7];     % indices of an = 1,2,4,5,8
cmap = colormap(parula(length(an_to_show)+1));
hold on
for idx = 1:length(an_to_show)
    errorbar(crslnkDensParams.crosslinker_dens_list, ...
        meanMod_crslnkDens(an_to_show(idx),:), ...
        (z/sqrt(crslnkDensParams.numRep)) * ...
        stdevMod_crslnkDens(an_to_show(idx),:), ...
        'Color',cmap(idx,:),'LineWidth',1.5)
end
set(gca,'FontName',fname,'FontSize',fsz)
lg = legend(string(crslnkDensParams.nFilPerAsterList(an_to_show)), ...
    'Location','south','FontSize',fszInset,'NumColumns',5);
title(lg,"Astral number $a_n$")
xlim([7e3,3.85e4])
xlabel({"Crosslinker density"; "[${\rm particles}\cdot\mu {\rm m}^{-2}$]"})
ylim([0,inf])
ylabel('Shear modulus $[{\rm pN}\cdot \mu{\rm m^{-1}}]$')
exportgraphics(crslnkDensShear_famOfan,fullfile(saveDir,figSubdir, ...
    ['crslnkDensShear_famOfan',filetype]),'ContentType','vector')

%% SUPP: CROSSLINKER STIFFNESS %%

[stdevMod_crslnkStiff,meanMod_crslnkStiff] = std(moduli_s_crslnkStiff,0,3);

crslnkStiffShear_famOfStiff = figure(82); clf;
set(crslnkStiffShear_famOfStiff,'units','centimeters','Position',[1,1,20,16])
cmap = colormap(parula(crslnkStiffParams.nCrslnkStiff+1));
hold on
for jdx = 1:crslnkStiffParams.nCrslnkStiff
    errorbar(crslnkStiffParams.nFilPerAsterList,meanMod_crslnkStiff(:,jdx), ...
        (z/sqrt(crslnkStiffParams.numRep))*stdevMod_crslnkStiff(:,jdx), ...
        'Color',cmap(jdx,:),'LineWidth',1.5)
end
hold off
lg = legend(string(crslnkStiffParams.crosslinker_stiff_list), ...
    'Location','northeast');
title(lg,{"$\:$Crosslinker stiffness$\:$"; "[${\rm pN}\cdot \mu {\rm m}^{-1}$]"})
xlabel('Astral number $a_n$')
xticks(crslnkStiffParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot\mu{\rm m}^{-1}$]')
set(gca,'FontName',fname,'FontSize',fsz)
exportgraphics(crslnkStiffShear_famOfStiff,fullfile(saveDir,figSubdir, ...
    ['crslnkStiffShear_famOfStiff',filetype]),'ContentType','vector')

crslnkStiffShear_famOfan = figure(81); clf;
set(crslnkStiffShear_famOfan,'units','centimeters','Position',[1,1,20,16])
an_to_show = [1,2,4,5,7];     % indices of an = 1,2,4,5,8
cmap = colormap(parula(length(an_to_show)+1));
hold on
for idx = 1:length(an_to_show)
    errorbar(crslnkStiffParams.crosslinker_stiff_list, ...
        meanMod_crslnkStiff(an_to_show(idx),:), ...
        (z/sqrt(crslnkStiffParams.numRep)) * ...
        stdevMod_crslnkStiff(an_to_show(idx),:), ...
        'Color',cmap(idx,:),'LineWidth',1.5)
end
reflineX = logspace(log10(55),2.1,5);
plot(reflineX,3.5 * reflineX.^(0.5),'--k','LineWidth',3)
text(70,37,"0.5",'FontName',fname,'FontSize',fsz2)
hold off
set(gca,'FontName',fname,'FontSize',fsz)
lg = legend(string(crslnkStiffParams.nFilPerAsterList(an_to_show)), ...
    'Location','northwest','FontSize',fszInset,'NumColumns',1);
title(lg,"$\:$Astral number $a_n\:$")
xscale('log')
xlim([10,450])
yscale('log')
ylim([2.5,60])
xlabel('Crosslinker stiffness [${\rm pN}\cdot \mu {\rm m}^{-1}$]')
ylabel('Shear modulus $[{\rm pN}\cdot \mu{\rm m^{-1}}]$')
exportgraphics(crslnkStiffShear_famOfan,fullfile(saveDir,figSubdir, ...
    ['crslnkStiffShear_famOfan',filetype]),'ContentType','vector')

%% SUPP: FILAMENT BENDING RIGIDITY %%

[stdevMod_kbend,meanMod_kbend] = std(moduli_s_kbend,0,3);

kbendShear_famOfkbend = figure(80); clf;
set(kbendShear_famOfkbend,'units','centimeters','Position',[1,1,20,16])
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
lg.FontSize = fszInset;
title(lg,{"Filament bending"; "rigidity [${\rm pN}\cdot \mu {\rm m}^2$]"})
xlabel('Astral number $a_n$')
xticks(kbendParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN} \cdot \mu{\rm m^{-1}}$]')
exportgraphics(kbendShear_famOfkbend,fullfile(saveDir,figSubdir, ...
    ['kbendShear_famOfkbend',filetype]),'ContentType','vector')

kbendShear_famOfan = figure(79); clf;
set(kbendShear_famOfan,'units','centimeters','Position',[1,1,20,16])
an_to_show = [1,2,4,5,7];     % indices of an = 1,2,4,8
cmap = colormap(parula(length(an_to_show)+1));
hold on
for idx = 1:length(an_to_show)
    errorbar(kbendParams.kbend_list,meanMod_kbend(an_to_show(idx),:),...
        (z/sqrt(kbendParams.numRep)) * ...
        stdevMod_kbend(an_to_show(idx),:),'Color', ...
        cmap(idx,:),'LineWidth',1.5)
end
set(gca,'FontName',fname,'FontSize',fsz2)
lg = legend(string(kbendParams.nFilPerAsterList(an_to_show)),'Location', ...
    'south','FontSize',fszInset,'NumColumns',5);
title(lg,"Astral number $a_n$")
xscale('log')
xlim([10^(-3),0.2])
xlabel('Filament bending rigidity $[{\rm pN}\cdot \mu {\rm m}^2]$')
ylim([0,inf])
ylabel('Shear modulus $[{\rm pN}\cdot \mu{\rm m^{-1}}]$')
exportgraphics(kbendShear_famOfan,fullfile(saveDir,figSubdir, ...
    ['kbendShear_famOfan',filetype]),'ContentType','vector')


%% SUPP: ANGULAR STIFFNESS AT ASTRAL CENTERS %%

[stdevMod_krot,meanMod_krot] = std(moduli_s_krot,0,3);

krotShear_famOfkrot = figure(78); clf;
set(krotShear_famOfkrot,'units','centimeters','Position',[1,1,20,16])
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
lg.FontSize = fszInset;
title(lg,{"$\:$ Angular stiffness at $\:$"; ...
    "$\:$ centers [${\rm pN}\cdot {\mu}{\rm m^{-1}}$] $\:$"})
xlabel('Astral number $a_n$')
xticks(krotParams.nFilPerAsterList)
ylabel('Shear modulus [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
exportgraphics(krotShear_famOfkrot,fullfile(saveDir,figSubdir, ...
    ['krotShear_famOfkrot',filetype]),'ContentType','vector')

krotShear_famOfan = figure(77); clf;
set(krotShear_famOfan,'units','centimeters','Position',[1,1,20,16])
an_to_show = [1,2,4,5,7];     % indices of an = 1,2,4,5,8
cmap = colormap(parula(length(an_to_show)+1));
hold on
for idx = 1:length(an_to_show)
    errorbar(krotParams.krot_list,meanMod_krot(an_to_show(idx),:),...
        (z/sqrt(krotParams.numRep)) * ...
        stdevMod_krot(an_to_show(idx),:),'Color', ...
        cmap(idx,:),'LineWidth',1.5)
end
set(gca,'FontName',fname,'FontSize',fsz2)
lg = legend(string(krotParams.nFilPerAsterList(an_to_show)),'Location', ...
    'south','FontSize',fszInset,'NumColumns',5);
title(lg,"Astral number $a_n$")
xlim([-10,410])
xlabel('Angular stiffness at centers [${\rm pN}\cdot {\mu}{\rm m^{-1}}$]')
ylim([0,inf])
ylabel('Shear modulus $[{\rm pN}\cdot \mu{\rm m^{-1}}]$')
exportgraphics(krotShear_famOfan,fullfile(saveDir,figSubdir, ...
    ['krotShear_famOfan',filetype]),'ContentType','vector')

%% SUPP: PERCOLATION PROBABILITIES & FITS - RANDOM ANGLE %%

an_to_show = [1,2,4,8];
nNetTypes = length(an_to_show);
nSizes = length(percParams_rnd_famOfD.Dlist);

% percolation data generated at l=1, D=10
% simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
% identical percolation scenario (same N)
lengthScaleFactor = sweep_s_rnd.len_fil / percParams_rnd_famOfD.l;
densityScaleFactor = 1/lengthScaleFactor;

% need to refresh fits to account for 10x density scaling
smsplFits_rnd_famOfD = cell(percParams_rnd_famOfD.numNetTypes, ...
    percParams_rnd_famOfD.numPercTypes,percParams_rnd_famOfD.numD);
for idx = 1:percParams_rnd_famOfD.numNetTypes
    netLabel = sprintf('an%02d',percParams_rnd_famOfD.allAstralNums(idx));
    for jdx = 1:percParams_rnd_famOfD.numPercTypes
        for kdx = 1:percParams_rnd_famOfD.numD
            if isfield(allCurves{kdx},netLabel)
                thisRawData = percDataUsed_rnd_famOfD{idx,jdx,kdx};
                x = thisRawData(:,1) * densityScaleFactor;
                y = thisRawData(:,2);
                smsplFits_rnd_famOfD{idx,jdx,kdx} = fit(x,y, ...
                    'smoothingspline');
            end
        end
    end
end

spanDataFits_rnd = figure(76); clf;
set(spanDataFits_rnd,'units','centimeters','Position',[1,1,35,12])
spanDataFits_rnd = tiledlayout(spanDataFits_rnd,1,4);
legendLabels = cell(2*nSizes,1);
for idx = 1:nNetTypes
    nexttile(idx)
    hold on
    for kdx = 1:nSizes
        thisData = percDataUsed_rnd_famOfD{an_to_show(idx),1,kdx};
        l = plot(smsplFits_rnd_famOfD{an_to_show(idx),1,kdx}, ...
            densityScaleFactor*thisData(:,1),thisData(:,2),'*');
        [l.Color] = deal(colors(kdx,:));
        [l.LineWidth] = deal(0.75);
        legendLabels{2*kdx - 1} = ''; % suppress data entry in legend
        legendLabels{2*kdx} = sprintf('%2.1f', ...
            percParams_rnd_famOfD.Dlist(kdx)*lengthScaleFactor);
    end
    yline(0.5,'--','LineWidth',1)
    hold off
    title(sprintf('$a_n = %i$',an_to_show(idx)))
    xlabel('')
    if idx == 1
        xlim([10,inf])
        xticks(10.^(1:2))
    end
    ylabel('')
    xscale('log')
    set(gca,'FontSize',fsz2,'FontName',fname)
    if idx == nNetTypes
        lg = legend(legendLabels,'FontSize',fszInset);
        lg.Layout.Tile = 'east';
        title(lg,{"System size"; "$s \; [{\rm \mu m}]$"},'FontSize',fszInset)
        xlabel(spanDataFits_rnd,'Density $\rho$ $[{\rm \mu m^{-1}}]$', ...
            'Interpreter','latex','FontSize',fsz2)
        ylabel(spanDataFits_rnd,'$\Pr\{$Spanning comp.$\}$', ...
            'Interpreter','latex','FontSize',fsz2)
    else
        legend('off')
    end
end
exportgraphics(spanDataFits_rnd,fullfile(saveDir,figSubdir, ...
    ['spanDataFits_rnd',filetype]),'ContentType','vector')

connDataFits_rnd = figure(75); clf;
set(connDataFits_rnd,'units','centimeters','Position',[1,1,35,12])
connDataFits_rnd = tiledlayout(connDataFits_rnd,1,4);
legendLabels = cell(2*nSizes,1);
for idx = 1:nNetTypes
    nexttile(idx)
    hold on
    for kdx = 1:nSizes
        thisData = percDataUsed_rnd_famOfD{an_to_show(idx),5,kdx};
        l = plot(smsplFits_rnd_famOfD{an_to_show(idx),5,kdx}, ...
            densityScaleFactor*thisData(:,1),thisData(:,2),'*');
        [l.Color] = deal(colors(kdx,:));
        [l.LineWidth] = deal(0.75);
        legendLabels{2*kdx - 1} = ''; % suppress data entry in legend
        legendLabels{2*kdx} = sprintf('%2.1f', ...
            percParams_rnd_famOfD.Dlist(kdx)*lengthScaleFactor);
    end
    yline(0.5,'--','LineWidth',1)
    hold off
    set(gca,'FontSize',fsz2,'FontName',fname)
    title(sprintf('$a_n = %i$',an_to_show(idx)))
    xlabel('')
    if idx == 1
        xticks(10.^(2:4))
        set(gca,'FontSize',fszInset)
    elseif idx == 2 || idx == 3
        xlim([10,inf])
        xticks(10.^(1:2))
    end
    ylabel('')
    xscale('log')
    if idx == nNetTypes
        lg = legend(legendLabels,'FontSize',fszInset);
        lg.Layout.Tile = 'east';
        title(lg,{"System size"; "$s \; [{\rm \mu m}]$"},'FontSize',fszInset)
        xlabel(connDataFits_rnd,'Density $\rho$ $[{\rm \mu m^{-1}}]$', ...
            'Interpreter','latex','FontSize',fsz2)
        ylabel(connDataFits_rnd,'$\Pr\{$Unique conn. comp.$\}$', ...
            'Interpreter','latex','FontSize',fsz2)
    else
        legend('off')
    end
end
exportgraphics(connDataFits_rnd,fullfile(saveDir,figSubdir, ...
    ['connDataFits_rnd',filetype]),'ContentType','vector')

%% test: MODULI BY PERCOLATION STATUS %%

boxPlotModuli_s = moduli_s_hiN;
boxPlotModuli_t = moduli_t_hiN;
for netIdx = 1:modSampling.nNetTypes
    for repIdx = 1:modSampling.numRep
        % for networks that failed spanning check, dig up what modulus
        % would have been assigned from linefit records
        % if negative, clip to 0
        % shear
        if ~percStatus_s_hiN(netIdx,repIdx,1)
            thisLine = shearLineFits_hiN.(networkLabel( ...
                modSampling.nFilPerAsterList(netIdx))){repIdx};
            boxPlotModuli_s(netIdx,repIdx) = max([1/thisLine.p1,0]);
        end
        % tension
        if ~percStatus_t_hiN(netIdx,repIdx,1)
            thisLine = tensLineFits_hiN.(networkLabel( ...
                modSampling.nFilPerAsterList(netIdx))){repIdx};
            boxPlotModuli_t(netIdx,repIdx) = max([1/thisLine.p1,0]);
        end
    end
end

% define groupings for box plots
notSpan_s = ~percStatus_s_hiN(:,:,1);
spanOnly_s = logical(percStatus_s_hiN(:,:,1) .* ~percStatus_s_hiN(:,:,2));
spanAndConn_s = logical(prod(percStatus_s_hiN,3));
notSpan_t = ~percStatus_t_hiN(:,:,1);
spanOnly_t = logical(percStatus_t_hiN(:,:,1) .* ~percStatus_t_hiN(:,:,2));
spanAndConn_t = logical(prod(percStatus_t_hiN,3));

groupNames = {'NS','SO','S&C'};

groups_s = cell(modSampling.nNetTypes,modSampling.numRep);
[groups_s{notSpan_s}] = deal(groupNames{1});     % group 1: not spanning
[groups_s{spanOnly_s}] = deal(groupNames{2});    % group 2: spanning only
[groups_s{spanAndConn_s}] = deal(groupNames{3}); % group 3: spanning and fully connected
groups_s = transpose(groups_s);

groups_t = cell(modSampling.nNetTypes,modSampling.numRep);
[groups_t{notSpan_t}] = deal(groupNames{1});     % group 1: not spanning
[groups_t{spanOnly_t}] = deal(groupNames{2});    % group 2: spanning only
[groups_t{spanAndConn_t}] = deal(groupNames{3}); % group 3: spanning and fully connected
groups_t = transpose(groups_t);

modByPercStat_s = figure(50); clf;
set(modByPercStat_s,'units','centimeters','Position',[1,1,15,22])
modByPercStat_s = tiledlayout(4,2,'TileSpacing','loose');
for netIdx = 1:modSampling.nNetTypes
    nexttile(netIdx)
    thisTable = table(groups_s(:,netIdx),boxPlotModuli_s(netIdx,:)', ...
        'VariableNames',{'percStatus','modulus'});
    grps = categorical(thisTable.percStatus,groupNames);
    % original idea: box plots
    % boxplot(thisTable.modulus,grps)
    % alternate idea: swarmchart
    swarmchart(grps,thisTable.modulus)
    title(sprintf('$a_n=%d$',modSampling.nFilPerAsterList(netIdx)))
end
title(modByPercStat_s,'Shear moduli','Interpreter','latex')
exportgraphics(modByPercStat_s,fullfile(saveDir,figSubdir, ...
    ['modByPercStat_s_swrm',filetype]),'ContentType','vector')

modByPercStat_t = figure(87); clf;
set(modByPercStat_t,'units','centimeters','Position',[1,1,15,22])
modByPercStat_t = tiledlayout(4,2,'TileSpacing','loose');
for netIdx = 1:modSampling.nNetTypes
    nexttile(netIdx)
    thisTable = table(groups_t(:,netIdx),boxPlotModuli_t(netIdx,:)', ...
        'VariableNames',{'percStatus','modulus'});
    grps = categorical(thisTable.percStatus,groupNames);
    % original idea: box plots
    % boxplot(thisTable.modulus,grps)
    % alternate idea: swarmchart
    swarmchart(grps,thisTable.modulus)
    title(sprintf('$a_n=%d$',modSampling.nFilPerAsterList(netIdx)))
end
title(modByPercStat_t,"Young's moduli",'Interpreter','latex')
exportgraphics(modByPercStat_t,fullfile(saveDir,figSubdir, ...
    ['modByPercStat_t_swrm',filetype]),'ContentType','vector')


