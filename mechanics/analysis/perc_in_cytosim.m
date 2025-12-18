% perc_in_cytosim.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "perc_in_cytosim251217";
% Ubuntu filepath
% saveDir = '~/Documents/astral-mikado-data';
% Windows filepath
saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';
dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';

newReports = false;
if newReports
    % Ubuntu
    % dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/" + ...
    %     "perc_in_cytosim/";
    % Windows
    dir = "C:\Users\bcber\Documents\AstralMikadoCYM\runs\" + ...
        "perc_in_cytosim\";
    
    load(fullfile(saveDir,dataSubdir,'percProbs_l01_D10.mat'), ...
    'l','D','actualDensities','percProbs')
    percParamsFromMatlab = struct('l',l,'D',D);
    clear('l','D')

    % ensure the following parameters match those used in
    % perc_in_cytosim.cym.tpl and its driver (for numRep)
    len_fil = 0.1;
    D = 1; % side length of square region
    astralNum = 2;

    % 4 selected parameter sets
    % (viscosity, kT, bind_rate, r_bind, crslnk_per_fil, crslnk_period)
    % units: (pN*s/um^2, pN*um, s^{-1}, um, n/a, s)
    paramSets = {[0.5, 2.1e-5, 1, 1e-3, 30, 50];
        [0.5, 4.2e-3, 1, 1e-3, 30, 50];
        [1e3, 2.1e-5, 1e2, 1e-4, 30, 5];
        [1e3, 2.1e-5, 1e2, 2.5e-5, 120, 5]};
    numParamSets = length(paramSets);

    % percolation data generated at l=1, D=10
    densFromMatlabSamp = actualDensities{astralNum};
    % simulation params were l=0.1, D=1 and 10*rho = 10*(N*l/D^2) to obtain
    % identical percolation scenario (same N)
    densityScaleFactor = percParamsFromMatlab.l / len_fil;
    densInCytosim = densityScaleFactor * densFromMatlabSamp(1:35);
    numDens = length(densInCytosim);
    nFil = densInCytosim * D^2 / len_fil;
    dataStartLine = [];  % irrelevant, not inspecting position data
    ptsPerFiber = [];   % irrelevant, not inspecting position data
    numRep = 100;

    % runsPerRep = numDens * numBindRates;  % runs before 251213
    % percStatus = false(numDens,numRep,numBindRates,2);
    % percProbsCYM = zeros(numDens,numBindRates,2);
    runsPerRep = numDens * numParamSets;
    percStatus = false(numDens,numParamSets,numRep,2);
    percProbsCYM = zeros(numDens,numParamSets,2);
    % "non-sweep" parameters set all at once at head of file, then density
    for kdx = 1:numRep
        repShift = (kdx-1) * runsPerRep;
        for jdx = 1:numParamSets
            paramShift = (jdx - 1) * numDens;
            for idx = 1:numDens
                thisRun = (idx-1) + paramShift + repShift;
                cytoparams = struct('nFil',nFil(idx),'nFilPerAster', ...
                    astralNum,'dataStartLine',dataStartLine, ...
                    'ptsPerFiber',ptsPerFiber);
                [percTF,~,~] = checkLinks(dir,thisRun,cytoparams,5);
                percStatus(idx,jdx,kdx,:) = reshape(percTF,[1,1,1,2]);
            end
        end
    end
    % spanning percolation probabilities
    percProbsCYM(:,:,1) = sum(percStatus(:,:,:,1),3) / numRep;
    % connectivity percolation probabilities
    percProbsCYM(:,:,2) = sum(percStatus(:,:,:,2),3) / numRep;

    save(fullfile(saveDir,dataSubdir,simName + ".mat"), ...
        'len_fil','D','astralNum','paramSets','numParamSets', ...
        'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
        'percStatus','percProbsCYM')
    clearvars -except saveDir dataSubdir figSubdir simName
end

load(fullfile(saveDir,dataSubdir,simName + ".mat"), ...
    'len_fil','D','astralNum','paramSets','numParamSets', ...
    'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
    'percStatus','percProbsCYM')


%% Summary plots

% option to also load smsplFits and to50 from percFits_nonparam_tzOnly.mat
% for comparison; would need to rescale densities
load(fullfile(saveDir,dataSubdir,'percProbs_l01_D10.mat'), ...
    'curves')

matlabCurve = curves.(sprintf('an%02d',astralNum));
matlabCurve = matlabCurve(:,1:35);
curveLabels = {"Matlab MC";
    "Cold";
    "Physiological";
    "Viscous set 1";
    "Viscous set 2"};

percCYMplots = figure(1); clf;
colors = colororder("gem12");
set(percCYMplots,'units','centimeters','Position',[1,1,25,20], ...
    'defaultLineLineWidth',1.5,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex','defaultLegendInterpreter', ...
    'latex')
percCYMplots = tiledlayout(2,1,"TileSpacing","compact");
nexttile(1)     % spanning perc. probs
hold on
plot(densityScaleFactor * matlabCurve(1,:),matlabCurve(2,:),'-ok', ...
    'DisplayName',curveLabels{1})
for kdx = 1:numParamSets
    % plot(densInCytosim,percProbsCYM(:,kdx,1),'-x','Color',colors(kdx,:), ...
    %     'DisplayName',sprintf('Bind=%.1g s^{-1}',crslnkBindRate(kdx)))
    % plot(densInCytosim,percProbsCYM(:,kdx,1),'-x','Color',colors(kdx,:), ...
    %     'DisplayName',sprintf('numCrslnk=%i * nFil, rBind=%.3g', ...
    %     crslnkPerFil(kdx),rbind(kdx)))
    plot(densInCytosim,percProbsCYM(:,kdx,1),'-x','Color',colors(kdx,:), ...
        'DisplayName',curveLabels{1+kdx})
end
hold off
xscale('log')
ylabel({'Spanning';'probability'})
set(gca,'FontName','CMU Serif','FontSize',22)
lg = legend();
lg.Layout.Tile = 'east';
nexttile(2)     % connectivity perc. probs
hold on
plot(densityScaleFactor * matlabCurve(1,:),matlabCurve(6,:),'-ok')
for kdx = 1:numParamSets
    plot(densInCytosim,percProbsCYM(:,kdx,2),'-x','Color',colors(kdx,:))
end
hold off
xscale('log')
ylabel({'Connectivity'; 'probability'})
set(gca,'FontName','CMU Serif','FontSize',22)
xlabel(percCYMplots,'Density $\rho\; [\mu {\rm m^{-1}}]$', ...
    'Interpreter','latex','FontName','CMU Serif','FontSize',22)
% title(percCYMplots,sprintf('a_n = %i',astralNum))
exportgraphics(percCYMplots,fullfile(saveDir,figSubdir,simName + ".pdf"), ...
    'ContentType','vector')

%% Synthesizing selected data

% clear;
% % Ubuntu filepath
% saveDir = '~/Documents/astral-mikado-data';
% % Windows filepath
% % saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';
% dataSubdir = 'mat_files';
% figSubdir = 'exploratory_figures';
% 
% % bindRate = [1,1e2], has r_bind = 1e-4, crosslinker amt = 30*nFil
% load(fullfile(saveDir,dataSubdir,"perc_in_cytosim251211-3" + ".mat"), ...
%     'len_fil','D','astralNum','crslnkBindRate','numCrslnkVals', ...
%     'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
%     'percProbsCYM')
% params251211_3 = struct('len_fil',len_fil,'D',D,'astralNum',astralNum, ...
%     'crslnkBindRate',crslnkBindRate,'numBindRates',numCrslnkVals, ...
%     'densityScaleFactor',densityScaleFactor,'densInCytosim',densInCytosim, ...
%     'numDens',numDens,'nFil',nFil,'numRep',numRep);
% percProbsCYM251211_3 = percProbsCYM;
% clear('len_fil','D','astralNum','crslnkBindRate','numCrslnkVals', ...
%     'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
%     'percProbsCYM')
% 
% % has bindRate = 1e2, r_bind = 2.5e-5, crosslinker amt = 120*nFil
% load(fullfile(saveDir,dataSubdir,"perc_in_cytosim251212" + ".mat"), ...
%     'len_fil','D','astralNum','crslnkBindRate','numCrslnkVals', ...
%     'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
%     'percProbsCYM')
% params251212 = struct('len_fil',len_fil,'D',D,'astralNum',astralNum, ...
%     'crslnkBindRate',crslnkBindRate,'numBindRates',numCrslnkVals, ...
%     'densityScaleFactor',densityScaleFactor,'densInCytosim',densInCytosim, ...
%     'numDens',numDens,'nFil',nFil,'numRep',numRep);
% percProbsCYM251212 = percProbsCYM;
% clear('len_fil','D','astralNum','crslnkBindRate','numCrslnkVals', ...
%     'densityScaleFactor','densInCytosim','numDens','nFil','numRep', ...
%     'percProbsCYM')
% 
% % Matlab MC data
% load(fullfile(saveDir,dataSubdir,'percFits_nonparam_tzOnly.mat'), ...
%     'percDataUsed')
% 
% percCYMplotsv2 = figure(2); clf;
% colors = colororder("gem12");
% set(percCYMplotsv2,'units','centimeters','Position',[1,1,15,20], ...
%     'defaultLineLineWidth',1)
% percCYMplotsv2 = tiledlayout(2,1,"TileSpacing","compact");
% nexttile(1)     % spanning perc. probs
% hold on
% plot(params251211_3.densityScaleFactor * ...
%     percDataUsed{params251211_3.astralNum,1}(:,1), ...
%     percDataUsed{params251211_3.astralNum,1}(:,2),'-ok', ...
%     'DisplayName','Matlab MC')
% plot(params251211_3.densInCytosim,percProbsCYM251211_3(:,2,1),'Color', ...
%     colors(1,:),'DisplayName','r_{bind}=1e-4, numCrslks=30*nFil')
% plot(params251212.densInCytosim,percProbsCYM251212(:,1,1),'Color', ...
%     colors(2,:),'DisplayName','r_{bind}=2.5e-5, numCrslks=120*nFil')
% hold off
% xscale('log')
% ylabel('Spanning perc. probability')
% lg = legend;
% lg.Layout.Tile = 'south';
% nexttile(2)     % connectivity perc. probs
% hold on
% plot(params251211_3.densityScaleFactor * ...
%     percDataUsed{params251211_3.astralNum,5}(:,1), ...
%     percDataUsed{params251211_3.astralNum,5}(:,2),'-ok')
% plot(params251211_3.densInCytosim,percProbsCYM251211_3(:,2,2),'Color', ...
%     colors(1,:)) % 'DisplayName','r_{bind}=1e-4, numCrslks=30*nFil'
% plot(params251212.densInCytosim,percProbsCYM251212(:,1,2),'Color', ...
%     colors(2,:)) %'DisplayName','r_{bind}=2.5e-5, numCrslks=120*nFil'
% hold off
% xscale('log')
% ylabel('Connectivity perc. probability')
% 
% xlabel(percCYMplotsv2,'Density \rho [{\mu}m^{-1}]')
% % title(percCYMplots,sprintf('a_n = %i',astralNum))
% title(percCYMplotsv2,[sprintf('a_n = %i',params251211_3.astralNum), ...
%     ', visc. = 1e3, Bind=1e2'])
% exportgraphics(percCYMplotsv2,fullfile(saveDir,figSubdir, ...
%     "perc_in_cytosim_rbindexploration" + ".pdf"), 'ContentType','vector')