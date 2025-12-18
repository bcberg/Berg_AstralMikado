% rf_shear_aster_actin_physkT.m
% Brady Berg, 10/24/2025
clear
close all
format short
format compact

%% Parse report files

nameprefix = 'rf_shear_aster_actin_physkT251113';

% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data/';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data\';
dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';
matFile = [nameprefix,'.mat'];

newReports = false;
if newReports
    filtering.TF = true;
    filtering.catchDrops = false;   % at higher k_B*T, lowest applied force
    % does not always significantly deform the network
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/" + ...
        "rf_shear_aster_actin_physkT/";

    % ensure the following parameters match those used in 
    % rf_shear_aster_actin_physkT.cym.tpl and its driver (for numRep)
    dens = 75;
    len_fil = 0.1;
    D = 1;  % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = (1:24)';
    nNetTypes = length(nFilPerAsterList);
    initFrames = 5:14;  % samples from 41s, 42s, ..., 50s
    numInit = length(initFrames);
    finalFrames = 39:48;    % samples from 291s, 292s, ..., 300s
    numFinal = length(finalFrames);
    dataStartLine = 6;
    ptsPerFiber = 6;
    forceVals = 0.25 * (0:5)';
    nForces = length(forceVals);
    numRep = 20;

    runsPerRep = nNetTypes * nForces;
    shearModuli = zeros(nNetTypes,numRep);
    percStatus = false(nNetTypes,numRep,2);
    for idx = 1:numel(shearModuli)
        [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
        nFilPerAster = nFilPerAsterList(netIdx);
        netLabel = sprintf('an%02i',nFilPerAster);

        if repIdx == 1
            initCoMs.(netLabel) = zeros(nForces,2,numRep);
            allInitCoMs.(netLabel) = zeros(nForces,2,numInit,numRep);
            finalCoMs.(netLabel) = zeros(nForces,2,numRep);
            allFinalCoMs.(netLabel) = zeros(nForces,2,numFinal,numRep);
            disps.(netLabel) = zeros(nForces,2,numRep);
            lineFits.(netLabel) = cell(numRep,1);
            fitStats.(netLabel) = cell(numRep,1);
        end

        cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
            'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);
        unshiftedRunVals = (netIdx - 1)*nForces:(netIdx * nForces - 1);
        theseRunVals = (repIdx - 1) * runsPerRep + unshiftedRunVals;
        % check network percolation status
        [percTF,~,~] = checkLinks(dir,theseRunVals(1),cytoparams,4);
        filtering.spanCheck = percTF(1);
        percStatus(netIdx,repIdx,:) = reshape(percTF,[1,1,2]);
        % analyze network position data
        [thisInitCoM,thisFinalCoM,thisAllInit,thisAllFinal] = ...
            parseForceSweep_manySamp(dir,theseRunVals,initFrames, ...
            finalFrames,cytoparams);
        initCoMs.(netLabel)(:,:,repIdx) = thisInitCoM;
        allInitCoMs.(netLabel)(:,:,:,repIdx) = thisAllInit;
        finalCoMs.(netLabel)(:,:,repIdx) = thisFinalCoM;
        allFinalCoMs.(netLabel)(:,:,:,repIdx) = thisAllFinal;
        disps.(netLabel)(:,:,repIdx) = thisFinalCoM - thisInitCoM;
        % estimate a modulus
        [thisModulus, linefit, gof] = getModulus(thisInitCoM,thisFinalCoM,...
            forceVals,"shear",filtering);
        shearModuli(netIdx,repIdx) = thisModulus;
        lineFits.(netLabel){repIdx} = linefit;
        fitStats.(netLabel){repIdx} = gof;
    end
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
        'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
        'nNetTypes','initCoMs','allInitCoMs','finalCoMs','allFinalCoMs',...
        'forceVals','disps','dens','percStatus','filtering')
    clearvars -except nameprefix matFile saveDir dataSubdir figSubdir
end

load(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
    'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
    'nNetTypes','initCoMs','allInitCoMs','finalCoMs','allFinalCoMs',...
    'forceVals','disps','dens','percStatus','filtering')

%% Refit data

refit = false;
if refit
    clearvars -except saveDir dataSubdir figSubdir matFile nameprefix ...
        checkDispData
    filtering.TF = true;
    filtering.catchDrops = false;   % at higher k_B*T, lowest applied force
    % does not always significantly deform the network
    load(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
        'nFilPerAsterList','numRep',...
        'nNetTypes','initCoMs','allInitCoMs','finalCoMs','allFinalCoMs',...
        'forceVals','disps','dens','percStatus')
    shearModuli = zeros(nNetTypes,numRep);
    for idx = 1:numel(shearModuli)
        [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
        nFilPerAster = nFilPerAsterList(netIdx);
        netLabel = sprintf('an%02d',nFilPerAster);

        if repIdx == 1
            lineFits.(netLabel) = cell(numRep,1);
            fitStats.(netLabel) = cell(numRep,1);
        end
        filtering.spanCheck = percStatus(netIdx,repIdx,1);
        [thisModulus,linefit,gof] = getModulus(initCoMs.(netLabel)(:,:,repIdx), ...
            finalCoMs.(netLabel)(:,:,repIdx),forceVals,"shear", ...
            filtering);
        shearModuli(netIdx,repIdx) = thisModulus;
        lineFits.(netLabel){repIdx} = linefit;
        fitStats.(netLabel){repIdx} = gof;
    end
    save(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
        'nFilPerAsterList','numRep','nNetTypes','initCoMs','allInitCoMs', ...
        'finalCoMs','allFinalCoMs','forceVals','disps','dens', ...
        'shearModuli','lineFits','fitStats','percStatus','filtering')
end

%% Summary plots

[stdevs,means] = std(shearModuli,0,2); % compute summary stats across cols
z = norminv(0.975);
fig1 = figure(1); clf;
hold on
for idx = 1:nNetTypes
    plot(nFilPerAsterList(idx),shearModuli(idx,:),'*k','MarkerSize',4)
end
errorbar(nFilPerAsterList,means,z*stdevs/sqrt(numRep),'-ob','LineWidth',1.5)
hold off
xlabel('Astral number')
ylabel('2D shear modulus [pN/{\mu}m]')
% ylim([0, inf])
exportgraphics(fig1,fullfile(saveDir,figSubdir,[nameprefix,'-moduli.png']), ...
    'Resolution',300)

% summarizing R^2 stats
Rsquares = zeros(nNetTypes,numRep);
for idx = 1:numel(Rsquares)
    [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
    networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
    thisGOF = fitStats.(networkLabel){repIdx};
    Rsquares(netIdx,repIdx) = thisGOF.rsquare;
end
fig2 = figure(2); clf;
hold on
for repIdx = 1:numRep
    plot(nFilPerAsterList,Rsquares(:,repIdx),'DisplayName',...
        sprintf('t%02i',repIdx),'LineWidth',1)
end
hold off
xlabel('Number of filaments per aster')
ylabel('R^2 value for linear fit')
legend('Location','eastoutside','NumColumns',2)

%% Displacement vs. force plots

plotTiling = [5,4];
dispPlots = figure(3);
set(dispPlots,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
    'defaultTextInterpreter','tex')
dispPlots = tiledlayout(dispPlots,plotTiling(1),plotTiling(2), ...
    'TileSpacing','compact');
numInit = size(allInitCoMs.an01,3);
numFinal = size(allFinalCoMs.an01,3);
z = norminv(0.975);
for netIdx = 1:nNetTypes
    astralNum = nFilPerAsterList(netIdx);
    netLabel = sprintf('an%02d',astralNum);
    for repIdx = 1:numRep
        % approximate displacement distribution as (final position samples)
        % - (average initial position)
        theseFinalPos = squeeze(allFinalCoMs.(netLabel)(:,1,:,repIdx));
        % after squeeze, shaped as theseFinalPos(forceVal,timeSample)
        thisInitPos = repmat(initCoMs.(netLabel)(:,1,repIdx),[1,numFinal]);
        finalDispSamples = theseFinalPos - thisInitPos;
        % compute stdev across timesamples
        dispStdev = std(finalDispSamples,0,2);
        dispMean = disps.(netLabel)(:,1,repIdx);
        nexttile(repIdx)
        hold on
        l1 = plot(forceVals,finalDispSamples,'xk');
        if shearModuli(netIdx,repIdx) == 0
            l2 = errorbar(forceVals,dispMean,dispStdev * z / sqrt(numFinal), ...
                '-o','LineWidth',1,'Color',[1,0,0],'DisplayName','disps');
            if ~percStatus(netIdx,repIdx,1) % if no spanning comp.
                text(forceVals(3),0.1*max(dispMean),'no span')
            end
        else
            l2 = errorbar(forceVals,dispMean,dispStdev * z / sqrt(numFinal), ...
                '-o','LineWidth',1,'Color',[0,0,1],'DisplayName','disps');
            l3 = plot(forceVals,forceVals / shearModuli(netIdx,repIdx),'--k', ...
                'LineWidth',1,'DisplayName','linefit');
        end
        hold off
        title(sprintf('rep%02d',repIdx))
    end
    xlabel(dispPlots,'Shear force [pN]')
    ylabel(dispPlots,'x-Displacement [{\mu}m]')
    title(dispPlots,netLabel)
    if netIdx == 1
        % export first page
        if ishandle(l3)
            lg = legend([l2,l3],'NumColumns',2);
        else
            lg = legend(l2,'NumColumns',1);
        end
        lg.Layout.Tile = 'south';
        exportgraphics(dispPlots,fullfile(saveDir,figSubdir, ...
            [nameprefix,'-linefits.pdf']),'Resolution',300)
        dispPlots = clf(dispPlots);
        set(dispPlots,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
            'defaultTextInterpreter','tex')
        dispPlots = tiledlayout(dispPlots,plotTiling(1),plotTiling(2), ...
            'TileSpacing','compact');
    else
        % export following pages
        if ishandle(l3)
            lg = legend([l2,l3],'NumColumns',2);
        else
            lg = legend(l2,'NumColumns',1);
        end
        lg.Layout.Tile = 'south';
        exportgraphics(dispPlots,fullfile(saveDir,figSubdir, ...
            [nameprefix,'-linefits.pdf']),'Resolution',300,'Append',true)
        dispPlots = clf(dispPlots);
        set(dispPlots,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
            'defaultTextInterpreter','tex')
        dispPlots = tiledlayout(dispPlots,plotTiling(1),plotTiling(2), ...
            'TileSpacing','compact');
    end
end