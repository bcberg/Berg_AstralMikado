% modulus_distr.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

nameprefix = 'mod_distr251020';

% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data/';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\AstralMikadoCYM\data';
dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';
matFile = [nameprefix,'.mat'];

newReports = false;
if newReports
    filtering_s.TF = true;
    filtering_t.TF = true;
    rundir = "/home/bcberg/Documents/AstralMikadoCYM/runs/";
    sheardir = fullfile(rundir,"modulus_distr_shear/");
    tensdir = fullfile(rundir,"modulus_distr_tens/");
    
    % ensure the following parameters match those used in
    % modulus_distr_****.cym.tpl and its driver (for numRep)
    dens = 75;  % N * len_fil / A
    len_fil = 0.1;
    D = 1; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = transpose([1:3,4:4:16]);
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    shearForces = transpose((0:5) * 0.05);
    tensForces = transpose((0:5) * 0.2);
    nShearForces = length(shearForces);
    nTensForces = length(tensForces);
    numRep = 100;
    
    runsPerRep = nNetTypes * nShearForces;
    shearModuli = zeros(nNetTypes,numRep);
    tensModuli = zeros(nNetTypes,numRep);
    percStatus_s = false(nNetTypes,numRep,2);
    percStatus_t = false(nNetTypes,numRep,2);
    for idx = 1:(nNetTypes*numRep)
        [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
        nFilPerAster = nFilPerAsterList(netIdx);
        networkLabel = sprintf('an%02i',nFilPerAster);

        if repIdx == 1
            shearInitCoMs.(networkLabel) = zeros(nShearForces,2,numRep);
            shearFinalCoMs.(networkLabel) = zeros(nShearForces,2,numRep);
            shearDisps.(networkLabel) = zeros(nShearForces,2,numRep);
            shearLineFits.(networkLabel) = cell(numRep,1);
            shearFitStats.(networkLabel) = cell(numRep,1);
            tensInitCoMs.(networkLabel) = zeros(nTensForces,2,numRep);
            tensFinalCoMs.(networkLabel) = zeros(nTensForces,2,numRep);
            tensDisps.(networkLabel) = zeros(nTensForces,2,numRep);
            tensLineFits.(networkLabel) = cell(numRep,1);
            tensFitStats.(networkLabel) = cell(numRep,1);
        end

        cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
            'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);

        unshiftedRunVals_s = (netIdx - 1)*nShearForces:(netIdx * nShearForces - 1);
        theseRunVals_s = (repIdx - 1) * runsPerRep + unshiftedRunVals_s;
        % check network percolation status
        [percTF_s,~,~] = checkLinks(sheardir,theseRunVals_s(1),cytoparams);
        filtering_s.spanCheck = percTF_s(1);
        percStatus_s(netIdx,repIdx,:) = reshape(percTF_s,[1,1,2]);
        % analyze network position data
        [thisInitCoM_s,thisFinalCoM_s] = parseForceSweep(sheardir, ...
            theseRunVals_s,cytoparams);
        shearInitCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM_s;
        shearFinalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM_s;
        shearDisps.(networkLabel)(:,:,repIdx) = thisFinalCoM_s - thisInitCoM_s;
        % estimate a modulus
        [thisShearModulus, linefit_s, gof_s] = getModulus(thisInitCoM_s, ...
            thisFinalCoM_s,shearForces,"shear",filtering_s);
        shearModuli(netIdx,repIdx) = thisShearModulus;
        shearLineFits.(networkLabel){repIdx} = linefit_s;
        shearFitStats.(networkLabel){repIdx} = gof_s;
        

        unshiftedRunVals_t = (netIdx - 1)*nTensForces:(netIdx * nTensForces - 1);
        theseRunVals_t = (repIdx - 1) * runsPerRep + unshiftedRunVals_t;
        % check network percolation status
        [percTF_t,~,~] = checkLinks(tensdir,theseRunVals_t(1),cytoparams);
        filtering_t.spanCheck = percTF_t(1);
        percStatus_t(netIdx,repIdx,:) = reshape(percTF_t,[1,1,2]);
        % analyze network position data
        [thisInitCoM_t,thisFinalCoM_t] = parseForceSweep(tensdir, ...
            theseRunVals_t,cytoparams);
        tensInitCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM_t;
        tensFinalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM_t;
        tensDisps.(networkLabel)(:,:,repIdx) = thisFinalCoM_t - thisInitCoM_t;
        % estimate a modulus
        [thisTensModulus, linefit_t, gof_t] = getModulus(thisInitCoM_t, ...
            thisFinalCoM_t,tensForces,"extension",filtering_t);
        tensModuli(netIdx,repIdx) = thisTensModulus;
        tensLineFits.(networkLabel){repIdx} = linefit_t;
        tensFitStats.(networkLabel){repIdx} = gof_t;
        
    end
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
        'nFilPerAsterList','dens','nNetTypes','numRep',...
        'filtering_s','filtering_t',...
        'shearForces','shearModuli','shearLineFits','shearFitStats',...
        'shearInitCoMs','shearFinalCoMs','shearDisps','percStatus_s',...
        'tensForces','tensModuli','tensLineFits','tensFitStats', ...
        'tensInitCoMs','tensFinalCoMs','tensDisps','percStatus_t')
    clearvars -except nameprefix matFile saveDir dataSubdir figSubdir
end

load(fullfile(saveDir,dataSubdir,matFile),'D','len_fil', ...
    'nFilPerAsterList','dens','nNetTypes','numRep', ...
    'filtering_s','filtering_t',...
    'shearForces','shearModuli','shearLineFits','shearFitStats',...
    'shearInitCoMs','shearFinalCoMs','shearDisps','percStatus_s',...
    'tensForces','tensModuli','tensLineFits','tensFitStats', ...
    'tensInitCoMs','tensFinalCoMs','tensDisps','percStatus_t')

%% Histograms of shear/tensile moduli

fig1 = figure(1);
plotTiling = [4,2];
set(fig1,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
    'defaultTextInterpreter','tex');
fig1 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
tileCount = prod(plotTiling);

meanShearMod = mean(shearModuli,2);
meanTensMod = mean(tensModuli,2);

tileIdx = 0;
% shear moduli on left | tensile moduli on right
for netIdx = 1:nNetTypes
    tileIdx = tileIdx + 1;
    nexttile(tileIdx)
    % plot shear modulus histogram
    h = histogram(shearModuli(netIdx,:),'Normalization','probability');
    % h.BinWidth = 5;
    % xlim([0,70])
    % ylim([0,0.5])
    xlabel('Shear modulus [pN/{\mu}m]')
    ylabel('Empirical probability')
    title(sprintf('a_n = %i (mean = %2.2f)',nFilPerAsterList(netIdx), ...
        meanShearMod(netIdx)))
    tileIdx = tileIdx + 1;

    nexttile(tileIdx)
    % plot tensile modulus histogram
    h = histogram(tensModuli(netIdx,:),'Normalization','probability');
    % h.BinWidth = 50;
    % xlim([0,500])
    % ylim([0,0.5])
    xlabel('Tensile modulus [pN/{\mu}m]')
    ylabel('Empirical probability')
    title(sprintf('a_n = %i (mean = %2.2f)',nFilPerAsterList(netIdx), ...
        meanTensMod(netIdx)))

    if tileIdx == tileCount && netIdx < nNetTypes % end of first page
        exportgraphics(fig1,fullfile(saveDir,figSubdir, ...
            [nameprefix,'-hist.pdf']),'Resolution',300)
        fig1 = clf(fig1);
        set(fig1,'units','inches','Position',[0,0,8.5,11],'visible', ...
            'off','defaultTextInterpreter','tex');
        fig1 = tiledlayout(plotTiling(1),plotTiling(2), ...
            'TileSpacing','compact');
        tileIdx = 0;
    elseif netIdx == nNetTypes % end of second (last) page
        exportgraphics(fig1,fullfile(saveDir,figSubdir,[nameprefix,'-hist.pdf']), ...
            'Resolution',300,'Append',true)
        fig1 = clf(fig1);
        set(fig1,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultTextInterpreter','tex');
        fig1 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    end
end


