% modulus_distr.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

nameprefix = 'mod_distr241120';
newReports = false;
if newReports
    rundir = "/home/bcberg/Documents/AstralMikadoCYM/runs/";
    sheardir = fullfile(rundir,"modulus_distr_shear/");
    tensdir = fullfile(rundir,"modulus_distr_tens/");
    
    % ensure the following parameters match those used in
    % modulus_distr_****.cym.tpl and its driver (for numRep)
    dens = 7.5;  % N * len_fil / A
    len_fil = 1;
    D = 10; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = transpose([1:3,4:4:24]);
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    shearForces = transpose((0:3) * 5/3);
    tensForces = transpose((0:3) * 20/3);
    nShearForces = length(shearForces);
    nTensForces = length(tensForces);
    numRep = 100;
    
    runsPerRep = nNetTypes * nShearForces;
    shearModuli = zeros(nNetTypes,numRep);
    tensModuli = zeros(nNetTypes,numRep);
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
        [thisInitCoM_s,thisFinalCoM_s] = parseForceSweep(sheardir, ...
            theseRunVals_s,cytoparams);
        shearInitCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM_s;
        shearFinalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM_s;
        shearDisps.(networkLabel)(:,:,repIdx) = thisFinalCoM_s - thisInitCoM_s;
        [thisShearModulus, linefit_s, gof_s] = getModulus(thisInitCoM_s, ...
            thisFinalCoM_s,shearForces,"shear");
        shearModuli(netIdx,repIdx) = thisShearModulus;
        shearLineFits.(networkLabel){repIdx} = linefit_s;
        shearFitStats.(networkLabel){repIdx} = gof_s;

        unshiftedRunVals_t = (netIdx - 1)*nShearForces:(netIdx * nShearForces - 1);
        theseRunVals_t = (repIdx - 1) * runsPerRep + unshiftedRunVals_t;
        [thisInitCoM_t,thisFinalCoM_t] = parseForceSweep(tensdir, ...
            theseRunVals_t,cytoparams);
        tensInitCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM_t;
        tensFinalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM_t;
        tensDisps.(networkLabel)(:,:,repIdx) = thisFinalCoM_t - thisInitCoM_t;
        [thisTensModulus, linefit_t, gof_t] = getModulus(thisInitCoM_t, ...
            thisFinalCoM_t,tensForces,"extension");
        tensModuli(netIdx,repIdx) = thisTensModulus;
        tensLineFits.(networkLabel){repIdx} = linefit_t;
        tensFitStats.(networkLabel){repIdx} = gof_t;
    end
    filename = sprintf([nameprefix,'_l%02i_D%i.mat'],len_fil,D);
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(['~/Documents/AstralMikadoCYM/data/',filename], ...
        'nFilPerAsterList','dens','nNetTypes','numRep',...
        'shearForces','shearModuli','shearLineFits','shearFitStats',...
        'shearInitCoMs','shearFinalCoMs','shearDisps',...
        'tensForces','tensModuli','tensLineFits','tensFitStats', ...
        'tensInitCoMs','tensFinalCoMs','tensDisps')
    clearvars -except filename nameprefix
else
    len_fil = 1;
    D = 10;
    filename = sprintf([nameprefix,'_l%02i_D%i.mat'],len_fil,D);
end
% Ubuntu filepath
% saveDir = '~/Documents/AstralMikadoCYM/data';
% Windows filepath
saveDir = 'C:\Users\bcber\Documents\AstralMikadoCYM\data';

load(fullfile(saveDir,filename), ...
    'nFilPerAsterList','dens','nNetTypes','numRep',...
    'shearForces','shearModuli','shearLineFits','shearFitStats',...
    'shearInitCoMs','shearFinalCoMs','shearDisps',...
    'tensForces','tensModuli','tensLineFits','tensFitStats', ...
    'tensInitCoMs','tensFinalCoMs','tensDisps')

%% Histograms of shear/tensile moduli

fig1 = figure(1);
plotTiling = [3,2];
set(fig1,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
    'defaultTextInterpreter','tex');
fig1 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
tileCount = prod(plotTiling);

meanShearMod = mean(shearModuli,2);
meanTensMod = mean(tensModuli,2);

netIdx = 1;
for pageIdx = 1:(ceil(nNetTypes/plotTiling(1)))
    % shear moduli on left | tensile moduli on right
    for tileIdx = 1:tileCount
        nexttile(tileIdx)
        if mod(tileIdx,2) == 1
            % plot shear modulus histogram
            h = histogram(shearModuli(netIdx,:),'Normalization','probability');
            h.BinWidth = 5;
            xlim([0,70])
            % ylim([0,0.5])
            xlabel('Shear modulus [pN/{\mu}m]')
            ylabel('Empirical probability')
            title(sprintf('a_n = %i (mean = %2.2f)',nFilPerAsterList(netIdx), ...
                meanShearMod(netIdx)))
        elseif mod(tileIdx,2) == 0
            % plot tensile modulus histogram
            h = histogram(tensModuli(netIdx,:),'Normalization','probability');
            h.BinWidth = 50;
            xlim([0,500])
            % ylim([0,0.5])
            xlabel('Tensile modulus [pN/{\mu}m]')
            ylabel('Empirical probability')
            title(sprintf('a_n = %i (mean = %2.2f)',nFilPerAsterList(netIdx), ...
                meanTensMod(netIdx)))

            % move on to next network type
            netIdx = netIdx + 1;
        end
    end
    if tileIdx == tileCount && pageIdx == 1
        exportgraphics(fig1,fullfile(saveDir,[nameprefix,'-hist.pdf']), ...
            'Resolution',300)
        fig1 = clf(fig1);
        set(fig1,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultTextInterpreter','tex');
        fig1 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    elseif tileIdx == tileCount && pageIdx > 1
        exportgraphics(fig1,fullfile(saveDir,[nameprefix,'-hist.pdf']), ...
            'Resolution',300,'Append',true)
        fig1 = clf(fig1);
        set(fig1,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultTextInterpreter','tex');
        fig1 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    end
end


