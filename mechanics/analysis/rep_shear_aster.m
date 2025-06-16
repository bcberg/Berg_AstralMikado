% rep_shear_aster.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

nameprefix = 'rep_shear241105';
newReports = false;
if newReports
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/rep_full_shear_aster/";
    
    % ensure the following parameters match those used in
    % rep_full_shear_aster.cym.tpl and its driver (for numRep)
    dens = 7.5;  % N * len_fil / A
    len_fil = 1;
    D = 10; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = (1:24)';
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    forceVals = (0:5)';
    nForces = length(forceVals);
    numRep = 30;
    
    runsPerRep = nNetTypes * nForces;
    shearModuli = zeros(nNetTypes,numRep);
    for idx = 1:numel(shearModuli)
        [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
        nFilPerAster = nFilPerAsterList(netIdx);
        networkLabel = sprintf('an%02i',nFilPerAster);

        if repIdx == 1
            initCoMs.(networkLabel) = zeros(nForces,2,numRep);
            finalCoMs.(networkLabel) = zeros(nForces,2,numRep);
            disps.(networkLabel) = zeros(nForces,2,numRep);
            lineFits.(networkLabel) = cell(numRep,1);
            fitStats.(networkLabel) = cell(numRep,1);
        end

        cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
            'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);
        unshiftedRunVals = (netIdx - 1)*nForces:(netIdx * nForces - 1);
        theseRunVals = (repIdx - 1) * runsPerRep + unshiftedRunVals;
        [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,theseRunVals,...
            cytoparams);
        initCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM;
        finalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM;
        disps.(networkLabel)(:,:,repIdx) = thisFinalCoM - thisInitCoM;
        [thisModulus, linefit, gof] = getModulus(thisInitCoM,thisFinalCoM,...
            forceVals,"shear");
        shearModuli(netIdx,repIdx) = thisModulus;
        lineFits.(networkLabel){repIdx} = linefit;
        fitStats.(networkLabel){repIdx} = gof;
    end
    filename = sprintf([nameprefix,'_l%02i_D%i.mat'],len_fil,D);
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(['~/Documents/astral-mikado-data/mat_files/',filename], ...
        'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
        'nNetTypes','initCoMs','finalCoMs','forceVals','disps','dens')
    clearvars -except filename nameprefix
else
    len_fil = 1;
    D = 10;
    filename = sprintf([nameprefix,'_l%02i_D%i.mat'],len_fil,D);
end
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data/';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\AstralMikadoCYM\data';

load(fullfile(saveDir,'mat_files',filename), ...
    'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
    'nNetTypes','initCoMs','finalCoMs','forceVals','disps','dens')


%% Summary plots

[stdevs,means] = std(shearModuli,0,2); % compute summary stats across cols
z = norminv(0.975);
fig1 = figure(1);
hold on
errorbar(nFilPerAsterList,means,z*stdevs/sqrt(numRep),'-ob')
for idx = 1:nNetTypes
    plot(nFilPerAsterList(idx),shearModuli(idx,:),'*k','MarkerSize',4)
end
hold off
xlabel('Astral number')
ylabel('2D shear modulus [pN/{\mu}m]')
exportgraphics(fig1,fullfile(saveDir,'exploratory_figures', ...
    [nameprefix,'-moduli.png']),'Resolution',300)

% summarizing R^2 stats
Rsquares = zeros(nNetTypes,numRep);
for idx = 1:numel(Rsquares)
    [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
    networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
    thisGOF = fitStats.(networkLabel){repIdx};
    Rsquares(netIdx,repIdx) = thisGOF.rsquare;
end
fig2 = figure(2);
hold on
for repIdx = 1:numRep
    plot(nFilPerAsterList,Rsquares(:,repIdx),'DisplayName',...
        sprintf('t%02i',repIdx),'LineWidth',1)
end
hold off
xlabel('Number of filaments per aster')
ylabel('R^2 value for linear fit')
legend('Location','eastoutside')

%% Displacement vs. force plots

fig3 = figure(3);
set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
    'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
plotTiling = [6,5];
fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');

for netIdx = 1:nNetTypes
    tileIdx = 1;
    networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
    for repIdx = 1:numRep
        nexttile(fig3,tileIdx)
        hold on
        plot(forceVals, disps.(networkLabel)(:,1,repIdx), '*-b')
        plot(forceVals, forceVals / shearModuli(netIdx,repIdx), '-k')
        hold off
        xlim('tight')
        yl = ylim;
        ylim([0,yl(2)])
        xlabel('Shear force [pN]')
        ylabel('x-Displ. [{\mu}m]')
        tileIdx = tileIdx + 1;
    end
    title(fig3,networkLabel)
    if netIdx == 1
        % export first page
        lg = legend("data","OLS fit");
        lg.NumColumns = 2;
        lg.Layout.Tile = 'south';
        exportgraphics(fig3,fullfile(saveDir,'exploratory_figures', ...
            [nameprefix,'-linefits.pdf']),'Resolution',300)
        fig3 = clf(fig3);
        set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
        fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    else
        % export next (or final) page(s), append to first file
        lg = legend("data","OLS fit");
        lg.NumColumns = 2;
        lg.Layout.Tile = 'south';
        exportgraphics(fig3,fullfile(saveDir,'exploratory_figures', ...
            [nameprefix,'-linefits.pdf']),'Resolution',300,'Append',true)
        fig3 = clf(fig3);
        set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
        fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    end
end
