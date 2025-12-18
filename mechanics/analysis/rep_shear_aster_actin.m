% rep_shear_aster_actin.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

nameprefix = 'rep_shear_actin251015';
newReports = false;
if newReports
    filtering.TF = true;
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/" + ...
        "rep_full_shear_aster_actin/";
    
    % ensure the following parameters match those used in
    % rep_full_shear_aster_actin.cym.tpl and its driver (for numRep)
    dens = 75;  % N * len_fil / A
    len_fil = 0.1;
    D = 1; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = (1:24)';
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    forceVals = 0.05 * (0:5)';
    nForces = length(forceVals);
    numRep = 30;
    
    runsPerRep = nNetTypes * nForces;
    shearModuli = zeros(nNetTypes,numRep);
    percStatus = false(nNetTypes,numRep,2);
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
        % check network percolation status
        [percTF,~,~] = checkLinks(dir,theseRunVals(1),cytoparams,4);
        filtering.spanCheck = percTF(1);
        percStatus(netIdx,repIdx,:) = reshape(percTF,[1,1,2]);
        % analyze network position data
        [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,theseRunVals,...
            cytoparams);
        initCoMs.(networkLabel)(:,:,repIdx) = thisInitCoM;
        finalCoMs.(networkLabel)(:,:,repIdx) = thisFinalCoM;
        disps.(networkLabel)(:,:,repIdx) = thisFinalCoM - thisInitCoM;
        % estimate a modulus
        [thisModulus, linefit, gof] = getModulus(thisInitCoM,thisFinalCoM,...
            forceVals,"shear",filtering);
        shearModuli(netIdx,repIdx) = thisModulus;
        lineFits.(networkLabel){repIdx} = linefit;
        fitStats.(networkLabel){repIdx} = gof;
    end
    filename = sprintf([nameprefix,'_l%1.1f_D%i.mat'],len_fil,D);
    % CHECK FILENAME to avoid overwriting analysis of past runs
    save(['~/Documents/astral-mikado-data/mat_files/',filename], ...
        'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
        'nNetTypes','initCoMs','finalCoMs','forceVals','disps','dens', ...
        'percStatus','filtering')
    clearvars -except filename nameprefix
else
    len_fil = 0.1;
    D = 1;
    filename = sprintf([nameprefix,'_l%1.1f_D%i.mat'],len_fil,D);
end
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data/';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data\';

load(fullfile(saveDir,'mat_files',filename), ...
    'nFilPerAsterList','shearModuli','lineFits','fitStats','numRep',...
    'nNetTypes','initCoMs','finalCoMs','forceVals','disps','dens', ...
    'percStatus','filtering')


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
ylim([0, inf])
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

fig3 = figure(3); clf;
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
        if shearModuli(netIdx,repIdx) == 0
            l1 = plot(forceVals, disps.(networkLabel)(:,1,repIdx), '*-r', ...
                'DisplayName','data');
            if ~percStatus(netIdx,repIdx,1) % if no spanning comp.
                text(forceVals(3),0.1*max(disps.(networkLabel)(:,1,repIdx)), ...
                    'no span')
            end
        else
            l1 = plot(forceVals, disps.(networkLabel)(:,1,repIdx), '*-b', ...
                'DisplayName','data');
            l2 = plot(forceVals, forceVals / shearModuli(netIdx,repIdx), ...
                '-k','DisplayName','linefit');
        end
        hold off
        xlim('tight')
        % yl = ylim;
        % ylim([0,yl(2)])

        tileIdx = tileIdx + 1;
    end
    xlabel(fig3,'Shear force [pN]')
    ylabel(fig3,'x-Displ. [{\mu}m]')
    title(fig3,networkLabel)
    if netIdx == 1
        % export first page
        if ishandle(l2)
            lg = legend([l1,l2],'NumColumns',2);
        else
            lg = legend(l1,'NumColumns',1);
        end
        lg.Layout.Tile = 'south';
        exportgraphics(fig3,fullfile(saveDir,'exploratory_figures', ...
            [nameprefix,'-linefits.pdf']),'Resolution',300)
        fig3 = clf(fig3);
        set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
        fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    else
        % export next (or final) page(s), append to first file
        if ishandle(l2)
            lg = legend([l1,l2],'NumColumns',2);
        else
            lg = legend(l1,'NumColumns',1);
        end
        lg.Layout.Tile = 'south';
        exportgraphics(fig3,fullfile(saveDir,'exploratory_figures', ...
            [nameprefix,'-linefits.pdf']),'Resolution',300,'Append',true)
        fig3 = clf(fig3);
        set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
            'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
        fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
    end
end
