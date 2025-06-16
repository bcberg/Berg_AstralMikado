% rep_shear_aster_otherRhos.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "rep_shear_otherRhos250324";
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';

newReports = false;
if newReports
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/rep_full_shear_aster_otherRhos/";
    
    % ensure the following parameters match those used in
    % rep_full_shear_aster.cym.tpl and its driver (for numRep)
    dens = [6,9];  % N * len_fil / A
    numDens = length(dens);
    len_fil = 1;
    D = 10; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = (1:24)';
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    forceVals = 5/3 * (0:3)';
    nForces = length(forceVals);
    numRep = 15;
    
    runsPerDens = nForces * nNetTypes;
    runsPerRep = nForces * nNetTypes * numDens;
    shearModuli = zeros(nNetTypes,numDens,numRep);
    for repIdx = 1:numRep
        for densIdx = 1:numDens
            for netIdx = 1:nNetTypes
                nFilPerAster = nFilPerAsterList(netIdx);
                netLabel = sprintf('an%02i',nFilPerAster);
                if repIdx == 1
                    initCoMs.(netLabel) = zeros(nForces,2,numDens,numRep);
                    finalCoMs.(netLabel) = zeros(nForces,2,numDens,numRep);
                    disps.(netLabel) = zeros(nForces,2,numDens,numRep);
                    lineFits.(netLabel) = cell(numDens,numRep);
                    fitStats.(netLabel) = cell(numDens,numRep);
                end
                cytoparams = struct('nFil',nFil(densIdx),'nFilPerAster', ...
                    nFilPerAster,'dataStartLine',dataStartLine, ...
                    'ptsPerFiber',ptsPerFiber);
                unshiftedRunVals = (netIdx - 1)*nForces:(netIdx * nForces - 1);
                densShift = (densIdx - 1) * runsPerDens;
                repShift = (repIdx - 1) * runsPerRep;
                theseRunVals = densShift + repShift + unshiftedRunVals;
                [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,theseRunVals,...
                    cytoparams);
                initCoMs.(netLabel)(:,:,densIdx,repIdx) = thisInitCoM;
                finalCoMs.(netLabel)(:,:,densIdx,repIdx) = thisFinalCoM;
                disps.(netLabel)(:,:,densIdx,repIdx) = thisFinalCoM - ...
                    thisInitCoM;
                [thisModulus, linefit, gof] = getModulus(thisInitCoM,thisFinalCoM,...
                    forceVals,"shear");
                shearModuli(netIdx,densIdx,repIdx) = thisModulus;
                lineFits.(netLabel){densIdx,repIdx} = linefit;
                fitStats.(netLabel){densIdx,repIdx} = gof;
            end
        end
    end
    save(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
        'D','nFil','nFilPerAsterList','forceVals','numRep','nNetTypes', ...
        'numDens','initCoMs','finalCoMs','disps','lineFits','fitStats', ...
        'shearModuli')
    clearvars -except simName saveDir
end

load(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
    'D','nFil','nFilPerAsterList','forceVals','numRep','nNetTypes', ...
    'numDens','initCoMs','finalCoMs','disps','lineFits','fitStats', ...
    'shearModuli')

lgray = [197,190,181]/255;

%% Modulus vs astral number at each density

[stdevMod,meanMod] = std(shearModuli,0,3); % compute summary stats across reps
z = norminv(0.975);
markerList = {'o','^'};
fig1 = figure(1);
hold on
for densIdx = 1:numDens
    for netIdx = 1:nNetTypes
        plot(nFilPerAsterList(netIdx),squeeze(shearModuli(netIdx,densIdx,:)), ...
            'Marker',markerList{densIdx},'MarkerSize',3,'MarkerFaceColor',lgray, ...
            'MarkerEdgeColor', lgray)
    end
end
l1 = errorbar(nFilPerAsterList,meanMod(:,1), ...
    z*stdevMod(:,1)/sqrt(numRep),'-ob');
l2 = errorbar(nFilPerAsterList,meanMod(:,2), ...
        z*stdevMod(:,2)/sqrt(numRep),'-^r');
hold off
xlabel('Astral number a_n','Interpreter','tex')
ylabel('2D shear modulus [pN/{\mu}m]','Interpreter','tex')
legend([l1,l2],{sprintf('\\rho=%1.1f',dens(1)), ...
    sprintf('\\rho=%1.1f',dens(2))},'Location','eastoutside', ...
    'Interpreter','tex')
exportgraphics(fig1,fullfile(saveDir,'exploratory_figures', ...
    simName + ".pdf"),'ContentType','vector')

% % summarizing R^2 stats
% Rsquares = zeros(nNetTypes,numRep);
% for netIdx = 1:numel(Rsquares)
%     [netIdx,repIdx] = ind2sub([nNetTypes,numRep],netIdx);
%     networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
%     thisGOF = fitStats.(networkLabel){repIdx};
%     Rsquares(netIdx,repIdx) = thisGOF.rsquare;
% end
% fig2 = figure(2);
% hold on
% for repIdx = 1:numRep
%     plot(nFilPerAsterList,Rsquares(:,repIdx),'DisplayName',...
%         sprintf('t%02i',repIdx),'LineWidth',1)
% end
% hold off
% xlabel('Number of filaments per aster')
% ylabel('R^2 value for linear fit')
% legend('Location','eastoutside')

%% Displacement vs. force plots

% fig3 = figure(3);
% set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
%     'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
% plotTiling = [6,5];
% fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
% 
% for netIdx = 1:nNetTypes
%     tileIdx = 1;
%     networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
%     for repIdx = 1:numRep
%         nexttile(fig3,tileIdx)
%         hold on
%         plot(forceVals, disps.(networkLabel)(:,1,repIdx), '*-b')
%         plot(forceVals, forceVals / shearModuli(netIdx,repIdx), '-k')
%         hold off
%         xlim('tight')
%         yl = ylim;
%         ylim([0,yl(2)])
%         xlabel('Shear force [pN]')
%         ylabel('x-Displ. [{\mu}m]')
%         tileIdx = tileIdx + 1;
%     end
%     title(fig3,networkLabel)
%     if netIdx == 1
%         % export first page
%         lg = legend("data","OLS fit");
%         lg.NumColumns = 2;
%         lg.Layout.Tile = 'south';
%         exportgraphics(fig3,fullfile(saveDir,[nameprefix,'-linefits.pdf']), ...
%             'Resolution',300)
%         fig3 = clf(fig3);
%         set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
%             'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
%         fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
%     else
%         % export next (or final) page(s), append to first file
%         lg = legend("data","OLS fit");
%         lg.NumColumns = 2;
%         lg.Layout.Tile = 'south';
%         exportgraphics(fig3,fullfile(saveDir,[nameprefix,'-linefits.pdf']), ...
%             'Resolution',300,'Append',true)
%         fig3 = clf(fig3);
%         set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
%             'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
%         fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
%     end
% end
