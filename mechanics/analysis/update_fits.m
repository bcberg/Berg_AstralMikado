% update_fits.m
% Brady Berg, 10/13/2025
clear
close all
format short
format compact

%% Specify data to refit

nameprefix = 'rep_shear_actin251015';
len_fil = 0.1;
D = 1; % side length of square region

% Ubuntu path
saveDir = '~/Documents/astral-mikado-data/';
% Windows path
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data\';

dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';
refit = false;
if refit
    filtering.TF = true;
    matfileName = [nameprefix,sprintf('_l%1.1f_D%i.mat',len_fil,D)];
    load(fullfile(saveDir,dataSubdir,matfileName), ...
        'nFilPerAsterList','numRep','nNetTypes','initCoMs','finalCoMs', ...
        'forceVals','disps','dens','percStatus')

    moduli = zeros(nNetTypes,numRep);
    % intercepts = zeros(nNetTypes,numRep);
    for idx = 1:numel(moduli)
        [netIdx,repIdx] = ind2sub([nNetTypes,numRep],idx);
        nFilPerAster = nFilPerAsterList(netIdx);
        netLabel = sprintf('an%02d',nFilPerAster);

        if repIdx == 1
            lineFits.(netLabel) = cell(numRep,1);
            fitStats.(netLabel) = cell(numRep,1);
        end
        filtering.spanCheck = percStatus(netIdx,repIdx,1);
        if contains(nameprefix,'shear')
            [thisModulus,linefit,gof] = getModulus(initCoMs.(netLabel)(:,:,repIdx), ...
                finalCoMs.(netLabel)(:,:,repIdx),forceVals,"shear", ...
                filtering);
        elseif contains(nameprefix,'tens')
            [thisModulus,linefit,gof] = getModulus(initCoMs.(netLabel)(:,:,repIdx), ...
                finalCoMs.(netLabel)(:,:,repIdx),forceVals,"extension", ...
                filtering);
        end
        moduli(netIdx,repIdx) = thisModulus;
        % intercepts(netIdx,repIdx) = linefit.p2;
        lineFits.(netLabel){repIdx} = linefit;
        fitStats.(netLabel){repIdx} = gof;
    end

    % giving specific names to the moduli and saving
    if contains(nameprefix,'shear')
        shearModuli = moduli;
        save(fullfile(saveDir,dataSubdir,matfileName), ...
            'nFilPerAsterList','numRep','nNetTypes','initCoMs','finalCoMs', ...
            'forceVals','disps','dens','shearModuli', ...
            'lineFits','fitStats','percStatus','filtering')
    elseif contains(nameprefix,'tens')
        tensModuli = moduli;
        save(fullfile(saveDir,dataSubdir,matfileName), ...
            'nFilPerAsterList','numRep','nNetTypes','initCoMs','finalCoMs', ...
            'forceVals','disps','dens','tensModuli','intercepts', ...
            'lineFits','fitStats','percStatus','filtering')
    end
    clear('moduli','idx','netIdx','repIdx','thisModulus','linefit','gof', ...
        'netLabel','nFilPerAster')
end

%% Remake plots: modulus vs. astral number, R^2 overview

if refit && contains(nameprefix,'shear')
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
    exportgraphics(fig1,fullfile(saveDir,figSubdir, ...
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
    legend('Location','eastoutside')

elseif refit && contains(nameprefix,'tens')
    [stdevs,means] = std(tensModuli,0,2); % compute summary stats across cols
    z = norminv(0.975);
    fig1 = figure(1); clf;
    hold on
    for idx = 1:nNetTypes
        plot(nFilPerAsterList(idx),tensModuli(idx,:),'*k','MarkerSize',4)
    end
    errorbar(nFilPerAsterList,means,z*stdevs/sqrt(numRep),'-ob','LineWidth',1.5)
    hold off
    xlabel('Astral number')
    ylabel('2D tensile modulus [pN/{\mu}m]')
    % ylim([0, inf])
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
    legend('Location','eastoutside')

end

%% Remake plots: displacement vs. force

if refit
    fig3 = figure(3); clf;
    set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
        'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
    plotTiling = [6,5];
    fig3 = tiledlayout(fig3,plotTiling(1),plotTiling(2),'TileSpacing','compact');
end

if refit && contains(nameprefix,'shear')
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
            exportgraphics(fig3,fullfile(saveDir,figSubdir, ...
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
            exportgraphics(fig3,fullfile(saveDir,figSubdir, ...
                [nameprefix,'-linefits.pdf']),'Resolution',300,'Append',true)
            fig3 = clf(fig3);
            set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
                'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
            fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
        end
    end
elseif refit && contains(nameprefix,'tens')
    for netIdx = 1:nNetTypes
        tileIdx = 1;
        networkLabel = sprintf('an%02i',nFilPerAsterList(netIdx));
        for repIdx = 1:numRep
            nexttile(fig3,tileIdx)
            hold on
            if tensModuli(netIdx,repIdx) == 0
                l1 = plot(forceVals, disps.(networkLabel)(:,2,repIdx), '*-r', ...
                    'DisplayName','data');
                if ~percStatus(netIdx,repIdx,1) % if no spanning comp.
                    text(forceVals(3),0.1*max(disps.(networkLabel)(:,2,repIdx)), ...
                        'no span')
                end
            else
                l1 = plot(forceVals, disps.(networkLabel)(:,2,repIdx), '*-b', ...
                    'DisplayName','data');
                l2 = plot(forceVals, forceVals / tensModuli(netIdx,repIdx), ...
                    '-k','DisplayName','linefit');
            end
            hold off
            xlim('tight')
            yl = ylim;
            % ylim([0,yl(2)])
            tileIdx = tileIdx + 1;
        end
        xlabel(fig3,'Tensile force [pN]')
        ylabel(fig3,'y-Displ. [{\mu}m]')
        title(fig3,networkLabel)
        if netIdx == 1
            % export first page
            if ishandle(l2)
                lg = legend([l1,l2],'NumColumns',2);
            else
                lg = legend(l1,'NumColumns',1);
            end
            lg.Layout.Tile = 'south';
            exportgraphics(fig3,fullfile(saveDir,figSubdir, ...
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
            exportgraphics(fig3,fullfile(saveDir,figSubdir, ...
                [nameprefix,'-linefits.pdf']),'Resolution',300,'Append',true)
            fig3 = clf(fig3);
            set(fig3,'units','inches','Position',[0,0,8.5,11],'visible','off', ...
                'defaultLineLineWidth',0.5,'defaultTextInterpreter','tex')
            fig3 = tiledlayout(plotTiling(1),plotTiling(2),'TileSpacing','compact');
        end
    end
end