% percCurveFitting.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Importing data

% Ubuntu path
saveDir = "~/Documents/astral-mikado-data/";
% Windows path
% saveDir = "C:\Users\bcber\Documents\astral-mikado-data";
matFileSubDir = "mat_files";
figSubDir = "exploratory_figures";
filePattern = "percProbs_l%02i_D%02i";

l = 1;
Dlist = [5,10,15,20];
numD = numel(Dlist);
allCurves = cell(1,numD);
allAstralNums = [];
for idx = 1:numD
    load(fullfile(saveDir,matFileSubDir, ...
        sprintf(filePattern,l,Dlist(idx))),'curves','astralNumList')
    allCurves{idx} = curves;
    allAstralNums = union(allAstralNums,astralNumList);
    clear curves
end
numNetTypes = length(allAstralNums);
percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
numPercTypes = length(percTypes); % i.e., = 5

%% Overplotting different system sizes (family of curves)

makeRawDataFig = false;
if makeRawDataFig
    rawDataFig = figure('units','inches','Position',[1,1,8.5,11]); clf;
    set(rawDataFig,'visible','off','defaultLineLineWidth',0.6)
    for idx = 1:numNetTypes
        tl = tiledlayout(3,2,'TileSpacing','compact');
        netLabel = sprintf('an%02i',allAstralNums(idx));
        tileIdx = 1;
        for jdx = 1:numPercTypes
            nexttile(tl,tileIdx)
            hold on
            for kdx = 1:numD
                curves = allCurves{kdx};
                if isfield(curves,netLabel)
                    plot(curves.(netLabel)(1,:),curves.(netLabel)(1+jdx,:), ...
                        'DisplayName',sprintf('D=%2.1f',Dlist(kdx)))
                end
            end
            hold off
            xscale('log')
            xlabel('Filament density [$\mu m^{-1}$]')
            ylabel('Percolation probability')
            title(percTypes{jdx})
            if jdx == 1
                lg = legend('FontSize',12);
                lg.Layout.Tile = 6;
            end
            tileIdx = tileIdx + 1;
        end
        title(tl,netLabel)
        if idx == 1
            % export first page
            exportgraphics(rawDataFig,fullfile(saveDir,'percMCRG_data.pdf'), ...
                'Resolution',300)
            rawDataFig = clf(rawDataFig);
        else
            % export next pages
            exportgraphics(rawDataFig,fullfile(saveDir,'percMCRG_data.pdf'), ...
                'Resolution',300,'Append',true)
            rawDataFig = clf(rawDataFig);
        end
        set(rawDataFig,'units','inches','Position',[1,1,8.5,11],'visible','off', ...
            'defaultLineLineWidth',0.6)
    end
end

%% Fit sigmoidal models

% default logistic model has form f(x) = a / (1+e^{-b*(x-c)})
% Currently fixing a=1
% Gompertz model has form g(x) = d + (a-d)e^{-e^{-b*(x-c)}}
% Currently fixing a=1, d=0

makeSigFitsFig = false;
if makeSigFitsFig
    sigFits = cell(numNetTypes,numPercTypes,numD);
    fitopts_logis = fitoptions('Method','NonlinearLeastSquares','Lower', ...
        [1,-5,-5],'Upper',[1,20,100],'StartPoint',[1,1,10]);
    fitopts_gomp = fitoptions('Method','NonlinearLeastSquares','Lower', ...
        [1,0,-5,0],'Upper',[1,25,100,0],'StartPoint',[1,1,10,0]);
    sigFitsFig = figure('units','inches','Position',[1,1,8.5,11]); clf;
    set(sigFitsFig,'defaultLineLineWidth',0.4,'defaultLineMarkerSize',3,...
        'visible','off')
    C = colororder("gem");
    for idx = 1:numNetTypes
        legendLabels = cell(2*numD,1);
        t = tiledlayout(3,2,'TileSpacing','compact');
        netLabel = sprintf('an%02i',allAstralNums(idx));
        tileIdx = 1;
        for jdx = 1:numPercTypes
            nexttile(t,tileIdx)
            hold on
            for kdx = 1:numD
                curves = allCurves{kdx};
                if isfield(curves,netLabel)
                    x = transpose(curves.(netLabel)(1,:));
                    y = transpose(curves.(netLabel)(1+jdx,:));
                    % restricting fit windows to exclude flat lines in data
                    % at p=0 and p=1
                    greaterThanZero = find(y>0,1);
                    % discard decreasing portions of connectivity perc
                    % curves that appear at low density & high astral num.
                    if greaterThanZero==1
                        lowerLim = 1;
                        while y(lowerLim+1) < y(lowerLim)
                            lowerLim = lowerLim + 1;
                        end
                    else
                        lowerLim = greaterThanZero - 1;
                    end
                    % upper limit is where the curve increases to 1 for the
                    % first time, otherwise the last data point
                    if any(y==1)
                        upperLimOptions = find(y==1);
                        reachedOne = find(upperLimOptions>lowerLim,1);
                        upperLim = upperLimOptions(reachedOne);
                    else
                        upperLim = length(y);
                    end
                    x = x(lowerLim:upperLim);
                    y = y(lowerLim:upperLim);
                    if jdx <= 4
                        sigFits{idx,jdx,kdx} = fit(x,y,'logistic',fitopts_logis);
                    elseif jdx == 5
                        sigFits{idx,jdx,kdx} = fit(x,y,'gompertz',fitopts_gomp);
                    end
                    l = plot(sigFits{idx,jdx,kdx},x,y,'*');
                    [l.Color] = deal(C(kdx,:));
                    legendLabels{2*kdx - 1} = sprintf('D=%2.1f data', ...
                        Dlist(kdx));
                    legendLabels{2*kdx} = sprintf('D=%2.1f fit', ...
                        Dlist(kdx));
                end
            end
            hold off
            xscale('log')
            xlabel('Filament density [$\mu m^{-1}$]')
            ylabel('Percolation probability')
            title(percTypes{jdx})
            if jdx == 1
                legendLabels = legendLabels(~cellfun('isempty', ...
                    legendLabels));
                lg = legend(legendLabels,'FontSize',12);
                lg.Layout.Tile = 6;
            else
                legend('off')
            end
            tileIdx = tileIdx + 1;
        end
        title(t,netLabel)
        if idx == 1
            % export first page
            exportgraphics(sigFitsFig,fullfile(saveDir,figSubDir, ...
                'percMCRG_paramFits_tzOnly.pdf'),'Resolution',300)
            sigFitsFig = clf(sigFitsFig);
        else
            % export next pages
            exportgraphics(sigFitsFig,fullfile(saveDir,figSubDir, ...
                'percMCRG_paramFits_tzOnly.pdf'),'Resolution',300,'Append',true)
            sigFitsFig = clf(sigFitsFig);
        end
        set(sigFitsFig,'units','inches','Position',[1,1,8.5,11], ...
            'defaultLineLineWidth',0.4,'defaultLineMarkerSize',3, ...
            'visible','off')
    end
    save(fullfile(saveDir,matFileSubDir,'percMCRG_paramFits_tzOnly.mat'), ...
        'sigFits','Dlist')
end

%% Nonparametric fitting

% current scheme: smoothing spline with automatic smoothing parameter
makeSmsplFitsFig = true;
if makeSmsplFitsFig
    smsplFits = cell(numNetTypes,numPercTypes,numD);
    percDataUsed = cell(numNetTypes,numPercTypes,numD);
    smsplFitsFig = figure('units','inches','Position',[1,1,8.5,11]); clf;
    set(smsplFitsFig,'defaultLineLineWidth',0.4,'defaultLineMarkerSize',4, ...
        'visible','off')
    C = colororder("gem");
    for idx = 1:numNetTypes
        legendLabels = cell(2*numD,1);
        t = tiledlayout(3,2,'TileSpacing','compact');
        netLabel = sprintf('an%02i',allAstralNums(idx));
        tileIdx = 1;
        for jdx = 1:numPercTypes
            nexttile(t,tileIdx)
            hold on
            for kdx = 1:numD
                curves = allCurves{kdx};
                if isfield(curves,netLabel)
                    x = transpose(curves.(netLabel)(1,:));
                    y = transpose(curves.(netLabel)(1+jdx,:));
                    % restricting fit windows to exclude flat lines in data
                    % at p=0 and p=1
                    greaterThanZero = find(y>0,1);
                    % discard decreasing portions of connectivity perc
                    % curves that appear at low density & high astral num.
                    if greaterThanZero==1
                        lowerLim = 1;
                        while y(lowerLim+1) < y(lowerLim)
                            lowerLim = lowerLim + 1;
                        end
                    else
                        lowerLim = max([greaterThanZero - 3,1]);
                    end
                    % upper limit is just past where the curve increases 
                    % to 1 for the first time, otherwise last data point
                    if any(y==1)
                        upperLimOptions = find(y==1);
                        reachedOne = find(upperLimOptions>lowerLim,1);
                        upperLim = min([upperLimOptions(reachedOne) + 2,...
                            length(y)]);
                    else
                        upperLim = length(y);
                    end
                    x = x(lowerLim:upperLim);
                    y = y(lowerLim:upperLim);
                    % smsplFits{idx,jdx,kdx} = fit(x,y, ...
                    %     'smoothingspline','SmoothingParam',0.99);
                    smsplFits{idx,jdx,kdx} = fit(x,y, ...
                        'smoothingspline');
                    percDataUsed{idx,jdx,kdx} = cat(2,x,y);
                    l = plot(smsplFits{idx,jdx,kdx},x,y,'*');
                    [l.Color] = deal(C(kdx,:));
                    legendLabels{2*kdx - 1} = sprintf('D=%2.1f data', ...
                        Dlist(kdx));
                    legendLabels{2*kdx} = sprintf('D=%2.1f fit', ...
                        Dlist(kdx));
                end
            end
            hold off
            xscale('log')
            xlabel('Filament density [$\mu m^{-1}$]')
            ylabel('Percolation probability')
            title(percTypes{jdx})
            if jdx == 1
                legendLabels = legendLabels(~cellfun('isempty', ...
                    legendLabels));
                lg = legend(legendLabels,'FontSize',12);
                lg.Layout.Tile = 6;
            else
                legend('off')
            end
            tileIdx = tileIdx + 1;
        end
        title(t,netLabel)
        if idx == 1
            % export first page
            exportgraphics(smsplFitsFig,fullfile(saveDir,figSubDir, ...
                'percMCRG_nonparamFits_tzOnly.pdf'),'Resolution',300)
            smsplFitsFig = clf(smsplFitsFig);
        else
            % export next pages
            exportgraphics(smsplFitsFig,fullfile(saveDir,figSubDir, ...
                'percMCRG_nonparamFits_tzOnly.pdf'),'Resolution',300,'Append',true)
            smsplFitsFig = clf(smsplFitsFig);
        end
        set(smsplFitsFig,'units','inches','Position',[1,1,8.5,11], ...
            'defaultLineLineWidth',0.4,'defaultLineMarkerSize',4, ...
            'visible','off')
    end
    save(fullfile(saveDir,matFileSubDir,'percMCRG_nonparamFits_tzOnly.mat'), ...
        'smsplFits','Dlist','percDataUsed','numNetTypes','numPercTypes', ...
        'numD','allAstralNums','allCurves')
end
