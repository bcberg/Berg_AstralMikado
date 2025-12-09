% percCrit_eqSpaced.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Parameters

% Ubuntu path
dataDirectory = "/home/bcberg/Documents/astral-mikado-data/mat_files";
figDirectory = "/home/bcberg/Documents/astral-mikado-data/exploratory_figures";
% Windows path
% dataDirectory = "C:\Users\bcber\Documents\astral-mikado-data\mat_files";
% figDirectory = "C:\Users\bcber\Documents\astral-mikado-data\exploratory_figures";

D = 10;
l = 1;
filename = sprintf("percProbs_eqSpaced_l%02i_D%02i",l,D);
load(fullfile(dataDirectory,filename + ".mat"),'actualDensities', ...
    'percProbs','curves','densityRange','densSpec','numNetTypes', ...
    'numUniqueDensVals_byrow','astralNumList')

percTypes = {"Top-to-Bottom (TB)", "Left-to-Right (LR)", "TB OR LR", ...
    "TB AND LR", "Single connected component"};
numPercTypes = length(percTypes); % i.e., = 5

%% Fitting smoothing splines

% current scheme: smoothing spline with automatic smoothing parameter
makeSmsplFitsFig = true;
if makeSmsplFitsFig
    smsplFits = cell(numNetTypes,numPercTypes);
    percDataUsed = cell(numNetTypes,numPercTypes);
    smsplFitsFig = figure('units','inches','Position',[1,1,8.5,11]); clf;
    set(smsplFitsFig,'defaultLineLineWidth',0.4,'defaultLineMarkerSize',4, ...
        'visible','off')
    C = colororder("gem");
    legendLabels = {sprintf('l/D=%1.2f data',l/D), ...
        sprintf('l/D=%1.2f fit',l/D)};
    for idx = 1:numNetTypes
        t = tiledlayout(3,2,'TileSpacing','compact');
        netLabel = sprintf('an%02d',astralNumList(idx));
        tileIdx = 1;
        for jdx = 1:numPercTypes
            nexttile(t,tileIdx)
            x = transpose(curves.(netLabel)(1,:));
            y = transpose(curves.(netLabel)(1+jdx,:));
            % restricting fit windows to exclude flat lines in data at
            % p=0,1
            greaterThanZero = find(y>0,1);
            % discard decreasing portions of connectivity perc curves that
            % may appear at low density and high astral number
            if greaterThanZero == 1
                lowerLim = 1;
                while y(lowerLim+1) < y(lowerLim)
                    lowerLim = lowerLim + 1;
                end
            else
                lowerLim = max([greaterThanZero - 3,1]);
            end
            % upper limit is just past where the curve increases to 1 for
            % the first time, otherwise last data point
            if any(y==1)
                upperLimOptions = find(y==1);
                reachedOne = find(upperLimOptions>lowerLim,1);
                upperLim = min([upperLimOptions(reachedOne) + 2,length(y)]);
            else
                upperLim = length(y);
            end
            x = x(lowerLim:upperLim);
            y = y(lowerLim:upperLim);
            smsplFits{idx,jdx} = fit(x,y,'smoothingspline');
            percDataUsed{idx,jdx} = cat(2,x,y);
            line = plot(smsplFits{idx,jdx},x,y,'*');
            [line.Color] = deal(C(1,:));   % in past, color distinguished multiple l/D vals
            xscale('log')
            xlabel('Filament density [$\mu m^{-1}$]')
            ylabel('Percolation probability')
            title(percTypes{jdx})
            if jdx == 1
                lg = legend(legendLabels,'FontSize',12);
                lg.Layout.Tile = 6;
            else
                legend('off')
            end
            tileIdx = tileIdx + 1;
        end
        title(t,netLabel)
        if idx == 1
            exportgraphics(smsplFitsFig,fullfile(figDirectory, ...
                'percFits_nonparam_tzOnly_eqSpaced.pdf'),'Resolution',300)
            smsplFitsFig = clf(smsplFitsFig);
        else
            exportgraphics(smsplFitsFig,fullfile(figDirectory, ...
                'percFits_nonparam_tzOnly_eqSpaced.pdf'),'Resolution', ...
                300,'Append',true)
            smsplFitsFig = clf(smsplFitsFig);
        end
        set(smsplFitsFig,'defaultLineLineWidth',0.4,'defaultLineMarkerSize',4, ...
        'visible','off')
    end
    % save(fullfile(dataDirectory,'percFits_nonparam_tzOnly_eqSpaced.mat'), ...
    %     'smsplFits','l','D','percDataUsed')
end

%% Estimating critical densities

percTypeChoice = [1,5];
numPercTypesChosen = length(percTypeChoice);
to50 = zeros(numPercTypesChosen,numNetTypes);
opts = optimoptions("fsolve","Display","off");
for idx = 1:numPercTypesChosen
    for jdx = 1:numNetTypes
        thisCurve = smsplFits{jdx,percTypeChoice(idx)};
        if percTypeChoice(idx) == 5 && jdx == 1
            initGuess = 50;
        else
            initGuess = 10;
        end
        to50(idx,jdx) = fsolve(@(x) thisCurve(x) - 0.5, initGuess, opts);
    end
end
% save estimates of critical densities to streamline plotting later
save(fullfile(dataDirectory,'percFits_nonparam_tzOnly_eqSpaced.mat'), ...
    'smsplFits','l','D','percDataUsed','percTypeChoice', ...
    'numPercTypesChosen','numNetTypes','astralNumList','to50')

threshPlotEQ = figure(2); clf;
set(threshPlotEQ,'units','centimeters','Position',[1,1,15,12], ...
    'defaultLineLineWidth',2.5)
hold on
for idx = 1:numPercTypesChosen
    plot(astralNumList',to50(idx,:),'Color',C(idx,:),'LineStyle','-')
end
yline(7.5,'-','7.5','Interpreter','latex','LineWidth',2, ...
    'FontSize',16,'LabelVerticalAlignment','middle')
hold off
xticks([1,4:4:24])
xlim([1,24])
xlabel('Astral number $a_n$')
ylim([3,51])
yticks(2.^(2:5))
yticklabels(arrayfun(@(x) sprintf("$2^%i$",x),2:5))
yscale('log')
ylabel('Critical density $\rho_c$ [$\mu{\rm m^{-1}}$]')
set(gca,'FontSize',22,'FontName','CMU Serif')
legend({'Spanning','Connectivity'},'Location','northeast', ...
    'FontSize',16,'FontName','CMU Serif')
exportgraphics(threshPlotEQ,fullfile(figDirectory,['percCrit_eqSpaced', ...
    '.pdf']),'ContentType','Vector')

%% Correlation between deltaRho and moduli

% shear modulus sweep, equally spaced filaments
load(fullfile(dataDirectory,'rf_shear_eq_actin251001_l0.1_D1.mat'), ...
    'shearModuli','forceVals','nFilPerAsterList','nNetTypes','numRep', ...
    'dens')
sweep_s_eq = struct('forceVals',forceVals,'nFilPerAsterList', ...
    nFilPerAsterList,'nNetTypes',nNetTypes,'numRep',numRep,'dens',dens, ...
    'len_fil',0.1,'D',1);
% disps_s_eq = disps;
moduli_s_eq = shearModuli;
clear('shearModuli','forceVals','nFilPerAsterList','nNetTypes','numRep',...
    'dens','len_fil','D')

connRhoc = transpose(to50(2,:));
deltaRho = sweep_s_eq.dens - 10 * connRhoc; % mult by 10 since rhoc were
% estimated from l = 1 & D = 10 but simulation done at l = 0.1 and D = 1
% (recall rho = N * l / D^2 where N is number of filaments)

corrPlot_eq = figure(3); clf;
set(corrPlot_eq,'units','centimeters','Position',[1,1,18,12])
scatter(repmat(deltaRho,[1,sweep_s_eq.numRep]),moduli_s_eq,'Marker','*', ...
    'MarkerEdgeColor','k')
yscale('log')
xlabel('Density rel. to critical density $\rho - \rho^{\rm conn}_c(a_n)$')
ylabel('Shear modulus $G$ [${\rm pN}\cdot \mu{\rm m^{-1}}$]')
title('equally spaced filaments')
