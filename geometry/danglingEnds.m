% danglingEnds.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Parameters

rho = 7.5;
l = 1;
D = 10;
targetFilNum = rho * D^2 / l;

% Ubuntu path
saveDirectory = '~/Documents/astral-mikado-data/';
% Windows path
% saveDirectory = 'C:\Users\bcber\Documents\astral-mikado-data\';
datafilename = 'danglingEnds_250319';
histfilename = 'endsHist_250319';
summaryfilename = 'endsSummary_250319';
filetype = '.pdf';

uciBlue = [0,100,164]/255;
uciGold = [255,210,0]/255;
uciOrange = [247,141,45]/255;
fname = 'CMU Serif';

%% Generating networks

generateNew = false;
if generateNew
    numNets = 10;   % number of networks to generate per astral number
    nFilPerAsterList = 1:24;    % astral number
    numNetTypes = length(nFilPerAsterList);
    endsLengths = struct();
    % mean of dangling ends distribution (over [0,l])
    endsMeans = zeros(numNetTypes,numNets);
    % mean of dangling ends distribution truncated to [0,l)
    endsMeansTrunc = zeros(numNetTypes,numNets);
    % proportion of total length NOT in dangling ends
    useFraction = zeros(numNetTypes,numNets);
    labelPat = "an%02i";
    for idx = 1:numNetTypes
        astralNum = nFilPerAsterList(idx);
        netLabel = sprintf(labelPat,astralNum);
        numAsters = round(targetFilNum/astralNum,TieBreaker="tozero");
        % (attempt to match Python rounding behavior)
        % preallocate arrays to store dangling ends data
        if astralNum == 1
            endsPerNet = numAsters * 2; % "left" and "right" ends
            totalLength = numAsters * l;
            endsLengths.(netLabel) = zeros(endsPerNet * numNets,1);
            for n = 1:numNets
                [thisNet,~,~] = generateAstralNetwork(numAsters,l,D, ...
                    astralNum,false);
                putHere = (1:endsPerNet) + (n-1)*endsPerNet;
                endsVec = reshape(thisNet.ends,[endsPerNet,1]);
                endsLengths.an01(putHere) = endsVec;
                endsMeans(idx,n) = mean(endsVec);
                endsMeansTrunc(idx,n) = mean(endsVec(endsVec < l));
                useFraction(idx,n) = 1 - sum(endsVec)/totalLength;
            end
        elseif astralNum >=2
            endsPerNet = numAsters * astralNum;
            totalLength = numAsters * astralNum * l;
            endsLengths.(netLabel) = zeros(endsPerNet * numNets,1);
            for n = 1:numNets
                [thisNet,~,~] = generateAstralNetwork(numAsters,l,D, ...
                    astralNum,false);
                putHere = (1:endsPerNet) + (n-1)*endsPerNet;
                % record only the "distal" dangling ends (others are 0)
                endsVec = thisNet.ends(:,2);
                endsLengths.(netLabel)(putHere) = endsVec;
                endsMeans(idx,n) = mean(endsVec);
                endsMeansTrunc(idx,n) = mean(endsVec(endsVec < l));
                useFraction(idx,n) = 1 - sum(endsVec)/totalLength;
            end
        end
    end
    stdOfMeans = std(endsMeans,0,2);
    stdOfMeansTrunc = std(endsMeansTrunc,0,2);
    save(fullfile(saveDirectory,'mat_files',datafilename),'numNets', ...
        'nFilPerAsterList','numNetTypes','endsLengths','endsMeans', ...
        'stdOfMeans','endsMeansTrunc','stdOfMeansTrunc','rho','l','D', ...
        'targetFilNum','labelPat','useFraction')
else
    load(fullfile(saveDirectory,'mat_files',datafilename),'numNets', ...
        'nFilPerAsterList','numNetTypes','endsLengths','endsMeans', ...
        'stdOfMeans','endsMeansTrunc','stdOfMeansTrunc','rho','l','D', ...
        'targetFilNum','labelPat','useFraction')
end

%% Dangling ends distributions

endsHists = figure(1);
set(endsHists,'units','centimeters','Position',[1,1,18,15],'visible','on')
endsHists = tiledlayout(2,2,'TileSpacing','compact');
tilesPerPage = 4;

selectedAstralNums = [1,4,8,16];
for idx = 1:length(selectedAstralNums)
    % nexttile(endsHists,mod(idx,tilesPerPage) + tilesPerPage * (mod(idx,tilesPerPage)==0))
    nexttile(idx)
    astralNum = selectedAstralNums(idx);
    netLabel = sprintf(labelPat,astralNum);
    theseEnds = endsLengths.(netLabel);
    hold on
    h = histogram(theseEnds,'Normalization','count', ...
        'BinWidth',0.1);
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
    hTrunc = histogram(theseEnds(theseEnds<l),'Normalization', ...
        'count','BinWidth',0.1);
    hTrunc.FaceColor = uciBlue;
    hTrunc.FaceAlpha = 1;
    title(sprintf('$a_n = %i$',astralNum))
    set(gca,'FontName',fname,'FontSize',18)
    hold off
    if astralNum == 1
        ylim([0,6000])
    else
        ylim([0,2500])
    end
end
xlabel(endsHists,'Dangling end length [$\mu\rm{m}$]','Interpreter', ...
    'latex','FontSize',18)
ylabel(endsHists,'Counts','FontName',fname,'FontSize',18)
lg = legend('All ends','Ends $<\ell$','Interpreter','latex');
lg.Layout.Tile = 'east';
exportgraphics(endsHists,fullfile(saveDirectory,'exploratory_figures',...
            [histfilename,filetype]),'ContentType','vector')

%% Mean and total dangling ends

linewidth = 1.5;
z = norminv(0.975);
% useFraction statistics
[useSTD,useMean] = std(useFraction,0,2);

fig2 = figure(2);
colororder([uciBlue; 0.8500 0.3250 0.0980])
yyaxis left
hold on
l1 = errorbar(nFilPerAsterList,mean(endsMeans,2), ...
    z*stdOfMeans/sqrt(numNets),'LineStyle','-','Color','k', ...
    'LineWidth',linewidth);
l2 = errorbar(nFilPerAsterList,mean(endsMeansTrunc,2), ...
    z*stdOfMeansTrunc/sqrt(numNets),'LineStyle','-','Color',uciBlue, ...
    'LineWidth',linewidth);
hold off
yl = ylim;
ylim([0,yl(2)])
ylabel('Dangling end length [$\mu\rm{m}$]')
yyaxis right
l3 = errorbar(nFilPerAsterList,useMean,z*useSTD/sqrt(numNets), ...
    'LineWidth',linewidth);
ylim([0,1])
ylabel('Use fraction')
xlabel('Astral number $a_n$')
xlim([nFilPerAsterList(1)-0.25,nFilPerAsterList(end)+0.25])
xticks([1,4:4:24])
legend([l1,l2,l3],{'All ends','Ends $<\ell$','Use fraction'}, ...
    'Location','south')
set(gca,'FontName',fname,'FontSize',18)
exportgraphics(fig2,fullfile(saveDirectory,'exploratory_figures', ...
    [summaryfilename,filetype]),'ContentType','vector')


