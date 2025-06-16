% springsPerNode.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% Define directories, network parameters

% Ubuntu path
saveDir = "~/Documents/astral-mikado-data/";
% Windows path
% saveDir = "C:\Users\bcber\Documents\astral-mikado-data";
matFileSubDir = "mat_files";
figSubDir = "exploratory_figures";

l = 1;
D = 10;
rho = 7.5;
targetFilNum = rho * D^2 / l;
an_list = (1:24)';
nNetTypes = length(an_list);
numAsters_list = round(targetFilNum ./ an_list,TieBreaker="tozero");
N = 1000;

%% Counting springs per node

recompute = false;
if recompute
    springPerNodeCounts = cell(nNetTypes,N);
    % ^entry {i,j} contains an array listing the number of springs at each node
    springPerNodeMeans = zeros(nNetTypes,N);

    p = parpool(8);
    parfor idx = 1:nNetTypes
        for jdx = 1:N
            [net,~,~] = generateAstralNetwork_mex(numAsters_list(idx), ...
                l,D,an_list(idx),false);
            numAugNodes = size(net.augNodes,1);
            theseCounts = zeros(numAugNodes,1);
            springEndpoints = [net.springs(:,1); net.springs(:,2)];
            [B,BG] = groupcounts(springEndpoints);
            theseCounts(BG) = B;
            springPerNodeCounts{idx,jdx} = cast(theseCounts,"int8");
            springPerNodeMeans(idx,jdx) = mean(theseCounts);
        end
    end
    delete(p)

    % compute mean across networks & 95% CIs
    z = norminv(0.975);
    [stdev,means_by_net_type] = std(springPerNodeMeans,0,2);
    ci_95 = stdev * z / sqrt(N);

    selected_an = [1,2,4,8,16,24];
    num_selected = length(selected_an);
    aggregatedCounts = cell(num_selected,1);
    for idx = 1:num_selected
        numAugNodes_byNet = cellfun(@length,springPerNodeCounts(selected_an(idx),:));
        aggregatedCounts{idx} = cast(zeros(sum(numAugNodes_byNet),1),"int8");
        addHere = 1;
        for jdx = 1:N
            aggregatedCounts{idx}(addHere:(addHere+numAugNodes_byNet(jdx)-1))=...
                springPerNodeCounts{selected_an(idx),jdx};
            addHere = addHere + numAugNodes_byNet(jdx);
        end
    end

    save(fullfile(saveDir,matFileSubDir,"springPerNode.mat"),'l','D', ...
        'rho','targetFilNum','an_list','nNetTypes','numAsters_list','N',...
        'springPerNodeCounts','springPerNodeMeans','stdev', ...
        'means_by_net_type','ci_95','selected_an','num_selected', ...
        'aggregatedCounts')
else
    load(fullfile(saveDir,matFileSubDir,"springPerNode.mat"),'l','D', ...
        'rho','targetFilNum','an_list','nNetTypes','numAsters_list','N',...
        'springPerNodeCounts','springPerNodeMeans','stdev', ...
        'means_by_net_type','ci_95','selected_an','num_selected', ...
        'aggregatedCounts')
end

%% Plotting

meanPlot = figure(1); clf;
set(meanPlot,'units','centimeters','Position',[1,1,10,10])
errorbar(an_list,means_by_net_type,ci_95,"vertical",'LineWidth',1)
xticks([1,4:4:24])
xlabel("Astral number $a_n$")
ylim([3,4])
ylabel({"Productive"; "segments per node"})
set(gca,"FontName","CMU Serif","FontSize",18)
exportgraphics(meanPlot,fullfile(saveDir,figSubDir,'springsPerNode.png'), ...
    'Resolution',300)

histPlots = figure(2); clf;
set(histPlots,'units','centimeters','Position',[1,1,20,16])
histPlots = tiledlayout(2,3);
for idx = 1:num_selected
    nexttile(idx)
    histogram(categorical(aggregatedCounts{idx}),'Normalization','count')
    title(sprintf("$a_n = %i$",selected_an(idx)))
    yscale('log')
    ylim([10^2,10^6])
    set(gca,'FontSize',12)
    if selected_an(idx) == 16
        labels = cellstr(string(0:16));
        [labels{2:2:16}] = deal('');
        ax = gca;
        ax.XAxis.TickLabels = labels;
    elseif selected_an(idx) == 24
        labels = cellstr(string(0:24));
        [labels{2:2:24}] = deal('');
        [labels{3:4:24}] = deal('');
        ax = gca;
        ax.XAxis.TickLabels = labels;
    end
end
xlabel(histPlots,'Productive segments per node','Interpreter','latex', ...
    'FontSize',18)
ylabel(histPlots,'Aggregate counts ($N=1000$)','Interpreter','latex','FontSize',18)
exportgraphics(histPlots,fullfile(saveDir,figSubDir,'springsPerNodeDistr.png'), ...
    'Resolution',300)
