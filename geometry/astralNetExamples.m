% astralNetExamples.m
% Brady Berg
clear; close all;

%% Save directories & graphics parameters

% Ubuntu path
saveDir = "~/Documents/astral-mikado-data";
% Windows path
% saveDir = "C:\Users\bcber\Documents\astral-mikado-data";

filename = "astralNetEx250423";

dataSubfolder = "mat_files";
figSubfolder = "subfigures";
figFiletype = ".pdf";

boundaryWidth = 4;
filamentWidth = 1.5;
springYellow = [255,230,0]/255;

%% Network parameters

regenerate = false;
if regenerate
    rho = 7.5;
    l = 1;
    D = 4;
    targetFilNum = rho * D^2 / l;

    an_classic = 1;
    [net_classic,crossings_classic,asters_classic] = generateAstralNetwork( ...
        targetFilNum,l,D,an_classic,false);

    an_low = 2;
    numAsters_low = round(targetFilNum/an_low,TieBreaker="tozero"); % attempt to match Python rounding behavior
    [net_low,crossings_low,asters_low] = generateAstralNetwork(numAsters_low, ...
        l,D,an_low,false);

    an_hi = 8;
    numAsters_hi = round(targetFilNum/an_hi,TieBreaker="tozero"); % attempt to match Python rounding behavior
    [net_hi,crossings_hi,asters_hi] = generateAstralNetwork(numAsters_hi, ...
        l,D,an_hi,false);

    save(fullfile(saveDir,dataSubfolder,filename + ".mat"),'rho','l','D', ...
        'targetFilNum','an_low','numAsters_low','net_low', ...
        'crossings_low','asters_low','an_hi','numAsters_hi','net_hi', ...
        'crossings_hi','asters_hi','an_classic','net_classic', ...
        'crossings_classic','asters_classic')
else
    load(fullfile(saveDir,dataSubfolder,filename + ".mat"),'rho','l','D', ...
        'targetFilNum','an_low','numAsters_low','net_low', ...
        'crossings_low','asters_low','an_hi','numAsters_hi','net_hi', ...
        'crossings_hi','asters_hi','an_classic','net_classic', ...
        'crossings_classic','asters_classic')
end

centerMarkSz = 16;
crslnkMarkSz = 12;

%% Non-astral example

anclassicEx = figure(1); clf;
set(anclassicEx,'Units','centimeters','Position',[1,1,6,6])
set(gca,'Color','k')
hold on
plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
for idx = 1:targetFilNum
    for jdx = 1:an_classic
        plot(asters_classic.centers(idx,1) + l * cos(asters_classic.orients(idx,jdx)) * (0:1), ...
            asters_classic.centers(idx,2) + l * sin(asters_classic.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',filamentWidth,'Color','w')
    end
end
% OPTION 1: PLOT CROSSLINKS (NODES) ONLY IF THEY BOUND A SPRING
% for idx = 1:size(net_classic.springs,1)
%     augNodeA = net_classic.springs(idx,1);
%     augNodeB = net_classic.springs(idx,2);
%     coords = [net_classic.augNodes(augNodeA,1:2);
%         net_classic.augNodes(augNodeB,1:2)];
%     plot(coords(:,1),coords(:,2),'.','MarkerSize',16,'LineWidth', ...
%         filamentWidth,'Color',springYellow)
% end
% OPTION 2: PLOT ALL CROSSLINKS (NODES)
plot(net_classic.nodes(:,1),net_classic.nodes(:,2),'.','MarkerSize', ...
    crslnkMarkSz,'Color',springYellow)
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - D/2)) * [-1,1] + D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - D/2)) * [-1,1] + D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
exportgraphics(anclassicEx,fullfile(saveDir,figSubfolder,"anclassicEx" + ...
    figFiletype),'ContentType','vector')

%% Example at low astral number

anlowEx = figure(2); clf;
set(anlowEx,'Units','centimeters','Position',[1,1,6,6])
set(gca,'Color','k')
hold on
plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
for idx = 1:numAsters_low
    for jdx = 1:an_low
        plot(asters_low.centers(idx,1) + l * cos(asters_low.orients(idx,jdx)) * (0:1), ...
            asters_low.centers(idx,2) + l * sin(asters_low.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',filamentWidth,'Color','w')
    end
end
% OPTION 1: plot crosslinks (nodes) only if they bound a spring
% for idx = 1:size(net_low.springs,1)
%     augNodeA = net_low.springs(idx,1);
%     augNodeB = net_low.springs(idx,2);
%     coords = [net_low.augNodes(augNodeA,1:2);
%         net_low.augNodes(augNodeB,1:2)];
%     plot(coords(:,1),coords(:,2),'.','MarkerSize',16,'LineWidth', ...
%         filamentWidth,'Color',springYellow)
%     % OPTION A: plot astral centers only if they are also a spring boundary
%     % if augNodeA <= numAsters_low && an_low >= 2
%     %     % i.e., if augNodeA is an astral center
%     %     plot(net_low.augNodes(augNodeA,1),net_low.augNodes(augNodeA,2), ...
%     %         '.r','MarkerSize',24)
%     % elseif augNodeB <= numAsters_low && an_low >= 2
%     %     % likewise for augNodeB
%     %     plot(net_low.augNodes(augNodeB,1),net_low.augNodes(augNodeB,2), ...
%     %         '.r','MarkerSize',24)
%     % end
% end
% OPTION 2: plot all crosslinks (nodes)
plot(net_low.nodes(:,1),net_low.nodes(:,2),'.','MarkerSize', ...
    crslnkMarkSz,'Color',springYellow)
% OPTION B: plot all astral centers
if an_low >= 2
    plot(asters_low.centers(:,1),asters_low.centers(:,2),'.r', ...
        'MarkerSize',centerMarkSz)
    % could use "scatter" to enable transparency for center marks
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - D/2)) * [-1,1] + D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - D/2)) * [-1,1] + D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
exportgraphics(anlowEx,fullfile(saveDir,figSubfolder,"anlowEx" + ...
    figFiletype),'ContentType','vector')

%% Example at high astral number

anhiEx = figure(3); clf;
set(anhiEx,'Units','centimeters','Position',[1,1,6,6])
set(gca,'Color','k')
hold on
plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
for idx = 1:numAsters_hi
    for jdx = 1:an_hi
        plot(asters_hi.centers(idx,1) + l * cos(asters_hi.orients(idx,jdx)) * (0:1), ...
            asters_hi.centers(idx,2) + l * sin(asters_hi.orients(idx,jdx)) * (0:1), ...
            '-','LineWidth',filamentWidth,'Color','w')
    end
end
% OPTION 1: plot crosslinks (nodes) only if they bound a spring
% for idx = 1:size(net_hi.springs,1)
%     augNodeA = net_hi.springs(idx,1);
%     augNodeB = net_hi.springs(idx,2);
%     coords = [net_hi.augNodes(augNodeA,1:2);
%         net_hi.augNodes(augNodeB,1:2)];
%     plot(coords(:,1),coords(:,2),'.','MarkerSize',16,'LineWidth', ...
%         filamentWidth,'Color',springYellow)
%     % OPTION A: plot astral centers only if they are also a spring boundary
%     % if augNodeA <= numAsters_hi && an_hi >= 2
%     %     % i.e., if augNodeA is an astral center
%     %     plot(net_hi.augNodes(augNodeA,1),net_hi.augNodes(augNodeA,2), ...
%     %         '.r','MarkerSize',24)
%     % elseif augNodeB <= numAsters_hi && an_hi >= 2
%     %     % likewise for augNodeB
%     %     plot(net_hi.augNodes(augNodeB,1),net_hi.augNodes(augNodeB,2), ...
%     %         '.r','MarkerSize',24)
%     % end
% end
% OPTION 2: plot all crosslinks (nodes)
plot(net_hi.nodes(:,1),net_hi.nodes(:,2),'.','MarkerSize', ...
    crslnkMarkSz,'Color',springYellow)
% OPTION B: plot all astral centers
if an_hi >= 2
    plot(asters_hi.centers(:,1),asters_hi.centers(:,2),'.r', ...
        'MarkerSize',centerMarkSz)
    % could use "scatter" to enable transparency for center marks
end
hold off
xlim('tight')
xl = xlim;
xlim(max(abs(xl - D/2)) * [-1,1] + D/2)
xticks([])
ylim('tight')
yl = ylim;
ylim(max(abs(yl - D/2)) * [-1,1] + D/2)
yticks([])
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
exportgraphics(anhiEx,fullfile(saveDir,figSubfolder,"anhiEx" + ...
    figFiletype),'ContentType','vector')

%% Individual asters

regenerateIndividual = false;
if regenerateIndividual
    asterTypes = 1:5;
    numAsterTypes = length(asterTypes);
    l = 1;
    centers = [2 * (1:numAsterTypes)', zeros(numAsterTypes,1)];
    orients = cell(numAsterTypes,1);
    for idx = 1:numAsterTypes
        orients{idx} = 2*pi*rand(asterTypes(idx),1);
    end
    save(fullfile(saveDir,dataSubfolder,"individualAsters.mat"), ...
        'centers','orients','l','numAsterTypes','asterTypes')
else
    load(fullfile(saveDir,dataSubfolder,"individualAsters.mat"), ...
        'centers','orients','l','numAsterTypes','asterTypes')
end

astralNumIndividual = figure(4); clf;
set(astralNumIndividual,'Units','centimeters','Position',[1,1,20,5])
set(gca,'Color','k')
axis equal
% axis padded
hold on
for idx = 1:numAsterTypes
    for jdx = 1:asterTypes(idx)
        plot(centers(idx,1) + (0:1)*l*cos(orients{idx}(jdx)), ...
            centers(idx,2) + (0:1)*l*sin(orients{idx}(jdx)), ...
            '-','LineWidth',2*filamentWidth,'Color','w')
        if asterTypes(idx) > 1
            plot(centers(idx,1),centers(idx,2),'.r','MarkerSize',2*centerMarkSz);
        end
        text(centers(idx,1),-1.5,num2str(asterTypes(idx)),'Color','w', ...
            'FontSize',18,'FontName','CMU Serif')
    end
end
text(0.75,-1.5,"$a_n = $",'Interpreter','latex','FontSize',18,'Color','w')
hold off
xticks([])
ylim([-1.75,1])
yticks([])
exportgraphics(astralNumIndividual,fullfile(saveDir,figSubfolder, ...
    "astral_num_ex" + figFiletype),'ContentType','vector')
