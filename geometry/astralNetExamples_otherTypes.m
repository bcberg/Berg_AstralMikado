% astralNetExamples.m
% Brady Berg
clear; close all;

%% Save directories & graphics parameters

% Ubuntu path
saveDir = "~/Documents/astral-mikado-data";
% Windows path
% saveDir = "C:\Users\bcber\Documents\astral-mikado-data";

filename = "astralNetEx_otherTypes250910";

dataSubfolder = "mat_files";
figSubfolder = "exploratory_figures";
figFiletype = ".pdf";

boundaryWidth = 4;
filamentWidth = 1.5;
springYellow = [255,230,0]/255;
centerMarkSz = 16;
crslnkMarkSz = 12;

%% Network parameters

regenerate = false;
if regenerate
    rho = 7.5;
    l = 1;
    D = 5;
    targetFilNum = rho * D^2 / l;

    % equally spaced filaments (radially)
    an_list_EQ = 1:6;
    numAsters_list_EQ = round(targetFilNum ./ an_list_EQ, TieBreaker="tozero");
    % ^ attempt to match Python rounding behavior

    for idx = 1:length(an_list_EQ)
        netLabel = sprintf('an%02d',an_list_EQ(idx));
        [netEQ.(netLabel),crossingsEQ.(netLabel),astersEQ.(netLabel)] = ...
            generateAstralNetwork_eqSpaced(numAsters_list_EQ(idx),l,D,...
            an_list_EQ(idx),false);
    end

    % exponentially distributed filament lengths
    an_list_EXP = [1,2,4,8,12,16];
    numAsters_list_EXP = round(targetFilNum ./ an_list_EXP, TieBreaker="tozero");
    % ^ attempt to match Python rounding behavior

    for idx = 1:length(an_list_EXP)
        netLabel = sprintf('an%02d',an_list_EXP(idx));
        [netEXP.(netLabel),crossingsEXP.(netLabel),astersEXP.(netLabel)] = ...
            generateAstralNetwork_exp(numAsters_list_EXP(idx),l,D,...
            an_list_EXP(idx),false);
    end

    save(fullfile(saveDir,dataSubfolder,filename + ".mat"),'rho','l',...
        'D','targetFilNum','an_list_EQ','numAsters_list_EQ','netEQ',...
        'crossingsEQ','astersEQ','an_list_EXP','numAsters_list_EXP',...
        'netEXP','crossingsEXP','astersEXP')
else
    load(fullfile(saveDir,dataSubfolder,filename + ".mat"),'rho','l',...
        'D','targetFilNum','an_list_EQ','numAsters_list_EQ','netEQ',...
        'crossingsEQ','astersEQ','an_list_EXP','numAsters_list_EXP',...
        'netEXP','crossingsEXP','astersEXP')
end

%% Equally spaced filaments

equalExamples = figure(1); clf;
set(equalExamples,'Units','centimeters','Position',[1,1,20,30])
equalExamples = tiledlayout(3,2,'TileSpacing','compact');

for netIdx = 1:length(an_list_EQ)
    nexttile(netIdx)
    netLabel = sprintf('an%02d',an_list_EQ(netIdx));
    set(gca,'Color','k')
    hold on
    plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
    plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
    for idx = 1:numAsters_list_EQ(netIdx)
        for jdx = 1:an_list_EQ(netIdx)
            plot(astersEQ.(netLabel).centers(idx,1) + ...
                l * cos(astersEQ.(netLabel).orients(idx,jdx)) * (0:1),...
                astersEQ.(netLabel).centers(idx,2) + ...
                l * sin(astersEQ.(netLabel).orients(idx,jdx)) * (0:1),...
                '-','LineWidth',filamentWidth,'Color','w')
        end
    end
    % plot all crosslinks
    plot(netEQ.(netLabel).nodes(:,1),netEQ.(netLabel).nodes(:,2),'.',...
        'MarkerSize',crslnkMarkSz,'Color',springYellow)
    % plot astral centers, if applicable
    if an_list_EQ(netIdx) >= 2
        plot(astersEQ.(netLabel).centers(:,1),...
            astersEQ.(netLabel).centers(:,2),'.r','MarkerSize',centerMarkSz)
        % could use "scatter" to enable transparency for center marks
    end
    hold off
    title(sprintf('$a_n = %d$',an_list_EQ(netIdx)),'Interpreter','latex')
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
end
exportgraphics(equalExamples,fullfile(saveDir,figSubfolder,"eqSpacedEx"+...
    figFiletype),'ContentType','vector')

%% Equally spaced: spanning example

filename = 'spanEx_eq';
goodExample = false;
timesTried = 0;
rho = 5;
l = 1;
D = 5;
targetFilNum = rho * D^2 / l;
astralNum = 5;
numAsters = round(targetFilNum/astralNum,TieBreaker="tozero");
maxAttempts = 100;
while ~goodExample && timesTried < maxAttempts
    [net_eq_span,crossings_eq_span,asters_eq_span] = ...
        generateAstralNetwork_eqSpaced(numAsters,l,D,astralNum,false);
    [percTF,~] = percCheck(crossings_eq_span,net_eq_span.nodes,D);
    timesTried = timesTried + 1;
    if percTF(1) && ~percTF(5)  % spanning, but not fully connected
        goodExample = true;
    end
end
if goodExample
    save(fullfile(saveDir,dataSubfolder,[filename,'.mat']),'astralNum', ...
        'rho','l','D','targetFilNum','numAsters','net_eq_span', ...
        'crossings_eq_span','asters_eq_span','percTF')
    spanEx_eq = figure(99);
    set(spanEx_eq,'units','centimeters','Position',[1,1,10,10])
    set(gca,'Color','k')
    hold on
    plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
    plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
    for idx = 1:numAsters
        for jdx = 1:astralNum
            plot(asters_eq_span.centers(idx,1) + ...
                l * cos(asters_eq_span.orients(idx,jdx)) * (0:1),...
                asters_eq_span.centers(idx,2) + ...
                l * sin(asters_eq_span.orients(idx,jdx)) * (0:1),...
                '-','LineWidth',filamentWidth,'Color','w')
        end
    end
    % plot all crosslinks
    plot(net_eq_span.nodes(:,1),net_eq_span.nodes(:,2),'.',...
        'MarkerSize',crslnkMarkSz,'Color',springYellow)
    % plot astral centers, if applicable
    if astralNum >= 2
        plot(asters_eq_span.centers(:,1),...
            asters_eq_span.centers(:,2),'.r','MarkerSize',centerMarkSz)
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
elseif timesTried == maxAttempts
    fprintf('Unable to find satisfactory spanning example (%d attempts)\n', ...
        maxAttempts)
end


%% Equally spaced: single connected component example

filename = 'connEx_eq';
goodExample = false;
timesTried = 0;
rho = 5;
l = 1;
D = 5;
targetFilNum = rho * D^2 / l;
astralNum = 5;
numAsters = round(targetFilNum/astralNum,TieBreaker="tozero");
maxAttempts = 100;
while ~goodExample && timesTried < maxAttempts
    [net_eq_conn,crossings_eq_conn,asters_eq_conn] = ...
        generateAstralNetwork_eqSpaced(numAsters,l,D,astralNum,false);
    [percTF,~] = percCheck(crossings_eq_conn,net_eq_conn.nodes,D);
    timesTried = timesTried + 1;
    if ~percTF(1) && percTF(5)  % single connected component, not spanning
        goodExample = true;
    end
end
if goodExample
    save(fullfile(saveDir,dataSubfolder,[filename,'.mat']),'astralNum', ...
        'rho','l','D','targetFilNum','numAsters','net_eq_conn', ...
        'crossings_eq_conn','asters_eq_conn','percTF')
    spanEx_eq = figure(99);
    set(spanEx_eq,'units','centimeters','Position',[1,1,10,10])
    set(gca,'Color','k')
    hold on
    plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
    plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
    for idx = 1:numAsters
        for jdx = 1:astralNum
            plot(asters_eq_conn.centers(idx,1) + ...
                l * cos(asters_eq_conn.orients(idx,jdx)) * (0:1),...
                asters_eq_conn.centers(idx,2) + ...
                l * sin(asters_eq_conn.orients(idx,jdx)) * (0:1),...
                '-','LineWidth',filamentWidth,'Color','w')
        end
    end
    % plot all crosslinks
    plot(net_eq_conn.nodes(:,1),net_eq_conn.nodes(:,2),'.',...
        'MarkerSize',crslnkMarkSz,'Color',springYellow)
    % plot astral centers, if applicable
    if astralNum >= 2
        plot(asters_eq_conn.centers(:,1),...
            asters_eq_conn.centers(:,2),'.r','MarkerSize',centerMarkSz)
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
elseif timesTried == maxAttempts
    fprintf('Unable to find satisfactory spanning example (%d attempts)\n', ...
        maxAttempts)
end

%% Exponentially distributed filament lengths

expExamples = figure(2); clf;
set(expExamples,'Units','centimeters','Position',[1,1,20,30])
expExamples = tiledlayout(3,2,'TileSpacing','compact');

for netIdx = 1:length(an_list_EXP)
    nexttile(netIdx)
    netLabel = sprintf('an%02d',an_list_EXP(netIdx));
    set(gca,'Color','k')
    hold on
    plot([0,D],[0,0],'-','LineWidth',4,'Color','blue')
    plot([0,D],[D,D],'-','LineWidth',4,'Color','blue')
    for idx = 1:numAsters_list_EXP(netIdx)
        for jdx = 1:an_list_EXP(netIdx)
            plot(astersEXP.(netLabel).centers(idx,1) + ...
                astersEXP.(netLabel).lengths(idx,jdx) * ...
                cos(astersEXP.(netLabel).orients(idx,jdx)) * (0:1),...
                astersEXP.(netLabel).centers(idx,2) + ...
                astersEXP.(netLabel).lengths(idx,jdx) * ...
                sin(astersEXP.(netLabel).orients(idx,jdx)) * (0:1),...
                '-','LineWidth',filamentWidth,'Color','w')
        end
    end
    % plot all crosslinks
    plot(netEXP.(netLabel).nodes(:,1),netEXP.(netLabel).nodes(:,2),'.',...
        'MarkerSize',crslnkMarkSz,'Color',springYellow)
    % plot astral centers, if applicable
    if an_list_EXP(netIdx) >= 2
        plot(astersEXP.(netLabel).centers(:,1),...
            astersEXP.(netLabel).centers(:,2),'.r','MarkerSize',centerMarkSz)
        % could use "scatter" to enable transparency for center marks
    end
    hold off
    title(sprintf('$a_n = %d$',an_list_EXP(netIdx)),'Interpreter','latex')
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
end
exportgraphics(expExamples,fullfile(saveDir,figSubfolder,"expEx"+...
    figFiletype),'ContentType','vector')

%% Individual asters (equally spaced)

regenerateIndividual = false;
if regenerateIndividual
    asterTypes = 1:5;
    numAsterTypes = length(asterTypes);
    l = 1;
    centers = [2 * (1:numAsterTypes)', zeros(numAsterTypes,1)];
    orients = cell(numAsterTypes,1);
    for idx = 1:numAsterTypes
        astralNum = asterTypes(idx);
        spacing = (2 * pi / astralNum) * (0:(astralNum-1));
        rotation = repmat(2 * pi * rand,[1,astralNum]);
        orients{idx} = spacing + rotation;
    end
    save(fullfile(saveDir,dataSubfolder,"individualAsters_eq.mat"), ...
        'centers','orients','l','numAsterTypes','asterTypes')
else
    load(fullfile(saveDir,dataSubfolder,"individualAsters_eq.mat"), ...
        'centers','orients','l','numAsterTypes','asterTypes')
end

astralNumIndividual_eq = figure(4); clf;
set(astralNumIndividual_eq,'Units','centimeters','Position',[1,1,20,5])
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
xlim([0.5,centers(end,1)+1.5])
xticks([])
ylim([-1.75,1.25])
yticks([])
exportgraphics(astralNumIndividual_eq,fullfile(saveDir,figSubfolder, ...
    "astral_num_ex_eq" + figFiletype),'ContentType','vector')
