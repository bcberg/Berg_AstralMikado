% test_end_force.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

newReports = true;  % only parse reports if data is new
if newReports
    dir = "~/Documents/AstralMikadoCYM/runs/test_end_force/";
    filePattern = strcat(dir, "run%04i", "/fiber_ends.txt");
    varNames = {'class', 'identity', 'length', 'stateM', 'posMX', 'posMY', ...
        'dirMX', 'dirMY', 'stateP', 'posPX', 'posPY', 'dirPX', 'dirPY'};
    opts = fixedWidthImportOptions('VariableNames',varNames,...
        'VariableWidths',10*ones(1,length(varNames)),'DataLines',[7,7],...
        'SelectedVariableNames',{'posPX','posPY'});
    opts = setvartype(opts,{'posPX','posPY'},'double');

    numRuns = 101;
    forceVals = (0:0.1:10)';   % ensure match with values in test_end_force.cym.tpl
    finalPos = zeros(length(forceVals),2);
    for idx = 1:numRuns
        filename = sprintf(filePattern,idx-1);
        thisPos = readmatrix(filename,opts);
        finalPos(idx,:) = thisPos;
    end
    save('~/Documents/AstralMikadoCYM/data/test_end_force.mat','forceVals','finalPos')
    clear
end
% for Ubuntu
basePath = '~/Documents/AstralMikadoCYM/data';
% for Windows
% basePath = 'C:\Users\bcber\Documents\AstralMikadoCYM\data';
filename = 'test_end_force.mat';
load(fullfile(basePath,filename),'forceVals','finalPos')

%% Plot deflection vs. end-loaded-beam analytics

R = 1e3;    % bending rigidity [pN*um^2]
L = 10;     % length [um]
% deflection from initial position x=0 [um]
prediction = forceVals * L^3 / (3 * R);
simulated = finalPos(:,1);

fig1 = figure(1);
hold on;
plot(forceVals,prediction,'-b','LineWidth',2)
plot(forceVals,simulated,'k*','MarkerSize',6)
hold off;
xlabel('Applied force [pN]')
ylabel('Fiber deflection [{\mu}m]')
legend('Cantilever beam','Simulated','Location','southeast')
title('Fiber rigidity 1e3, base rigidity 1e6 [pN\cdot{\mu}m^2]')
exportgraphics(fig1,fullfile(basePath,'test_end_force.png'), ...
    'Resolution',300)