% astralNetTestVisuals.m
% Brady Berg
clear; close all;
format long; format compact;
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

%% testing generation

% tiny test
% astralNum = 1;
% numAsters = 60;
% l = 5;
% D = 10;
% [network,crossings,asters] = generateAstralNetwork(numAsters,l,D, ...
%     astralNum,false);

% astralNum = 5;
% numAsters = 12;
% l = 5;
% D = 10;
% [network,crossings,asters] = generateAstralNetwork(numAsters,l,D, ...
%     astralNum,false);

% like cytosim
astralNum = 5;
rho = 1;
l = 5;
D = 50;
targetFilNum = rho * D^2 / l;
numAsters = round(targetFilNum/astralNum,TieBreaker="tozero");  % attempt to match Python rounding behavior;
[network,crossings,asters] = generateAstralNetwork(numAsters,l,D, ...
    astralNum,false);

%% testing connCheck

% [connTF, connStats] = connCheck(crossings);

%% testing percCheck

[percTF, percStats] = percCheck(crossings,network.nodes,D);

%% visualizing

domainBox = [D, D; 0, D; 0, 0; D, 0; D, D];
% aster plot
figure(1)
hold on
for idx = 1:numAsters
    if astralNum >= 2
        plot(asters.centers(idx,1),asters.centers(idx,2),'or')
    end
    for jdx = 1:astralNum
        plot(asters.centers(idx,1) + l * cos(asters.orients(idx,jdx)) * (0:1), ...
            asters.centers(idx,2) + l * sin(asters.orients(idx,jdx)) * (0:1), ...
            '-k')
    end
end
plot(domainBox(:,1),domainBox(:,2),'--b','LineWidth',0.75)
hold off
percStatus = sprintf('[percTB, percLR] = [%i,%i]',percTF(1),percTF(2));
title(percStatus)

% highlight nodes
figure(2)
hold on
for idx = 1:numAsters
    if astralNum >= 2
        plot(asters.centers(idx,1),asters.centers(idx,2),'or')
    end
    for jdx = 1:astralNum
        plot(asters.centers(idx,1) + l * cos(asters.orients(idx,jdx)) * (0:1), ...
            asters.centers(idx,2) + l * sin(asters.orients(idx,jdx)) * (0:1), ...
            '-k')
    end
end
plot(network.nodes(:,1),network.nodes(:,2),'*m')
plot(domainBox(:,1),domainBox(:,2),'--b','LineWidth',0.75)
hold off
title(percStatus)

% draw spring network
figure(3)
hold on
for idx = 1:size(network.springs,1)
    augNodeA = network.springs(idx,1);
    augNodeB = network.springs(idx,2);
    coords = [network.augNodes(augNodeA,1:2);
        network.augNodes(augNodeB,1:2)];
    plot(coords(:,1), coords(:,2),'.k-')
    if augNodeA <= numAsters && astralNum >= 2
        % i.e., if augNodeA is an astral center
        plot(network.augNodes(augNodeA,1),network.augNodes(augNodeA,2),'or')
    elseif augNodeB <= numAsters && astralNum >= 2
        % likewise for augNodeB
        plot(network.augNodes(augNodeB,1),network.augNodes(augNodeB,2),'or')
    end
end
plot(domainBox(:,1),domainBox(:,2),'--b','LineWidth',0.75)
hold off
title(percStatus)

% dangling end distribution
% figure(4)
% if astralNum == 1
%     histogram(reshape(network.ends,[2*numAsters,1]))
% elseif astralNum >= 2
%     histogram(network.ends(:,2))
% end
% xlabel('Dangling end length')
% ylabel('Counts')
