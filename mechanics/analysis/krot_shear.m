% krot_shear.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "krot_shear251202";
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';
dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';

newReports = false;
if newReports
    filtering.TF = true;
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/krot_shear/";
    
    % ensure the following parameters match those used in
    % krot_shear.cym.tpl and its driver (for numRep)
    dens = 75;  % N * len_fil / A
    len_fil = 0.1;
    D = 1; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = ([1:6,8:4:20])'; % "astral number"
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    forceVals = 0.05 * (0:5)';
    nForces = length(forceVals);
    numRep = 10;
    krot_list = 50 * (0:8);   % [pN*um/rad] (units according to Cytosim docs)
    nKrot = length(krot_list);

    runsPerKrot = nNetTypes * nForces;
    runsPerRep = nKrot * nNetTypes * nForces;
    shearModuli = zeros(nNetTypes,nKrot,numRep);
    percStatus = zeros(nNetTypes,nKrot,numRep,2);
    for repIdx = 1:numRep
        for kIdx = 1:nKrot
            for netIdx = 1:nNetTypes
                nFilPerAster = nFilPerAsterList(netIdx);
                netLabel = sprintf('an%02i',nFilPerAster);
                if repIdx == 1
                    initCoMs.(netLabel) = zeros(nForces,2,nKrot,numRep);
                    finalCoMs.(netLabel) = zeros(nForces,2,nKrot,numRep);
                    disps.(netLabel) = zeros(nForces,2,nKrot,numRep);
                    lineFits.(netLabel) = cell(nKrot,numRep);
                    fitStats.(netLabel) = cell(nKrot,numRep);
                end
                cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
                    'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);
                unshiftedRunVals = (netIdx-1)*nForces:(netIdx*nForces - 1);
                krotShift = (kIdx - 1) * runsPerKrot;
                repShift = (repIdx - 1) * runsPerRep;
                trueRunVals = unshiftedRunVals + krotShift + repShift;
                % check network percolation status
                [percTF,~,~] = checkLinks(dir,trueRunVals(1),cytoparams);
                filtering.spanCheck = percTF(1);
                percStatus(netIdx,kIdx,repIdx,:) = reshape(percTF,[1,1,1,2]);
                % analyze network position data
                [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,trueRunVals, ...
                    cytoparams);
                initCoMs.(netLabel)(:,:,kIdx,repIdx) = thisInitCoM;
                finalCoMs.(netLabel)(:,:,kIdx,repIdx) = thisFinalCoM;
                disps.(netLabel)(:,:,kIdx,repIdx) = thisFinalCoM ...
                    - thisInitCoM;
                % estimate a modulus
                [thisModulus,linefit,gof] = getModulus(thisInitCoM, ...
                    thisFinalCoM,forceVals,"shear",filtering);
                shearModuli(netIdx,kIdx,repIdx) = thisModulus;
                linefits.(netLabel){kIdx,repIdx} = linefit;
                fitStats.(netLabel){kIdx,repIdx} = gof;
            end
        end
    end
    save(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
        'D','nFil','nFilPerAsterList','forceVals','numRep','nKrot', ...
        'nNetTypes','krot_list','initCoMs','finalCoMs','disps', ...
        'linefits','fitStats','shearModuli','percStatus','filtering')
    clearvars -except simName saveDir dataSubdir figSubdir
end
load(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
    'D','nFil','nFilPerAsterList','forceVals','numRep','nKrot', ...
    'nNetTypes','krot_list','initCoMs','finalCoMs','disps', ...
    'linefits','fitStats','shearModuli','percStatus','filtering')

lgray = [197,190,181]/255;

[stdevMod,meanMod] = std(shearModuli,0,3);
z = norminv(0.975);

%% Modulus vs. krot, separate axes for each astral number

fig1 = figure(1); clf;
set(fig1,'units','inches','Position',[0,0,8.5,11])
fig1 = tiledlayout(2,2);
tileIdx = 1;
for netIdx = 1:nNetTypes
    nexttile(tileIdx)
    hold on
    for kIdx = 1:nKrot
        plot(krot_list(kIdx),squeeze(shearModuli(netIdx,kIdx,:)), ...
            'Marker','o','MarkerSize',3,'MarkerFaceColor',lgray, ...
            'MarkerEdgeColor',lgray)
    end
    errorbar(krot_list,meanMod(netIdx,:), ...
        (z/sqrt(numRep))*stdevMod(netIdx,:),'-ob','LineWidth',1.5)
    hold off
    xlabel('Angular stiffness at centers [pN*{\mu}m/rad]')
    ylabel('Shear modulus [pN/{\mu}m]')
    title(sprintf('a_n = %i',nFilPerAsterList(netIdx)))
    xscale('linear')

    if netIdx == 4
        exportgraphics(fig1,fullfile(saveDir,"exploratory_figures", ...
            simName + ".pdf"),"ContentType",'vector');
        clf(fig1)
        set(fig1,'units','inches','Position',[0,0,8.5,11])
        fig1 = tiledlayout(2,2);
        tileIdx = 0;
    elseif netIdx == 8 || netIdx == 11
        exportgraphics(fig1,fullfile(saveDir,"exploratory_figures", ...
            simName + ".pdf"),"ContentType",'vector','Append',true);
        clf(fig1)
        set(fig1,'units','inches','Position',[0,0,8.5,11])
        fig1 = tiledlayout(2,2);
        tileIdx = 0;
    end
    tileIdx = tileIdx + 1;
end
close(1)

%% Modulus vs. astral number, family of curves (krot)

fig2 = figure(2); clf;
set(fig2,'units','inches','Position',[0,0,6,6])
cmap = colormap(parula(nKrot));
hold on
for jdx = 1:nKrot
    errorbar(nFilPerAsterList,meanMod(:,jdx), ...
        (z/sqrt(numRep))*stdevMod(:,jdx),'Color',cmap(jdx,:), ...
        'LineWidth',0.75)
end
lg = legend(string(krot_list),'Location','northeast');
title(lg,{"Stiffness at centers"; "[pN/{\mu}m]"})
% cbar = colorbar('eastoutside','TickLabels',string(krot_list));
% cbar.Label.String = "Angular stiffness [pN*{\mu}m/rad]";
xlabel('Astral number a_n')
xticks(nFilPerAsterList)
ylabel('Shear modulus [pN/{\mu}m]')
yl = ylim;
ylim([0,yl(2)])
exportgraphics(fig2,fullfile(saveDir,"exploratory_figures", ...
    simName + "_krotfam.pdf"),'ContentType','vector')
