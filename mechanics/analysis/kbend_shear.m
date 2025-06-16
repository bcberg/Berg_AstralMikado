% kbend_shear.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "kbend_shear250228";
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';

newReports = false;
if newReports
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/kbend_shear/";
    
    % ensure the following parameters match those used in
    % kbend_shear.cym.tpl and its driver (for numRep)
    dens = 7.5;  % N * len_fil / A
    len_fil = 1;
    D = 10; % side length of square region
    nFil = dens * D^2 / len_fil;
    nFilPerAsterList = ([1:6,8:4:24])'; % "astral number"
    nNetTypes = length(nFilPerAsterList);
    dataStartLine = 6;  % line of first numerical data in report files
    ptsPerFiber = 6;
    forceVals = (5*(0:3)/3)';
    nForces = length(forceVals);
    numRep = 10;
    kbend_list = [1,5*(1:9)];   % [pN*um^2] (units according to Cytosim docs)
    nKbend = length(kbend_list);

    runsPerKbend = nNetTypes * nForces;
    runsPerRep = nKbend * nNetTypes * nForces;
    shearModuli = zeros(nNetTypes,nKbend,numRep);
    for repIdx = 1:numRep
        for kIdx = 1:nKbend
            for netIdx = 1:nNetTypes
                nFilPerAster = nFilPerAsterList(netIdx);
                netLabel = sprintf('an%02i',nFilPerAster);
                if repIdx == 1
                    initCoMs.(netLabel) = zeros(nForces,2,nKbend,numRep);
                    finalCoMs.(netLabel) = zeros(nForces,2,nKbend,numRep);
                    disps.(netLabel) = zeros(nForces,2,nKbend,numRep);
                    lineFits.(netLabel) = cell(nKbend,numRep);
                    fitStats.(netLabel) = cell(nKbend,numRep);
                end
                cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
                    'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);
                unshiftedRunVals = (netIdx-1)*nForces:(netIdx*nForces - 1);
                kbendShift = (kIdx - 1) * runsPerKbend;
                repShift = (repIdx - 1) * runsPerRep;
                trueRunVals = unshiftedRunVals + kbendShift + repShift;
                [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,trueRunVals, ...
                    cytoparams);
                initCoMs.(netLabel)(:,:,kIdx,repIdx) = thisInitCoM;
                finalCoMs.(netLabel)(:,:,kIdx,repIdx) = thisFinalCoM;
                disps.(netLabel)(:,:,kIdx,repIdx) = thisFinalCoM ...
                    - thisInitCoM;
                [thisModulus,linefit,gof] = getModulus(thisInitCoM, ...
                    thisFinalCoM,forceVals,"shear");
                shearModuli(netIdx,kIdx,repIdx) = thisModulus;
                linefits.(netLabel){kIdx,repIdx} = linefit;
                fitStats.(netLabel){kIdx,repIdx} = gof;
            end
        end
    end
    save(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
        'D','nFil','nFilPerAsterList','forceVals','numRep','nKbend', ...
        'nNetTypes','kbend_list','initCoMs','finalCoMs','disps', ...
        'linefits','fitStats','shearModuli')
    clearvars -except simName saveDir
end
load(fullfile(saveDir,"mat_files",simName + ".mat"),'dens','len_fil', ...
        'D','nFil','nFilPerAsterList','forceVals','numRep','nKbend', ...
        'nNetTypes','kbend_list','initCoMs','finalCoMs','disps', ...
        'linefits','fitStats','shearModuli')

lgray = [197,190,181]/255;

[stdevMod,meanMod] = std(shearModuli,0,3);
z = norminv(0.975);

%% Modulus vs. kbend, separate axes for each astral number

fig1 = figure(1);
set(fig1,'units','inches','Position',[0,0,8.5,11])
fig1 = tiledlayout(2,2);
tileIdx = 1;
for netIdx = 1:nNetTypes
    nexttile(tileIdx)
    hold on
    for kIdx = 1:nKbend
        plot(kbend_list(kIdx),squeeze(shearModuli(netIdx,kIdx,:)), ...
            'Marker','o','MarkerSize',3,'MarkerFaceColor',lgray, ...
            'MarkerEdgeColor',lgray)
    end
    errorbar(kbend_list,meanMod(netIdx,:), ...
        (z/sqrt(numRep))*stdevMod(netIdx,:),'-ob','LineWidth',1.5)
    hold off
    xlabel('Filament bending rigidity [pN*{\mu}m^2]')
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

%% Modulus vs. astral number, family of curves (kbend)

fig2 = figure(2);
set(fig2,'units','inches','Position',[0,0,6,6])
cmap = colormap(parula(nKbend));
hold on
for jdx = 1:nKbend
    errorbar(nFilPerAsterList,meanMod(:,jdx), ...
        (z/sqrt(numRep))*stdevMod(:,jdx),'Color',cmap(jdx,:), ...
        'LineWidth',0.75)
end
lg = legend(string(kbend_list),'Location','northeast');
title(lg,{"Bending rigidity"; "[pN*{\mu}m^2]"})
xlabel('Astral number a_n')
xticks(nFilPerAsterList)
ylabel('Shear modulus [pN/{\mu}m]')
yl = ylim;
ylim([0,yl(2)])
exportgraphics(fig2,fullfile(saveDir,"exploratory_figures", ...
    simName + "_kbendfam.pdf"),'ContentType','vector')
