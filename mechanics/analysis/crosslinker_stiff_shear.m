% crosslinker_stiff_shear.m
% Brady Berg
clear
close all
format short
format compact

%% Parse report files

simName = "crosslinker_stiff_shear251127";
% Ubuntu filepath
saveDir = '~/Documents/astral-mikado-data';
% Windows filepath
% saveDir = 'C:\Users\bcber\Documents\astral-mikado-data';
dataSubdir = 'mat_files';
figSubdir = 'exploratory_figures';

newReports = true;
if newReports
    filtering.TF = true;
    dir = "/home/bcberg/Documents/AstralMikadoCYM/runs/crosslinker_stiff_shear/";
    
    % ensure the following parameters match those used in
    % crosslinker_dens_shear.cym.tpl and its driver (for numRep)
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
    numRep = 15;
    crosslinker_stiff_list = 12.5 * 2.^(0:5);  % [pN/um]
    nCrslnkStiff = length(crosslinker_stiff_list);
    
    runsPerNetType = nCrslnkStiff * nForces;
    runsPerRep = nNetTypes * nCrslnkStiff * nForces;
    shearModuli = zeros(nNetTypes,nCrslnkStiff,numRep);
    percStatus = false(nNetTypes,nCrslnkStiff,numRep,2);
    for repIdx = 1:numRep
        for lnkIdx = 1:nCrslnkStiff
            for netIdx = 1:nNetTypes
                nFilPerAster = nFilPerAsterList(netIdx);
                netLabel = sprintf('an%02i',nFilPerAster);
                if repIdx == 1
                    initCoMs.(netLabel) = zeros(nForces,2,nCrslnkStiff, ...
                        numRep);
                    finalCoMs.(netLabel) = zeros(nForces,2,nCrslnkStiff, ...
                        numRep);
                    disps.(netLabel) = zeros(nForces,2,nCrslnkStiff, ...
                        numRep);
                    lineFits.(netLabel) = cell(nCrslnkStiff,numRep);
                    fitStats.(netLabel) = cell(nCrslnkStiff,numRep);
                end
                cytoparams = struct('nFil',nFil,'nFilPerAster',nFilPerAster,...
                    'dataStartLine',dataStartLine,'ptsPerFiber',ptsPerFiber);
                unshiftedRunVals = (lnkIdx - 1)*nForces:(lnkIdx*nForces - 1);
                netShift = (netIdx - 1) * runsPerNetType;
                repShift = (repIdx - 1) * runsPerRep;
                trueRunVals = unshiftedRunVals + netShift + repShift;
                % check network percolation status
                [percTF,~,~] = checkLinks(dir,trueRunVals(1),cytoparams);
                filtering.spanCheck = percTF(1);
                percStatus(netIdx,lnkIdx,repIdx,:) = reshape(percTF,[1,1,1,2]);
                % analyze network position data
                [thisInitCoM,thisFinalCoM] = parseForceSweep(dir,trueRunVals, ...
                    cytoparams);
                initCoMs.(netLabel)(:,:,lnkIdx,repIdx) = thisInitCoM;
                finalCoMs.(netLabel)(:,:,lnkIdx,repIdx) = thisFinalCoM;
                disps.(netLabel)(:,:,lnkIdx,repIdx) = thisFinalCoM ...
                    - thisInitCoM;
                % estimate a modulus
                [thisModulus,linefit,gof] = getModulus(thisInitCoM, ...
                    thisFinalCoM,forceVals,"shear",filtering);
                shearModuli(netIdx,lnkIdx,repIdx) = thisModulus;
                linefits.(netLabel){lnkIdx,repIdx} = linefit;
                fitStats.(netLabel){lnkIdx,repIdx} = gof;
            end
        end
    end
    save(fullfile(saveDir,dataSubdir,simName + ".mat"),'dens','len_fil', ...
        'D','nFil','nFilPerAsterList','forceVals','numRep','nNetTypes', ...
        'crosslinker_stiff_list','nCrslnkStiff','initCoMs','finalCoMs', ...
        'disps','linefits','fitStats','shearModuli','percStatus', ...
        'filtering')
    clearvars -except simName saveDir dataSubdir figSubdir
end

load(fullfile(saveDir,dataSubdir,simName + ".mat"),'dens','len_fil', ...
    'D','nFil','nFilPerAsterList','forceVals','numRep','nNetTypes', ...
    'crosslinker_stiff_list','nCrslnkStiff','initCoMs','finalCoMs', ...
    'disps','linefits','fitStats','shearModuli','percStatus', ...
    'filtering')

lgray = [197,190,181]/255;

%% Modulus vs. crosslinker stiffness, separate axes for each astral number

[stdevMod,meanMod] = std(shearModuli,0,3);
z = norminv(0.975);

fig1 = figure(1);
set(fig1,'units','inches','Position',[0,0,8.5,11])
fig1 = tiledlayout(2,2);
tileIdx = 1;
for netIdx = 1:nNetTypes
    nexttile(tileIdx)
    hold on
    for lnkIdx = 1:nCrslnkStiff
        plot(crosslinker_stiff_list(lnkIdx),squeeze(shearModuli(netIdx,lnkIdx,:)), ...
            'Marker','o','MarkerSize',3,'MarkerFaceColor',lgray, ...
            'MarkerEdgeColor',lgray)
    end
    errorbar(crosslinker_stiff_list,meanMod(netIdx,:), ...
        (z/sqrt(numRep))*stdevMod(netIdx,:),'-ob','LineWidth',1.5)
    hold off
    xlabel('Crosslinker stiffness [pN/{\mu}m')
    ylabel('Shear modulus [pN/{\mu}m]')
    title(sprintf('a_n = %i',nFilPerAsterList(netIdx)))
    xscale('log')
    
    if netIdx == 4
        exportgraphics(fig1,fullfile(saveDir,figSubdir, ...
            simName + ".pdf"),"ContentType",'vector');
        clf(fig1)
        set(fig1,'units','inches','Position',[0,0,8.5,11])
        fig1 = tiledlayout(2,2);
        tileIdx = 0;
    elseif netIdx == 8 || netIdx == 11
        exportgraphics(fig1,fullfile(saveDir,figSubdir, ...
            simName + ".pdf"),"ContentType",'vector','Append',true);
        clf(fig1)
        set(fig1,'units','inches','Position',[0,0,8.5,11])
        fig1 = tiledlayout(2,2);
        tileIdx = 0;
    end
    tileIdx = tileIdx + 1;
end
    
