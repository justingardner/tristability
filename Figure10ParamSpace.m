function [rivScores,hystSwitches] = Figure10ParamSpace(simulations,fattcs,noises,condNames)
%Plots the parameters for which the model works like we want...
simIDs = 300:358;%500:548;%590;
%noises = [0.0 0.001 0.01 0.1 1.0 10.0];
%fattcs = [0.0 .01 .02 .03 0.04 .06 0.08 .12 0.16 .24 0.32 .48 0.64];
%simIDs = reshape(simIDs,7,7)';%Make into the right shape for the params we tested


nNoises = size(noises,2);
nFattcs = size(fattcs,2);
nConds = 10;
p.dT = .5;

mSize = 7;
LTM = .2;
LT = .2;
fontSize = 10;
titleFontSize = 12;

rivScores = cell(6,4);%Six conditions where we look for rivalry, four scores about it
hystSwitches = zeros(nNoises,nFattcs,4);%For conditions with hysteresis - fast and slow first switch, fast and slow cumulative

for noise=1:nNoises
    for fattc = 1:nFattcs
        for cond = 1:5
            binSum = simulations{noise,fattc,cond};
            leftRight = binSum([1 3],:)';
            leftRightCorr = corr(leftRight);
            leftRightCorr = leftRightCorr(2,1);

            [~,dominanceRecord] = max(binSum,[],1);
            changeLocs = [find(diff([-1 dominanceRecord]) ~= 0)]; % where does dominant unit change
            changeLocs(end+1) = size(dominanceRecord,2); % Ensures last button is included
            durs = diff(changeLocs); % Find dominance durations
            if(size(durs,2)>2)%Remove the first and last durations for accuracy
                durs = durs(2:end-1);
            end

            %[mu,sigma,logLikelihood] = bootstrapLogNormalFit(durs,'',1,false);
            mu=1;
            sigma=1;
            logLikelihood=1;

            rivScores{cond,1}(noise,fattc) = leftRightCorr;
            rivScores{cond,2}(noise,fattc) = mu;
            rivScores{cond,3}(noise,fattc) = sigma;
            rivScores{cond,4}(noise,fattc) = logLikelihood;
        end
    end
end

for simInd=1:size(simIDs,2)

    simID = simIDs(simInd);
    simFileName = strcat(num2str(simID),'_hysteresis.mat');
    load(simFileName);
    noiseInd = find(currNoise==noises);
    fattcInd = find(currFattc==fattcs);
    [OutThreshs,InThreshs,outConfs,inConfs,slowHystConfs,fastHystConfs] = Figure3_FirstSwitchThresholds(0,0,0,simID);
    hystSwitches(noiseInd,fattcInd,1) = fastHystConfs(2);
    hystSwitches(noiseInd,fattcInd,2) = slowHystConfs(2);

    [slowHystConfs fastHystConfs] = Figure4_ProportionedCumulativePercepts(0,0,simID);
    hystSwitches(noiseInd,fattcInd,3) = fastHystConfs(2);
    hystSwitches(noiseInd,fattcInd,4) = slowHystConfs(2);
end


valueNames = {'Mean dominance durations','Mean depth of suppression','Proportion of time with a dominant (>1.5x) node','Smallest dominance proportion of any node','Left/Right corr','Left/all corr','Center/all corr','Right/all corr','mu','sigma','logLikelihood'};
condNames = {'Rivalrous attended','Rivalrous unattended','Plaid attended','Plaid unattended','Swaps with blanks','Tristable','Fast rot out','Fast rot in','Slow rot out','Slow rot in'};

fastHysts = hystSwitches(:,:,1)-hystSwitches(:,:,2);
slowHysts = hystSwitches(:,:,3)-hystSwitches(:,:,4);

figure();

subplot(2,3,1);
imagesc(rivScores{1,1});
colorbar;
title('Rivalrous attended L/R correlation');

subplot(2,3,2);
imagesc(rivScores{2,1});
colorbar;
title('Rivalrous unattended L/R correlation');

subplot(2,3,3);
imagesc(rivScores{3,1});
colorbar;
title('Plaid attended L/R correlation');

subplot(2,3,4);
imagesc(rivScores{4,1});
colorbar;
title('Plaid unattended L/R correlation');

subplot(2,3,5);
imagesc(hystSwitches(:,:,1));
colorbar;
title('First switch fast hysts');

subplot(2,3,6);
imagesc(hystSwitches(:,:,2));
colorbar;
title('First switch slow hysts');


figure();

subplot(2,3,1);
imagesc(rivScores{1,1}<0);
colorbar;
title('Rivalrous attended - NEGATIVE correlation');

subplot(2,3,2);
imagesc(rivScores{2,1}>0);
colorbar;
title('Rivalrous unattended - POSITIVE correlation');

subplot(2,3,3);
imagesc(rivScores{3,1}>0);
colorbar;
title('Plaid attended - POSITIVE correlation');

subplot(2,3,4);
imagesc(rivScores{4,1}>0);
colorbar;
title('Plaid unattended - POSITIVE correlation');

subplot(2,3,5);
imagesc(fastHysts>0);
colorbar;
title('Fast hysts - Hysteresis');

subplot(2,3,6);
imagesc(slowHysts<0);
colorbar;
title('Slow hysts - Reverse Hysteresis');


figure();
subplot(2,1,1);
imagesc((rivScores{1,1}<0) + (rivScores{2,1}>0) + (rivScores{3,1}>0) + (rivScores{4,1}>0) + (hystSwitches(:,:,1)>0) + (hystSwitches(:,:,2)<0));
colorbar;
title('Regions?');
xlabel('Attentional self-facilitation');
ylabel('Noise');
drawPublishAxis('xTickOffset',-.005,'yTickOffset',.5,'whichAxis','both','lineWidth',1,'labelFontSize',fontSize,'titleStr','Model behavior across noise and facilitation parameter values','xAxisMajorTickLen',-1/64,'xAxisMinorTickLen',-1/128,'xAxisOffset',-.025,'xAxisMin',1,'xAxisMax',size(fattcs,2),'xTick',1:size(fattcs,2),'xTickLabel',fattcs,'xTickLabelSigfigs',2,'yTickLabelSigfigs',2,'yAxisMin',1,'yAxisMax',size(noises,2),'yAxisMajorTickLen',-1/64,'yTick',1:size(noises,2),'yTickLabel',noises)


subplot(2,1,2);
imagesc(rivScores{1,1}<0 & rivScores{2,1}>0 & rivScores{3,1}>0 & rivScores{4,1}>0 & hystSwitches(:,:,1)>0 & hystSwitches(:,:,2)<0);
colorbar;
title('Good zone');

