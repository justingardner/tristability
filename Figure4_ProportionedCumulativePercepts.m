function [slowHystConfs fastHystConfs] = Figure4_ProportionedCumulativePercepts(isSimulated,isFlatModel,subjID)
%This function plots a cumulative view of the total types of each press at
%each time point

SIcolor = [50 75 106]./255;
SOcolor = [145 60 85]./255;
FIcolor = [88 118 154]./255;
FOcolor = [210 112 141]./255;
makeFigure = true;

if(nargin>2)
    subjIDs=subjID;
elseif(isSimulated)
    if(isFlatModel)
        %simulatedTristableFlatData = [1054:1063];%original run
        %simulatedTristableFlatData = [1030:1040];%New run with .22 fattc
        simulatedTristableFlatData = 523;%testing
        subjIDs = simulatedTristableFlatData;
    else
        simulatedTristableModAttnData = [900:909];%[1013:1023];
        subjIDs = simulatedTristableModAttnData;
    end
else
    %goodIQRSubjs = [1821 1822 1823 1825 1826 1827 1828 1830 1833 1835 1836];%Feb 23 good IQRs
    goodIQRSubjs = [1821 1822 1823 1825 1826 1827 1828 1830 1833 1835];%March 9 good IQRs without corrupted 1836
    subjIDs = goodIQRSubjs;
end

nSubj = size(subjIDs,2);
subjResults = cell(nSubj,1);

mSize = 7;
LTM = .2;
LT = .2;
fontSize = 10;
titleFontSize = 12;
nBins = 10;

for s=1:nSubj
    subjID = subjIDs(s);
    if(subjID>1100)
        load(strcat(num2str(subjID),'_hysteresis1.mat'));
        subjResults{s} = results;
        load(strcat(num2str(subjID),'_hysteresis2.mat'));
        subjResults{s} = [subjResults{s} ;results];
    else
        load(strcat(num2str(subjID),'_hysteresis.mat'));
        subjResults{s} = results;
    end
end

speeds = cell(nSubj,1);
rotDirs = cell(nSubj,1);
resps = cell(nSubj,1);

%Get responses
for s = 1:nSubj
    results = subjResults{s};
    subjConds = size(results,1);
    resps{s} = cell(subjConds,1);
    for cond = 1:subjConds
        resp = results{cond,2};
        %Get condition parameters
        speeds{s}(cond) = results{cond,1}{1};
        rotDirs{s}(cond) = results{cond,1}{2};
        %Flip responses so they're all aligned. When rotDir is -1, the
        %orientation index started at the end and moved backwards instead of
        %the beginning and going forward through orientations = [5:.5:25];
        %By flipping, they will all go from 5-25 by .5
        if(rotDirs{s}(cond) == -1)
            resp = fliplr(resp);
        end
        resps{s}{cond} = resp;
    end
end

slowCumulativeRespsOutBySubj = [];
fastCumulativeRespsOutBySubj = [];
slowCumulativeRespsInBySubj = [];
fastCumulativeRespsInBySubj = [];

for s=1:nSubj
    %Start a new collection for each subject
    slowCumulativeRespsOut = zeros(5,nBins);%Null, left, right, up, riv presses for each bin
    fastCumulativeRespsOut = zeros(5,nBins);%Null, left, right, up, riv presses for each bin

    slowCumulativeRespsIn = zeros(5,nBins);%Null, left, right, up, riv presses for each bin
    fastCumulativeRespsIn = zeros(5,nBins);%Null, left, right, up, riv presses for each bin

    subjConds = size(results,1);
    for cond=1:subjConds
        resp = resps{s}{cond};
        orientationBinSize = size(resp,2)/nBins;
        cumulativeResps = zeros(5,nBins);
        for b=1:size(resp,2)
            oriBin = ceil(b/orientationBinSize);
            button = resp(b);
            if(button==0)
                cumulativeResps(1,oriBin) = cumulativeResps(1,oriBin)+1;
            elseif(button==124)
                cumulativeResps(2,oriBin) = cumulativeResps(2,oriBin)+1;
            elseif(button==125)
                cumulativeResps(3,oriBin) = cumulativeResps(3,oriBin)+1;
            elseif(button==127)
                cumulativeResps(4,oriBin) = cumulativeResps(4,oriBin)+1;
            end
            if(button==124 || button==125)
                cumulativeResps(5,oriBin) = cumulativeResps(5,oriBin)+1;
            end
        end
        if(rotDirs{s}(cond)==1)%Rotating out...
            if(speeds{s}(cond)==4)%4 seconds per degree (slow) condition
                slowCumulativeRespsOut = slowCumulativeRespsOut+cumulativeResps;
            elseif(speeds{s}(cond)==.5)%.5 seconds per degree (fast) condition
                fastCumulativeRespsOut = fastCumulativeRespsOut+cumulativeResps;
            end
        elseif(rotDirs{s}(cond)==-1)%Rotating in...
            if(speeds{s}(cond)==4)%4 seconds per degree (slow) condition
                slowCumulativeRespsIn = slowCumulativeRespsIn+cumulativeResps;
            elseif(speeds{s}(cond)==.5)%.5 seconds per degree (fast) condition
                fastCumulativeRespsIn = fastCumulativeRespsIn+cumulativeResps;
            end
        end
    end
    slowCumulativeRespsOutBySubj = cat(3,slowCumulativeRespsOutBySubj,slowCumulativeRespsOut);
    fastCumulativeRespsOutBySubj = cat(3,fastCumulativeRespsOutBySubj,fastCumulativeRespsOut);
    slowCumulativeRespsInBySubj = cat(3,slowCumulativeRespsInBySubj,slowCumulativeRespsIn);
    fastCumulativeRespsInBySubj = cat(3,fastCumulativeRespsInBySubj,fastCumulativeRespsIn);
end

%Get total records across subjects
slowCumulativeRespsOut = sum(slowCumulativeRespsOutBySubj,3);
fastCumulativeRespsOut = sum(fastCumulativeRespsOutBySubj,3);
slowCumulativeRespsIn = sum(slowCumulativeRespsInBySubj,3);
fastCumulativeRespsIn = sum(fastCumulativeRespsInBySubj,3);

%Get proportion of rivalry out of the three possible reported percept types
slowCumulativeProportionOut = slowCumulativeRespsOut(5,:)./(slowCumulativeRespsOut(4,:)+slowCumulativeRespsOut(5,:));
fastCumulativeProportionOut = fastCumulativeRespsOut(5,:)./(fastCumulativeRespsOut(4,:)+fastCumulativeRespsOut(5,:));
slowCumulativeProportionIn = slowCumulativeRespsIn(5,:)./(slowCumulativeRespsIn(4,:)+slowCumulativeRespsIn(5,:));
fastCumulativeProportionIn = fastCumulativeRespsIn(5,:)./(fastCumulativeRespsIn(4,:)+fastCumulativeRespsIn(5,:));

%Bootstrap to find error bars
nBootstraps = 200;
if(nSubj==1)
    nBootstraps=10;
end
SOBoots = zeros(nBootstraps,nBins);
SIBoots = zeros(nBootstraps,nBins);
FOBoots = zeros(nBootstraps,nBins);
FIBoots = zeros(nBootstraps,nBins);

SOFits = zeros(nBootstraps,3);
SIFits = zeros(nBootstraps,3);
FOFits = zeros(nBootstraps,3);
FIFits = zeros(nBootstraps,3);

slowShifts = zeros(nBootstraps,1);
fastShifts = zeros(nBootstraps,1);

for b=1:nBootstraps
    currSOBootstrap = zeros(nSubj,nBins);
    currSIBootstrap = zeros(nSubj,nBins);
    currFOBootstrap = zeros(nSubj,nBins);
    currFIBootstrap = zeros(nSubj,nBins);
    for s=1:nSubj
        randSubjInd = min(ceil(rand()*nSubj),nSubj);
        randSubjRiv = slowCumulativeRespsOutBySubj(5,:,randSubjInd);
        randSubjFus = slowCumulativeRespsOutBySubj(4,:,randSubjInd);
        currSOBootstrap(s,:) = randSubjRiv./(randSubjRiv+randSubjFus);
        
        randSubjInd = min(ceil(rand()*nSubj),nSubj);
        randSubjRiv = slowCumulativeRespsInBySubj(5,:,randSubjInd);
        randSubjFus = slowCumulativeRespsInBySubj(4,:,randSubjInd);
        currSIBootstrap(s,:) = randSubjRiv./(randSubjRiv+randSubjFus);
        
        randSubjInd = min(ceil(rand()*nSubj),nSubj);
        randSubjRiv = fastCumulativeRespsOutBySubj(5,:,randSubjInd);
        randSubjFus = fastCumulativeRespsOutBySubj(4,:,randSubjInd);
        currFOBootstrap(s,:) = randSubjRiv./(randSubjRiv+randSubjFus);
        
        randSubjInd = min(ceil(rand()*nSubj),nSubj);
        randSubjRiv = fastCumulativeRespsInBySubj(5,:,randSubjInd);
        randSubjFus = fastCumulativeRespsInBySubj(4,:,randSubjInd);
        currFIBootstrap(s,:) = randSubjRiv./(randSubjRiv+randSubjFus);
    end
    
    SOBoots(b,:) = mean(currSOBootstrap,1);
    SIBoots(b,:) = mean(currSIBootstrap,1);
    FOBoots(b,:) = mean(currFOBootstrap,1);
    FIBoots(b,:) = mean(currFIBootstrap,1);
    
    [SOFits(b,1) SOFits(b,2) SOFits(b,3)] = getSigmoidFit(SOBoots(b,:),0);
    [SIFits(b,1) SIFits(b,2) SIFits(b,3)] = getSigmoidFit(SIBoots(b,:),0);
    [FOFits(b,1) FOFits(b,2) FOFits(b,3)] = getSigmoidFit(FOBoots(b,:),0);
    [FIFits(b,1) FIFits(b,2) FIFits(b,3)] = getSigmoidFit(FIBoots(b,:),0);
    
    SOCrosses = sum(SOBoots(b,:)>.5)./size(SOBoots,2);
    SOCrosses = (SOCrosses==1 | SOCrosses==0);
    SICrosses = sum(SIBoots(b,:)>.5)./size(SIBoots,2);
    SICrosses = (SICrosses==1 | SICrosses==0);
    FOCrosses = sum(FOBoots(b,:)>.5)./size(FOBoots,2);
    FOCrosses = (FOCrosses==1 | FOCrosses==0);
    FICrosses = sum(FIBoots(b,:)>.5)./size(FIBoots,2);
    FICrosses = (FICrosses==1 | FICrosses==0);
    
    if(SOCrosses)
        SOFits(b,:)=NaN;
    end
    if(SICrosses)
        SIFits(b,:)=NaN;
    end
    if(FOCrosses)
        FOFits(b,:)=NaN;
    end
    if(FICrosses)
        FIFits(b,:)=NaN;
    end
        
    
    slowShifts(b,1) = SOFits(b,2)-SIFits(b,2);
    fastShifts(b,1) = FOFits(b,2)-FIFits(b,2);
end

%Get 95% confidence intervals from bootstraps
SOconfs = zeros(2,nBins);%Lower, upper then degrees disp
SIconfs = zeros(2,nBins);
FOconfs = zeros(2,nBins);
FIconfs = zeros(2,nBins);

twoPointFive = ceil(nBootstraps*.025);
ninetySevenPointFive = floor(nBootstraps*.975);

slowHystConfs = zeros(1,3);%low conf, mean, high conf with 95%CI
fastHystConfs = zeros(1,3);%low conf, mean, high conf with 95%CI

sortedSlowShifts = sort(slowShifts);
sortedFastShifts = sort(fastShifts);

slowHystConfs(1) = sortedSlowShifts(twoPointFive);
slowHystConfs(2) = mean(sortedSlowShifts);
slowHystConfs(3) = sortedSlowShifts(ninetySevenPointFive);

fastHystConfs(1) = sortedFastShifts(twoPointFive);
fastHystConfs(2) = mean(sortedFastShifts);
fastHystConfs(3) = sortedFastShifts(ninetySevenPointFive);

for d=1:nBins
    bootStraps = SOBoots(:,d);
    sortedBoots = sort(bootStraps);
    SOconfs(1,d) = slowCumulativeProportionOut(d)-sortedBoots(twoPointFive);
    SOconfs(2,d) = sortedBoots(ninetySevenPointFive)-slowCumulativeProportionOut(d);

    bootStraps = SIBoots(:,d);
    sortedBoots = sort(bootStraps);
    SIconfs(1,d) = slowCumulativeProportionIn(d)-sortedBoots(twoPointFive);
    SIconfs(2,d) = sortedBoots(ninetySevenPointFive)-slowCumulativeProportionIn(d);

    bootStraps = FOBoots(:,d);
    sortedBoots = sort(bootStraps);
    FOconfs(1,d) = fastCumulativeProportionOut(d)-sortedBoots(twoPointFive);
    FOconfs(2,d) = sortedBoots(ninetySevenPointFive)-fastCumulativeProportionOut(d);

    bootStraps = FIBoots(:,d);
    sortedBoots = sort(bootStraps);
    FIconfs(1,d) = fastCumulativeProportionIn(d)-sortedBoots(twoPointFive);
    FIconfs(2,d) = sortedBoots(ninetySevenPointFive)-fastCumulativeProportionIn(d);
end

binDisps = 1+[0:nBins].*(39./nBins);
binDisps = mean([binDisps(1:end-1);binDisps(2:end)]);

if(makeFigure)

    %Plot
    h = figure('rend','painter','pos',[10 10 900 1650]);

    for subplotID = 1:2
        if(subplotID==1)
            rivalryProportionOut = slowCumulativeProportionOut';
            rivalryProportionIn = slowCumulativeProportionIn';
            outConfs = SOconfs;
            inConfs = SIconfs;
            outCol = SOcolor;
            inCol = SIcolor;
            title('Slow condition');
        elseif(subplotID==2)
            rivalryProportionOut = fastCumulativeProportionOut';
            rivalryProportionIn = fastCumulativeProportionIn';
            outConfs = FOconfs;
            inConfs = FIconfs;
            outCol = FOcolor;
            inCol = FIcolor;
            title('Fast condition');
        end

        subplot(1,2,subplotID);

        plot(binDisps,rivalryProportionOut,'-','color',outCol,'LineWidth', LT,'MarkerEdgeColor','w');
        hold on;
        plot(binDisps,rivalryProportionIn,'-','color',inCol,'LineWidth', LT,'MarkerEdgeColor','w');


        legendText = {' Rotating out',' Rotating in'};
        legendColors = {outCol,inCol};
        legendX = 30;
        legendY = 1.1;
        legendSpacing = .045;
        legendGap = 2;
        legendLineSegSize = 10;
        for i=1:2
            plot([legendX-legendLineSegSize-legendGap legendX-legendGap],[legendY-(legendSpacing*i) legendY-(legendSpacing*i)],'lineWidth',LT,'color',legendColors{i});
            text(legendX,legendY-(legendSpacing*i),legendText{i},'VerticalAlignment','middle','HorizontalAlignment','left','FontAngle','oblique','FontSize',fontSize);
        end

        myerrorbar(binDisps+.1,rivalryProportionOut,'yError',outConfs(1,:),'yErrorBarType','lower','Symbol','.','Color',outCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);
        myerrorbar(binDisps-.1,rivalryProportionIn,'yError',inConfs(1,:),'yErrorBarType','lower','Symbol','.','Color',inCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);

        myerrorbar(binDisps+.1,rivalryProportionOut,'yError',outConfs(2,:),'yErrorBarType','upper','Symbol','.','Color',outCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);
        myerrorbar(binDisps-.1,rivalryProportionIn,'yError',inConfs(2,:),'yErrorBarType','upper','Symbol','.','Color',inCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);

        myerrorbar(binDisps,rivalryProportionOut,'Symbol','o','Color',outCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);
        myerrorbar(binDisps,rivalryProportionIn,'Symbol','o','Color',inCol,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM,'yTeelen',0);

        xaxis(0,40);
        xlabel('Angle between Gabors (degrees)');
        ylabel('Proportion of rivalrous reports');
        drawPublishAxis('xTickOffset',-.005,'yTickOffset',.5,'whichAxis','both','lineWidth',1,'labelFontSize',fontSize,'titleStr','','xAxisMajorTickLen',-1/64,'xAxisMinorTickLen',-1/128,'xAxisOffset',-.025,'xAxisMin',0,'xAxisMax',40,'xTick',0.0:10.0:40.0,'xTickLabel',10:10:50,'xTickLabelSigfigs',0,'yTickLabelSigfigs',2,'yAxisMin',0,'yAxisMax',1,'yAxisMajorTickLen',-1/64,'yTick',0.0:.25:1.0,'yTickLabel',0.0:.25:1.0)

    end


    set(h,'PaperUnits','inches');
    x_width = 10;
    y_width = 10;
    set(h,'PaperPosition',[.5 .5 x_width y_width]);
    %saveas(h,'Figure4_CumulativePercepts_raw.pdf')
    if(isSimulated)
        if(isFlatModel)
            saveas(h,'Figure5b_CumulativePercepts_Flatb.pdf');
        else
            saveas(h,'Figure7b_CumulativePercepts_ModAttnb.pdf');
        end
    else
        saveas(h,'Figure4b_CumulativePercepts_Realb.pdf');
    end
end
end

