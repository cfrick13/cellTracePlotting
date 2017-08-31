function plot_tracked_cells_hypothesis(FileName,timeCutoff,sortTime,cmapchoice,groupBasedOn,colorBasedOn,subplotGroupBasedOn,labelBasedOn,basalLength,smadnormalizestr,reporternormalizestr,includeStr,excludeStr,SmadTotalOrMedian,ReporterTotalOrMedian,celltracestr)
% close all

%important parameters
%   medfiltnum
%   jittercut
medfiltnum = 3;

%determine the location of the matlab function and establish export
%directory in relation to that filepath
mdir = mfilename('fullpath');
[~,b] = regexp(mdir,'Tracking\w*/');
if isempty(b)
    [~,b] = regexp(mdir,'Tracking\w*\');
end
parentdir = mdir(1:b); %specifies folder in which all analysis is being done


addpath([parentdir 'Colormaps'])
loaddir = strcat(parentdir,'Export'); %specifies where data is exported
cd(loaddir);

[~,b] = regexp(mdir,'/');
if isempty(b)
    [~,b] = regexp(mdir,'\');
end

mfiledir =mdir(1:b(end)); %specifies location of matlab function file

%load the exported tracking structure
cd(loaddir)
load(FileName)

%load metadata associated with the experiment (requires manual input if there is ambiguity)
[a,~] = regexp(FileName,'_tracking');
fileDateStr = FileName(1:a-1);
datequery = strcat(fileDateStr,'*metaData.mat');
cd(loaddir)
filelist = dir(datequery);
if length({filelist.name}) ==1
    metaData = load(char(filelist.name));
else
    error('need to extract metadata')
    filename = uigetfile();
    metaData = load(filename);
end
%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
timeVec = metaData.timeVec;


%load information regarding doses and scenes and tgfbeta addition
datequery = strcat(fileDateStr,'*DoseAndScene*');
cd(loaddir)
filelist = dir(datequery);
if isempty(filelist)
    error('run makeDoseStruct')
    dosestruct = makeDoseStruct; %run function to make doseStruct
else
    dosestructstruct = load(char(filelist.name));
    dosestruct = dosestructstruct.dosestruct;
end


%determine the scenes present in the experiment in order to combine metadata
scenestr = 'scene';
sceneListArray = vertcat({exportStruct.(scenestr)});
sceneList = unique(sceneListArray);
sceneListArrayTwo = vertcat({dosestruct.(scenestr)});

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr = fnamesz{p};
        var = dosestruct(indicestwo).(fnstr);
        if isnumeric(var)
            [exportStruct(indices).(fnstr)] = deal(var);
        else
            [exportStruct(indices).(fnstr)] = deal(char(var));
        end
    end
end

%combine the exportStruct information with dosestruct information
for i=1:length(sceneList)
    sceneChoice=sceneList{i};
    indices = strcmp(sceneListArray,sceneChoice);
    indicestwo = strcmp(sceneListArrayTwo,sceneChoice);
    
    fnamesz = fieldnames(dosestruct);
    for p = 1:length(fnamesz)
        fnstr1 = 'dosestr';
        fnstr2 = 'conditions';
        fnstr3 = 'expdate';
        var1 = dosestruct(indicestwo).(fnstr1);
        var2 = dosestruct(indicestwo).(fnstr2);
        [exportStruct(indices).doseAndCondition] = deal([char(var1) ' ' char(var2)]);
        [exportStruct(indices).doseconddate] = deal([char(var1) ' ' fileDateStr ' ' char(var2)]);
        [exportStruct(indices).conddate] = deal([fileDateStr ' ' char(var2)]);
        [exportStruct(indices).(fnstr3)] = deal(fileDateStr);
    end
end


%determine details needed for plotting such as when Tgfbeta is added, etc
stimulationFrame = exportStruct(1).tgfFrame;
numberOfFrames = size(timeVec,2);
finalFrame = numberOfFrames;


%need to determine the number of scenes present and choose the time vector
%depending on the scene from which it was imaged
numberOfCells = length(indices);
timeMatrix = zeros(numberOfCells,finalFrame);

for i=1:numberOfCells
    sceneChoice=exportStruct(i).scene;
    idxtwo = strcmp(sceneListArrayTwo,sceneChoice);% sceneListArrayTwo = vertcat({dosestruct.(scenestr)});
    stimulationFrame = dosestruct(idxtwo).tgfFrame+1;
    timeMatrix(i,:) = timeVec(idxtwo,1:finalFrame)-timeVec(1,stimulationFrame); %subtract by time closest to Tgfbeta addition
    exportStruct(i).timeMatrix = timeMatrix(i,:);
    exportStruct(i).cellID = num2str(exportStruct(i).cellID);
    exportStruct(i).cellNum = num2str(i);
end





%% group based on desired criteria, set colormap, etc
%choose all cells that match specified chosen clone input strings
allList = {exportStruct.('doseconddate')};
[~,~,~,d] = regexp(allList,includeStr);
subidx = ~cellfun(@isempty,d,'UniformOutput',1);
exportStruct = exportStruct(subidx);

%exclude all cells that match specified exludeNames input strings
allList = {exportStruct.('doseconddate')};
[~,~,~,d] = regexp(allList,excludeStr);
subidx = ~cellfun(@isempty,d,'UniformOutput',1);
exportStruct = exportStruct(~subidx);




%% extract cell traces and plot
%function to exract the cell traces, normalized and not
indices = true(1,length(exportStruct)); %choose all cells initially
cTracks = struct();
cTracks.Smad = extractTraces(exportStruct,indices,'NucEGFP',finalFrame,stimulationFrame,basalLength,medfiltnum);
cTracks.Reporter = extractTraces(exportStruct,indices,'NucRFP',finalFrame,stimulationFrame,basalLength,medfiltnum);



if strcmp(SmadTotalOrMedian,'nuccytoXabundance')
    smadtraces1 = cTracks.('Smad').('nuccyto').('none');
    smadtraces2 = cTracks.('Smad').('median').('abundance');
    smadtraces = smadtraces1.*smadtraces2;
else
   smadtraces = cTracks.('Smad').(SmadTotalOrMedian).(smadnormalizestr); 
end

reportertraces = cTracks.('Reporter').(ReporterTotalOrMedian).(reporternormalizestr);
    

%% division identifying and interpolating function
cellMat = struct();
[tracemat,divtracemat,allmat,nanidx] = divisionFunction(cTracks,smadtraces,stimulationFrame,timeMatrix);
cellMat.Smad.division = divtracemat;
cellMat.Smad.traces = tracemat;
cellMat.Smad.all = allmat;
[tracemat,divtracemat,allmat,~] = divisionFunction(cTracks,reportertraces,stimulationFrame,timeMatrix);
cellMat.Reporter.division = divtracemat;
cellMat.Reporter.traces = tracemat;
cellMat.Reporter.all = allmat;

%% update the color values based on normalized or not
exportSubStruct = exportStruct(nanidx);
timeMatrix = timeMatrix(nanidx,:);

%% set cutoff based on time
tvec = timeMatrix(1,:);
timeidx = find(tvec>timeCutoff);
celltraces = cellMat.Smad.all;
celltime = celltraces(:,timeidx(1));
timenan = ~isnan(celltime);

exportStruct = exportSubStruct;
exportSubStruct = exportStruct(timenan);
newCellMat = updatealltraces(cellMat,timenan);
cellMat = newCellMat;
timeMatrix = timeMatrix(timenan,:);


%% sort traces
sortstr = 'Smad';
[sortedCellMat,sortidx] = sorttraces(cellMat,sortstr,timeMatrix,sortTime);
cellMat = sortedCellMat;

%% update the color values based on sorted or not
exportStruct = exportSubStruct;
timeMatrixOld = timeMatrix;
for i = 1:length(sortidx)
exportSubStruct(i) = exportStruct(sortidx(i));
timeMatrix(i,:) = timeMatrixOld(sortidx(i),:);
end

%% establish lists and listArrays for coloring and labeling
%establish the color map for plotting
%choose which field based upon which each cell trace will get colored
if strcmp(colorBasedOn,'individually')
    coloringChoice = 'scene';
else
    coloringChoice = colorBasedOn; %color traces based on what
end
[coloringList,coloringArrayList] = listSpitter(exportSubStruct,coloringChoice);
cmap = determineCMAP(cmapchoice,coloringList);

%assign a color array using the created colormap based on the choices above
colormapindividualMatrix = zeros(length(coloringArrayList),size(cmap,2));
for i=1:length(coloringArrayList)
    cA = coloringArrayList{i};
    idx = strcmp(coloringList,cA);
    colormapindividualMatrix(i,:) = cmap(idx,:); %this is used to color each individual trace according to the colorBasedOn criteria selected
end
randmat = normrnd(0,0.1,size(colormapindividualMatrix,1),3);
colormapindividualMatrix = colormapindividualMatrix + randmat;
colormapindividualMatrix(colormapindividualMatrix>1) = 1;
colormapindividualMatrix(colormapindividualMatrix<0) = 0;

%GROUP BASED ON
[groupList,groupArrayList] = listSpitter(exportSubStruct,groupBasedOn);

%assign a color for grouped median plots
groupColorMap = determineCMAP(cmapchoice,groupList);
groupColorMatrix = zeros(length(groupArrayList),size(groupColorMap,2));
for i=1:length(groupArrayList)
    cA = groupArrayList{i};
    idx = strcmp(groupList,cA);
    groupColorMatrix(i,:) = groupColorMap(idx,:); %this is used to color each individual trace according to the colorBasedOn criteria selected
end

%establish subplot groups
[subplotList,subplotArrayList] = listSpitter(exportSubStruct,subplotGroupBasedOn);

%establish the string array to use for labeling traces
[labelList,labelArrayList] = listSpitter(exportSubStruct,labelBasedOn);




%% how to make correlations

corrMat.Smad.abundance = cellMat.Smad.all;
corrMat.Reporter.abundance = cellMat.Reporter.all;
trueCorr = struct();
%1. abundance, foldchange, difference
    %2. rate, integral
fnames = fieldnames(corrMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    abundance = corrMat.(fnstr).abundance;
    
    valout = fcanddiff(abundance,stimulationFrame,basalLength);

    corrMat.(fnstr).foldchange = valout.foldchange;
    corrMat.(fnstr).difference = valout.difference;
    
    subfnames = fieldnames(corrMat.(fnstr));
    %abundance, foldchange, difference
    for subfnum = 1:length(subfnames)
        %scalar, rate, integral
        subfnstr = subfnames{subfnum};
        vals = corrMat.(fnstr).(subfnstr);
        [scalar,rate,integral] = rateintegralfunc(vals,stimulationFrame,timeMatrix);
        trueCorr.(fnstr).([subfnstr 'Xscalar']) = scalar;
        trueCorr.(fnstr).([subfnstr 'Xrate']) = rate;
        trueCorr.(fnstr).([subfnstr 'Xintegral']) = integral;
        
    end
end
    


%% find peaks
alltraces = corrMat.('Smad').foldchange;
% alltraces = cellMat.('Smad').all;
class1 = [];
class2 = [];
class3 = [];
tmat1 = [];
tmat2 = [];
tmat3 = [];
idx1 = [];
idx2 = [];
idx3 = [];
for j = 1:size(alltraces,1)
   cellvec = alltraces(j,:);
   tvec = timeMatrix(j,:);
   
   
   trange = tvec>0 & tvec<32;
   tval1 = tvec(trange);
   c1 = cellvec(trange);
   [fc1,i1] = max(c1);
   
   trange2 = tvec>60 & tvec<400;
   tval2 = tvec(trange2);
   c2 = cellvec(trange2);
   [fc2,i2] = max(c2);
   

   if ~isempty(fc1) && ~isempty(fc2)
      if fc1>2 && fc2>(fc1*1.1)
          class1 = vertcat(class1,cellvec);
          tmat1 = vertcat(tmat1,tvec);
          idx1 = vertcat(idx1,j);
      elseif ((fc1>2.7)&&(fc1<4.5)) && fc2>(fc1*0.8)
          class2 = vertcat(class2,cellvec);
          tmat2 = vertcat(tmat2,tvec);
          idx2 = vertcat(idx2,j);
      elseif fc1>2
          class3 = vertcat(class3,cellvec);
          tmat3 = vertcat(tmat3,tvec);
          idx3 = vertcat(idx3,j);
      end
       
   end
   
   
%    figure(44)
%    plot(tvec,cellvec);
%    hold on;
%    scatter(tval1(i1),fc1);
%    scatter(tval2(i2),fc2);hold off
%    stopherwe=1;
end
% 

reptraces = cellMat.('Reporter').all;

[maxrepvals,maxidx] = max(reptraces,[],2);
max1 = maxrepvals(idx1);
max2 = maxrepvals(idx2);
max3 = maxrepvals(idx3);
% subplot(1,3,1);
% plot(tmat1',class1');
% subplot(1,3,2);
% plot(tmat2',class2');
% subplot(1,3,3);
% plot(tmat3',class3');
% f=gcf;
% for h = f.Children'
%     h.YLim = [0 15];
% end

%% make figure


%set figure position
f = figure(22);
f.Units = 'pixels';
scrnsize = get(0,'screensize');
figPos = scrnsize;
figPos(1) = 1;
figPos(2) = 1;
figPos(3) = figPos(3);
figPos(4) = figPos(4);
f.Position = figPos;
f.Units = 'normalized';
f.Color = 'w';
    

tmat = timeMatrix;
toteOrMed = {SmadTotalOrMedian,ReporterTotalOrMedian};
traceArray = {'Smad','Reporter'};
tracenormArray = {smadnormalizestr,reporternormalizestr};
for tracenum = 1:length(traceArray)
    figure(22)
    tracestr = traceArray{tracenum};
    toteOrMedstr = toteOrMed{tracenum};
    tracenormstr = tracenormArray{tracenum};
    
    alltraces = cellMat.(tracestr).all;
    celltraces = cellMat.(tracestr).traces;
    divtraces = cellMat.(tracestr).division;
    tracesforSort = cellMat.('Smad').all;
    
    celltracemat = alltraces;
    divtracemat = divtraces;
%     [smoothtraces,smoothdivtraces] = smoothforplot(celltraces,divtraces,tmat);
    [smoothtraces,smoothdivtraces] = smoothforplot(alltraces,divtraces,tmat);
    celltracemat = smoothtraces;
    divtracemat = smoothdivtraces;
    
    
    subplot(2,3,1 + 3*(tracenum-1));
    cellmat = celltracemat(idx1,:);
    timemat = tmat(idx1,:);
    plot(timemat',cellmat')
    
    subplot(2,3,2 + 3*(tracenum-1));
    cellmat = celltracemat(idx2,:);
    timemat = tmat(idx2,:);
    plot(timemat',cellmat')
    
    subplot(2,3,3 + 3*(tracenum-1));
    cellmat = celltracemat(idx3,:);
    timemat = tmat(idx3,:);
    plot(timemat',cellmat')
    
    cellmat1 = celltracemat(idx1,:);
    timemat1 = tmat(idx1,:);
    cellmat2 = celltracemat(idx2,:);
    timemat2 = tmat(idx2,:);
    cellmat3 = celltracemat(idx3,:);
    timemat3 = tmat(idx2,:);
    

    
    f111 = figure(111);
    subplot(1,2,tracenum)
    [top1,bot1,med1] = booterfunc(cellmat1);
    [top2,bot2,med2] = booterfunc(cellmat2);
    [top3,bot3,med3] = booterfunc(cellmat3);
    
    cmap = lines(3);

    p3 = plot(timemat3(1,:),med3);hold on
    p3.LineWidth = 3;
    p3.Color = cmap(3,:);
    p3.DisplayName = 'class3';
    p1 = plot(timemat1(1,:),med1);
    p1.LineWidth = 3;
    p1.DisplayName = 'class1';
    p1.Color = cmap(1,:);
    p2 = plot(timemat2(1,:),med2);
    p2.LineWidth = 3;
    p2.DisplayName = 'class2';
    p2.Color = cmap(2,:);

    
    
    e = errorbar(timemat3(1,:),med3,bot3,top3);hold on
    e.Color = cmap(3,:);
    e.LineWidth = 0.5;
    e = errorbar(timemat1(1,:),med1,bot1,top1);hold on
    e.Color = cmap(1,:);
    e.LineWidth = 0.5;
    e = errorbar(timemat2(1,:),med2,bot2,top2);hold on
    e.Color = cmap(2,:);
    e.LineWidth = 0.5;


    nval1 = num2str(length(idx1));
    nval2 = num2str(length(idx2));
    nval3 = num2str(length(idx3));
    
    h=gca;
    h.Title.String = tracestr;
    h.XLabel.String = 'minutes';
    h.XLim = [-120 900];
    h.YLabel.String = {[tracestr ' ' tracenormstr],'median +/- 80%CI'};
    l = legend(h,[p1 p2 p3],['class1 (n=' nval1 ')'],['class2 (n=' nval2 ')'],['class3 (n=' nval3 ')']);
    
    
    %% determine pvalue difference
%     figure(999)
%     for k = 1:size(cellmat2,2)
%         x = cellmat1(:,k);
%         y = cellmat2(:,k);
%         z = cellmat3(:,k);
%         bootiter=100;
%         for bi = 1:bootiter
%             xmat = datasample(x,length(x),1,'Replace',true);
%             ymat = datasample(y,length(y),1,'Replace',true);
%             zmat = datasample(z,length(z),1,'Replace',true);
%             [~,pi12(bi)] = ttest2(xmat,ymat,'Vartype','unequal');
%             [~,pi13(bi)] = ttest2(xmat,zmat,'Vartype','unequal');
%             [~,pi23(bi)] = ttest2(ymat,zmat,'Vartype','unequal');
%         end
%         p12(k) = mean(pi12);
%         p13(k) = mean(pi13);
%         p23(k) = mean(pi23);
%     end
%     
% %     [h,p] = ttest2(cellmat1,cellmat2,'Vartype','unequal');
%     subplot(1,2,tracenum);
%     plot(timemat1(1,:),p12,'DisplayName','c1 vs c2');hold on
%     plot(timemat1(1,:),p13,'DisplayName','c1 vs c3');hold on
%     plot(timemat1(1,:),p23,'DisplayName','c2 vs c3');hold on
%     xlabel('minutes')
%     ylabel('Welchs ttest pvalue')
%     l = legend('show');
end






tmat = timeMatrix;

end


function cmaplz = colorcubemodified(cmaplength,cmapstr)
%default darkmin = 0.2
%default brightmax = 2
darkmin = 0.5;
colorationmin = 0.4;
brightmax = 2;

%generate colormap based on number of cells tracked
cnew=[];
%                 ccc = vertcat(colormap('summer'),colormap('autumn'),colormap('winter'),colormap('spring'));
%                 ccc = vertcat(colormap('hsv'),colormap('hot'));
ccc = flipud(colormap(cmapstr));
cmapccc = zeros(1000,3);
for j = 1:size(cmapccc,2)
    x = linspace(0,1,size(ccc,1));
    v = ccc(:,j);
    xq = linspace(0,1,1000);
    cmapccc(1:length(xq),j) = interp1(x,v,xq);
end

%                 ccc = colormap(jet(1000));
cccyc = 0;
for k = 1:size(cmapccc,1)
    cvec = cmapccc(k,:);
    if sum(cvec)>darkmin && sum(abs(diff(cvec)))>colorationmin && sum(cvec)<brightmax
        cccyc = cccyc+1;
        cnew(cccyc,:) = cvec;
    end
end

ccnew = zeros(cmaplength,3);
for j = 1:size(ccnew,2)
    x = linspace(0,1,size(cnew,1));
    v = cnew(:,j);
    xq = linspace(0,1,cmaplength);
    ccnew(1:length(xq),j) = interp1(x,v,xq);
end

cmaplz = ccnew;
end

function channelinputs =channelregexpmaker(channelstoinput)
channelinputs = '(';
for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
    if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i});
    elseif i < length(channelstoinput)
        channelinputs = strcat(channelinputs,'|',channelstoinput{i});
    else
        channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
    end
end
end

function cT = extractTraces(exportStruct,indices,speciesstr,finalFrame,stimulationFrame,basalLength,medfiltnum)
cT = struct();
xtstrArray = {'nuccyto','median','mean','total'};
for xtnum = 1:length(xtstrArray)
    %define how (median, mean, total) and what to quantify (EGFP, RFP)
    xtstr = xtstrArray{xtnum};
    xTracesString = [xtstr speciesstr];
    
    % extract the cell traces for the desired number of frames
    cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
    cellTracesUnfilt = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]
    % cellTraces = medfilt1(cellTracesUnfilt,medfiltnum,[],2,'omitnan','truncate'); %use median filtering
    % cellTraces = smooth(cellTracesUnfilt,0.1,'sgolay'); %use median filtering
    cellTraces = movmean(cellTracesUnfilt,medfiltnum,2,'Endpoints','shrink');
    cellTraces = cellTracesUnfilt;
    
    % normalize by basal values
    if (stimulationFrame-basalLength)<1
        basalLength = stimulationFrame-1;
    end
    %determine cell traces abundance (normalize to basal population median)
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    cellTracesAbundance = cellTraces./nanmedian(basalValues(:));
    
    %determine cell traces difference (subtract basal values of individual cells)
    basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
    basalVector = nanmedian(basalValues,2);
    invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
    invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
    basalMatrix = ones(size(cellTraces,2),1)*basalVector';
    cellTracesDifference = cellTracesAbundance - basalMatrix';
    
    %determine cell traces fold-change
    cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');
    
    %assign cellTraces to structure for easy access to data.
    cT.(xtstr).none = cellTraces;
    cT.(xtstr).abundance = cellTracesAbundance;
    cT.(xtstr).difference = cellTracesDifference;
    cT.(xtstr).foldchange = cellTracesFoldChange;
end
end
function cT = extractTracesNC(exportStruct,indices,speciesstr,finalFrame,stimulationFrame,basalLength,medfiltnum)
cT = struct();
xtstrArray = {'median'};
for xtnum = 1:length(xtstrArray)
    %define how (median, mean, total) and what to quantify (EGFP, RFP)
    xtstr = xtstrArray{xtnum};
    xTracesString = [xtstr speciesstr];
    
    % extract the cell traces for the desired number of frames
    cellTracesFull = vertcat(exportStruct(indices).(xTracesString));
    cellTracesUnfilt = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]
    % cellTraces = medfilt1(cellTracesUnfilt,medfiltnum,[],2,'omitnan','truncate'); %use median filtering
    % cellTraces = smooth(cellTracesUnfilt,0.1,'sgolay'); %use median filtering
    cellTraces = movmean(cellTracesUnfilt,medfiltnum,2,'Endpoints','shrink');
    cellTraces = cellTracesUnfilt;
    
    % normalize by basal values
    if (stimulationFrame-basalLength)<1
        basalLength = stimulationFrame-1;
    end
    %determine cell traces abundance (normalize to basal population median)
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    cellTracesAbundance = cellTraces./nanmedian(basalValues(:));
    
    %determine cell traces difference (subtract basal values of individual cells)
    basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
    basalVector = nanmedian(basalValues,2);
    invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
    invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
    basalMatrix = ones(size(cellTraces,2),1)*basalVector';
    cellTracesDifference = cellTracesAbundance - basalMatrix';
    
    %determine cell traces fold-change
    cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');
    
    %assign cellTraces to structure for easy access to data.
    cT.(xtstr).none = cellTraces;
    cT.(xtstr).abundance = cellTracesAbundance;
    cT.(xtstr).difference = cellTracesDifference;
    cT.(xtstr).foldchange = cellTracesFoldChange;
end
end

function allidx = jittercorrfunc(pmat,jittercut)
%jitter correction in reporter
jittercorrectMat = pmat;
diffrmatz = diff(jittercorrectMat,1,2);
diffrmat = [diffrmatz(:,1) diffrmatz];
upidx = diffrmat>jittercut;
downidx = diffrmat<(-1.*jittercut);
allidx = upidx | downidx;
end

function cmap = determineCMAP(colormapChoice,uniqueColoring)
cmapBig = feval(colormapChoice,length(uniqueColoring));
cmapidx = floor((linspace(1,size(cmapBig,1),length(uniqueColoring))));
cmap = cmapBig(cmapidx,:);
end

function [valList,valArrayList] = listSpitter(exportSubStruct,valBasedOn)
    indices = true(1,length(exportSubStruct));
    if ~iscell(valBasedOn)
        valArrayList = vertcat({exportSubStruct(indices).(valBasedOn)});
        valList = unique(valArrayList);
    else
        for i = 1:length(indices)
            concatstr = [];
            for j = 1:length(valBasedOn)
                cstr = valBasedOn{j};
                newstr = exportSubStruct(i).(cstr);
                concatstr = [concatstr ' - ' newstr];
            end
            valArrayList{i} = concatstr;
        end
        valList = unique(valArrayList);
    end
end

function [tracemat,divtracemat,allmat,nanidx] = divisionFunction(cTracks,smadtraces,stimulationFrame,timeMatrix)
celldivtraces = cTracks.Smad.total.foldchange;
repdivtraces = cTracks.Reporter.total.foldchange;

celldivtracesMed = cTracks.Smad.median.foldchange;
repdivtracesMed = cTracks.Reporter.median.foldchange;


nanidx = ~isnan(celldivtraces(:,stimulationFrame));

celltraces = smadtraces(nanidx,:);
tmat = timeMatrix(nanidx,:);
rmat = repdivtraces(nanidx,:);
pmat = celldivtraces(nanidx,:);
rmatmed = repdivtracesMed(nanidx,:);
pmatmed = celldivtracesMed(nanidx,:);
divtracemat = nan(size(rmat));
tracemat = nan(size(rmat));
allmat = nan(size(rmat));

for i = 1:size(rmat,1)
    tvec = round(tmat(i,:));
    pvec = pmat(i,:);
    rvec = rmat(i,:);
    cellvec = celltraces(i,:);
    cellnan = isnan(cellvec);
    ps = smooth(pvec,0.05,'sgolay');
    rs = smooth(rvec,0.05,'sgolay');
    pvec = ps;
    rvec = rs;
    pvec(cellnan) = NaN;
    rvec(cellnan) = NaN;
    
    rdiff = zeros(1,size(rmat,2));
    pdiff = zeros(1,size(rmat,2));
    for j = 1:size(rmat,2)
        afound = find(tvec > tvec(j)-20,1,'first');
        bfound = find(tvec < tvec(j)+20,1,'last');
        a = max([1 afound]);
        b = min([bfound size(rmat,2)]);
        rdiff(1,j) = rvec(b)./rvec(a);
        pdiff(1,j) = pvec(b)./pvec(a);
    end
    pdiff(cellnan) = NaN;
    rdiff(cellnan) = NaN;
    
    
    
    se = strel('disk',2);
    pdiffless = imdilate(pdiff<0.7,se);
    rdiffless = imdilate(rdiff<0.7,se);
    divpoint = pdiffless & rdiffless;
    CC = bwconncomp(divpoint);
    PX = CC.PixelIdxList;
    imglog = false(size(rvec));
    for k = 1:length(PX)
        px = PX{k};
        pxlog = false(size(rvec));
        pxlog(px) = true;
        se = strel('disk',9);
        pxdilate = find(imdilate(pxlog,se));
        rtest = rdiff(pxdilate);
        idx = rtest>1.1 | rtest<0.9;
        pxnew = pxdilate(idx);
        imglog(pxnew) = true;
    end
    se =strel('disk',1);
    imglogc = imclose(imglog,se);
    CC = bwconncomp(imglogc);
    PX = CC.PixelIdxList;
    pxlog = cellfun(@(x) length(x)>1,PX,'UniformOutput',1);
    PXnew = PX(pxlog);
    CC.PixelIdxList = PXnew;
    CC.NumObjects = length(PXnew);
    
    
    %now interpolate through the dividing region
    cellvec = celltraces(i,:);
    cellvectest = cellvec;
    cellvectest(isnan(pdiff))=NaN;
    tracevec = nan(size(rvec));
    divvec = nan(size(rvec));
    pxstartmat = [1 cellfun(@(x) min([x(end)+1 length(rvec)]),PXnew,'UniformOutput',1)];
    pxendmat = [cellfun(@(x) max([1 x(1)-1]),PXnew,'UniformOutput',1) length(rvec)];
    for k = 1:length(pxstartmat)
        pxstart = pxstartmat(k);
        pxend = pxendmat(k);
        if k == 1
            adjustval =0;
        else
            adjustval = tracevec(pxendmat(k-1))-cellvec(pxstart);
            adjustval = 0;
        end
        tracevec(pxstart:pxend) = cellvec(pxstart:pxend)+adjustval;
        if k>1
            c = min([max([2 pxendmat(k-1)]) length(rvec)-1]);
            d = max([2 pxstartmat(k)]);
            x = [c-1 c d d+1];
            v = [cellvectest(c-1) cellvectest(c) cellvectest(d) cellvectest(d+1)];
            xq = c:1:d;
            interpp = interp1(x,v,xq);
            divvec(c:d) = interpp;
        end
    end
    allvec = tracevec;
    allvec(~isnan(divvec)) = divvec(~isnan(divvec));
    tracemat(i,:) = tracevec;
    divtracemat(i,:) = divvec;
    allmat(i,:) = allvec;
%             figure(2)
%             plot(tvec,cellvec,'LineWidth',1);hold on
%             plot(tvec,pvec,'LineWidth',1);hold on
%             plot(tvec,rvec,'LineWidth',1);hold on
%             plot(tvec,tracevec,'LineWidth',2);hold on
%             plot(tvec,divvec,'k:','LineWidth',2);hold off
%             xlim([-100 1400])
%             ylim([0 4])
%             if i>33
%                 stophere=1;
%             end
end
end

function [smoothtraces,smoothdivtraces] = smoothforplot(celltraces,divtraces,tmat)
   
forsort1 = celltraces;
forsort1(~isnan(divtraces)) = divtraces(~isnan(divtraces));
tracesforSort = forsort1;

sm = nan(size(tracesforSort));
for i = 1:size(sm,1)
    tvec = tmat(i,:);
    fsvec = tracesforSort(i,:);
    smvec = smooth(tvec,fsvec,0.1,'loess');
    sm(i,:) = smvec;
end
sm(isnan(tracesforSort))=NaN;

smoothdivtraces = nan(size(divtraces));
smoothdivtraces(~isnan(divtraces)) = sm(~isnan(divtraces));

smoothtraces = sm;
smoothtraces(~isnan(divtraces)) = NaN;

% figure(9)
% subplot(1,2,1);
% plot(tracesforSort');ylim([0 4])
% subplot(1,2,2);
% plot(sm'); ylim([0 4])
% sss=1;
end

function [sortedCellMat,sortidx] = sorttraces(cellMat,sortstr,tmat,sortTime)
% celltraces = cellMat.(sortstr).traces;
% divtraces = cellMat.(sortstr).division;
% forsort1 = celltraces;
% forsort1(~isnan(divtraces)) = divtraces(~isnan(divtraces));
alltraces = cellMat.(sortstr).all;
tracesforSort = alltraces;
tvec = tmat(1,:);


%determine sortidx
sortFrame = find(tvec(1,:)>sortTime);
if isempty(sortFrame)
    disp('sort frame error')
    sortFrame = stimulationFrame+1;
end
spvecForSorting = tracesforSort(:,sortFrame(1));
[~,sortidx] = sort(spvecForSorting);

sortedCellMat = struct();
fnames = fieldnames(cellMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    subMat = cellMat.(fnstr);
    subfnames = fieldnames(subMat);
    for subfnum = 1:length(subfnames)
        subfnstr = subfnames{subfnum};
        subtraces = subMat.(subfnstr);
        
        sortedsubcelltraces = zeros(size(tracesforSort));
        for k = 1:length(sortidx)
            sortedsubcelltraces(k,:) = subtraces(sortidx(k),:);
        end
        subtraces = sortedsubcelltraces;
        subMat.(subfnstr) = subtraces;
    end
    sortedCellMat.(fnstr) = subMat;
    
end
end


function newCellMat = updatealltraces(cellMat,timenan)

newCellMat = cellMat;
fnames = fieldnames(cellMat);
for fnum = 1:length(fnames)
    fnstr = fnames{fnum};
    cellfn = cellMat.(fnstr);
    
    fnames2 = fieldnames(cellfn);
    for fnum2 = 1:length(fnames2)
        fnstr2 = fnames2{fnum2};
        
        cellfn2 = cellfn.(fnstr2);
        cellnew = cellfn2(timenan,:);
        newCellMat.(fnstr).(fnstr2) = cellnew;
    end
end
end

%% new functions
function cT = fcanddiff(abundance,stimulationFrame,basalLength)
cT = struct();
cellTraces = abundance;

% normalize by basal values
if (stimulationFrame-basalLength)<1
    basalLength = stimulationFrame-1;
end
%determine cell traces abundance (normalize to basal population median)
basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
cellTracesAbundance = cellTraces./nanmedian(basalValues(:));

%determine cell traces difference (subtract basal values of individual cells)
basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
basalVector = nanmedian(basalValues,2);
invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
basalMatrix = ones(size(cellTraces,2),1)*basalVector';
cellTracesDifference = cellTracesAbundance - basalMatrix';

%determine cell traces fold-change
cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');

%assign cellTraces to structure for easy access to data.
cT.none = cellTraces;
cT.abundance = cellTracesAbundance;
cT.difference = cellTracesDifference;
cT.foldchange = cellTracesFoldChange;
end


function [scalar,rate,integral] = rateintegralfunc(vals,stimulationFrame,timeMatrix)
%scalar
scalar = vals;

%rate
[rate,~] = gradient(vals,timeMatrix(1,:),timeMatrix(:,1));

%integral
preintegral = cumsum(vals,2,'omitnan');
preintegral(isnan(vals)) = NaN;
basalValues = preintegral(:,stimulationFrame);
basalVector = nanmedian(basalValues,2);
basalMatrix = ones(size(preintegral,2),1)*basalVector';
integral = preintegral - basalMatrix';
end

function [rho,nvalue] = determineCorrelation(pMat,rMat,corrtype,stimulationFrame)
rho = zeros(size(pMat,2),size(rMat,2));
nvalue = zeros(size(pMat,2),size(rMat,2));
    for st = stimulationFrame:size(pMat,2) %number of frames
        x = pMat(:,st);
        for rt = st:size(rMat,2)
            y = rMat(:,rt);
            xn = ~isnan(x);
            yn = ~isnan(y);
            cidx = xn&yn;
            xi = x(cidx);
            yi = y(cidx);
            
%             xtf = ~isoutlierz(xi,'quartiles');
%             ytf = ~isoutlierz(yi,'quartiles');
%             oidx = xtf&ytf;
            oidx = true(size(xi));
            xu = xi(oidx);
            yu = yi(oidx);
            
            if isempty(yu) || isempty(xu)
                r =0;
            else
                r = corr(xu,yu,'type',corrtype,'rows','all');
            end
            rho(st,rt) = r;
            nvalue(st,rt) = length(xu);
        end
    end
%minimum number of data points for considering correlation
nidx  = nvalue<30;
rho(nidx) = 0;
end

function tf = isoutlierz(xi,funame)
    q1 = prctile(xi,25);
    q2 = prctile(xi,50);
    q3 = prctile(xi,75);
    iqr = q3-q1;
    outlierdist = iqr*1.5;
    lowcut = q1-outlierdist;
    highcut = q3+outlierdist;
    
    tf1 = xi>lowcut;
    tf2 = xi<highcut;
    tf = ~(tf1&tf2);
end

function [top1,bot1,med1] = booterfunc(cellmat1)
% nanmat = ~isnan(cellmat1);
% nval = sum(nanmat,1);
% nvalidx = nval>7;
% cellmat1(:,~nvalidx) = NaN;
% 
% botvals = bootstrp(100,@(x) prctile(x,10,1),cellmat1);
% topvals = bootstrp(100,@(x) prctile(x,90,1),cellmat1);
% medvals = bootstrp(100,@(x) prctile(x,50,1),cellmat1);
% 
% 
% med1 = nanmean(medvals,1);
% top1 = nanmean(topvals,1)-med1;
% bot1 = med1 - nanmean(botvals,1);


nanmat = ~isnan(cellmat1);
nval = sum(nanmat,1);
nvalidx = nval>6;
cellmat1(:,~nvalidx) = NaN;

% botvals = bootstrp(100,@(x) prctile(x,10,1),cellmat1);
% topvals = bootstrp(100,@(x) prctile(x,90,1),cellmat1);
% medvals = bootstrp(100,@(x) prctile(x,50,1),cellmat1);

botvals = prctile(cellmat1,10,1);
topvals = prctile(cellmat1,90,1);
medvals = prctile(cellmat1,50,1);


med1 = nanmean(medvals,1);
top1 = nanmean(topvals,1)-med1;
bot1 = med1 - nanmean(botvals,1);
end
