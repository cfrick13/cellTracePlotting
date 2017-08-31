function nucleiPlotterNew(FileName,groupBasedOn,colorBasedOn,subplotGroupBasedOn,normalizeBasedOn,labelBasedOn,basalLength,toteORmedReporter,normalizestr,includeStr,excludeStr)

%determine the location of the matlab function and establish export 
%directory in relation to that filepath
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b); %specifies folder in which all analysis is being done
    loaddir = strcat(parentdir,'Export'); %specifies where data is exported
    cd(loaddir);

    [~,b] = regexp(mdir,'/');
            if isempty(b)
                [~,b] = regexp(mdir,'\');
            end

    mfiledir =mdir(1:b(end)); %specifies location of matlab function file

   
    cd(mfiledir)
    tic
    load('exportNucleiStructCompiled.mat') 
    toc

            
   %choose all cells that match specified chosen clone input strings
    allList = {exportNucleiStructCompiled.('doseconddate')};
    [~,~,~,d] = regexp(allList,includeStr);
    subidx = ~cellfun(@isempty,d,'UniformOutput',1);
    exportNucleiStruct = exportNucleiStructCompiled(subidx);

    %exclude all cells that match specified exludeNames input strings
    allList = {exportNucleiStruct.('doseconddate')};
    [~,~,~,d] = regexp(allList,excludeStr);
    subidx = ~cellfun(@isempty,d,'UniformOutput',1);
    exportNucleiStruct = exportNucleiStruct(~subidx);

    %choose which field based upon which each cell trace will get colored
        coloringChoice = colorBasedOn; %color traces based on what
        colormapChoice = 'jet'; %colormap to be used
        darkenFactor = 1;
    %establish the color map for plotting
        coloringArray = vertcat({exportNucleiStruct.(coloringChoice)});
        indices = true(1,length(exportNucleiStruct));
        coloringArrayTrunc = vertcat({exportNucleiStruct(indices).(coloringChoice)});
        uniqueColoring = unique(coloringArrayTrunc);
        figure(1)
        cmapBig = colormap(colormapChoice)./darkenFactor;
        cmapidx = floor((linspace(1,size(cmapBig,1),length(uniqueColoring))));
        cmap = cmapBig(cmapidx,:);
        close 1


    %assign a color array using the created colormap based on the choices above
        colormapMatrix = zeros(length(coloringArray),size(cmap,2));
        for i=1:length(coloringArray)
           cA = coloringArray{i};
           idx = strcmp(uniqueColoring,cA);
           colormapMatrix(i,:) = cmap(idx,:);
        end

        cd(mfiledir)
%% 
%compile a new structure (newstruct) based on specified criteria, groupBasedOn
    str = groupBasedOn;
    criteriaArrayList = {exportNucleiStruct.(str)};
    criteriaList = unique(criteriaArrayList);
    criteriaStruct = struct();
    fnames = fieldnames(exportNucleiStruct);
    excludeNames = '(scene|cellID|dose|tgfFrame|tgfFramestr|flatfield|wells|conditions|flatfielddosestr|framestr|Centroid|expdate|timeMatrix|conddate)';
    [~,~,~,d] = regexpi(fnames,excludeNames);
    subidx = ~cellfun(@isempty,d,'UniformOutput',1);
    fnames(subidx)=[];
    framelength = zeros(1,length(criteriaList));
    for d = 1:length(criteriaList)
        dname = criteriaList{d};
        idx = strcmp(criteriaArrayList,dname);
        nidx = find(idx==1);
        
        %%%
        timeMatrixNew = exportNucleiStruct(nidx(1)).timeMatrix;
        finalFrame = length(timeMatrixNew);
        framelength(d) = length(timeMatrixNew);
        %%%
        
        
        arraycat = cell(length(nidx),finalFrame);

        for fi = 1:length(fnames)
            fnamestr = fnames{fi};
            for j = 1:length(nidx)
                jn = nidx(j);
                arraycat(j,:)  = exportNucleiStruct(jn).(fnamestr);
            end
            arraymatcat = cell(1,finalFrame);
            for j = 1:finalFrame
                arraymatcat{j} = [arraycat{:,j}];
            end
            criteriaStruct(d).(fnamestr) = arraymatcat;
        end
    end
    
maxlength = max(framelength);
finalFrame = maxlength;

%% compile all traces into one mat (based on grouping) //used for determining the axis limits
%assemble compiled matrices
    str = groupBasedOn;
    criteriaArrayList = {exportNucleiStruct.(str)};
    criteriaList = unique(criteriaArrayList);

    
    subplotArrayList = {exportNucleiStruct.(subplotGroupBasedOn)};
    subplotList = unique(subplotArrayList);
    
    labelArrayList = {exportNucleiStruct.(labelBasedOn)};
    labelList = unique(labelArrayList);
    
    str = 'conditions';
    conditionsArrayList = {exportNucleiStruct.(str)};
    conditionsList = unique(conditionsArrayList);
    
    str = 'doseAndCondition';
    doseAndConditionsArrayList = {exportNucleiStruct.(str)};
    doseAndConditionsList = unique(doseAndConditionsArrayList);
    
    str = 'doseconddate';
    doseconddateArrayList = {exportNucleiStruct.(str)};
    doseconddateList = unique(doseconddateArrayList);
    
    str = 'expdate';
    expdateArrayList = {exportNucleiStruct.(str)};
    normGroupList = unique(expdateArrayList);
    
    
    compiledStruct = struct();
    fnames = fieldnames(exportNucleiStruct);
    [~,~,~,d] = regexp(fnames,excludeNames);
    subidx = ~cellfun(@isempty,d,'UniformOutput',1);
    fnames(subidx)=[];
    for fi = 1:length(fnames)
        
        fnamestr = fnames{fi};
        compiledTime = nan(length(criteriaList),finalFrame);
        compiledMatrix = nan(length(criteriaList),finalFrame);
        compiledMatrixSTD = nan(length(criteriaList),finalFrame);
        colorbasedonmap = zeros(length(criteriaList),3);
        subplotMatchList = cell(1,length(criteriaList));
        labelMatchList = cell(1,length(criteriaList));
            for d = 1:length(criteriaList)
                dname = criteriaList{d};
                idx = strcmp(criteriaArrayList,dname);
                subplotMatchList{d} = unique(subplotArrayList(idx)); 
                labelMatchList{d} = unique(labelArrayList(idx)); 
                
                nidx = find(idx==1);
                
                %%%
                timeMatrixNew = exportNucleiStruct(nidx(1)).timeMatrix;
                finalFrame = length(timeMatrixNew);
                tframeidx = true(1,finalFrame);
                %%%
                        
                arraycat = cell(length(nidx),finalFrame);
                for j = 1:length(nidx)
                    jn = nidx(j);
                    criteriaArray  = exportNucleiStruct(jn).(fnamestr);
                    criteriaArrayT = criteriaArray(tframeidx);
                    arraycat(j,:)  = criteriaArrayT;
                end
                
                arraymatcat = cell(1,finalFrame);
                for j = 1:finalFrame
                    arraymatcat{j} = [arraycat{:,j}];
                end
                    
                    
                timeMatVal = timeMatrixNew;
                stimulationFramez = find(floor(timeMatVal)<0);
                stimulationFrame=stimulationFramez(end);
                
                criteriaMat = cellfun(@(x) nanmedian(single(x)),arraymatcat,'UniformOutput',1);
                if stimulationFrame-basalLength<0
                    basalLength = stimulationFrame-1;
                end
                
                basalValues = criteriaMat(:,(stimulationFrame-basalLength):stimulationFrame+1);
                basalVal = nanmedian(basalValues(:));
                
                if strcmp(normalizestr,'abundance')
                    criteriaMatVal = criteriaMat;
                    criteriaMatSTDval = cellfun(@(x) std(single(x))./sqrt(length(x)),arraymatcat,'UniformOutput',1);
                elseif strcmp(normalizestr,'difference')
                    criteriaMatVal = criteriaMat-nanmedian(basalVal);
                    criteriaMatSTDval = cellfun(@(x) std(x-basalVal)./sqrt(length(x)),arraymatcat,'UniformOutput',1);
                elseif strcmp(normalizestr,'foldchange')
                    criteriaMatVal = criteriaMat./nanmedian(basalVal);
                    criteriaMatSTDval = cellfun(@(x) std(x./basalVal)./sqrt(length(x)),arraymatcat,'UniformOutput',1);
                end

                compiledMatrix(d,1:finalFrame) = criteriaMatVal;
                compiledMatrixSTD(d,1:finalFrame) = criteriaMatSTDval;
                compiledTime(d,1:finalFrame) = timeMatVal;
                
                colorbasedonmap(d,:) = colormapMatrix(nidx(1),:);
            end
        compiledStruct.(fnamestr) = compiledMatrix;
        compiledStruct.([fnamestr 'std']) = compiledMatrixSTD;
        
    end
   



%%
toteOrMed = 'median';
smadTracesString = [toteOrMed 'Smad']; %value to plot
reporterTracesString = [toteORmedReporter 'Reporter'];

pMatrix = compiledStruct.(smadTracesString);
rMatrix = compiledStruct.(reporterTracesString);
pMatrixSTD = compiledStruct.([smadTracesString 'std']);
rMatrixSTD = compiledStruct.([reporterTracesString 'std']);
%determine the axis limits


    %determine the cells to group for normalization
    if isstr(normalizeBasedOn)
        str = normalizeBasedOn;
        normGroupArrayList = {exportNucleiStruct.(str)};
        normGroupList = unique(normGroupArrayList);
    else
        
        normGroupArrayList=[];
        for i = 1:length(normalizeBasedOn)
            str = normalizeBasedOn{i};
            normGroupArrayList = vertcat(normGroupArrayList,{exportNucleiStruct.(str)});
        end
        
        concatedGroupArrayList = cell(size(normGroupArrayList,2),1);
        for i = 1:size(normGroupArrayList,2)
            normGroup = normGroupArrayList(:,i);
            channelinputs = channelregexpmaker(normGroup);
            concatedGroupArrayList{i,:} = channelinputs;
        end
        
        normGroupList = unique(concatedGroupArrayList);
    end
    
    
    %normalize the cells to 0 to 1
    for enum = 1:length(normGroupList)
        normgroupstr = normGroupList{enum};
        [~,~,~,d] = regexp(criteriaList,normgroupstr);
        if isstr(normalizeBasedOn)
            subidx = ~cellfun(@isempty,d,'UniformOutput',1);
        else
            subidx = cellfun(@(x) length(x)==length(normalizeBasedOn),d,'UniformOutput',1);
        end
        
        %rescale data
        loperc = 0;
        hiperc = 100;
        rmat = rMatrix(subidx,:);
        rbottom = prctile(rmat(:),loperc);
        rtop = prctile(rmat(:),hiperc);
        scaleFactor = 1./(rtop - rbottom);
        rmatnew = rmat.*scaleFactor;
        rmatnew = rmatnew-(rbottom.*scaleFactor);
        rMatrix(subidx,:) = rmatnew;
        
        %rescale errorbars
        rmatstd = rMatrixSTD(subidx,:);
        rmatstdnew = rmatstd.*scaleFactor;
        rMatrixSTD(subidx,:) = rmatstdnew;
    end

    



    

    
%make the beautiful plot    
    plotCriteriaList = criteriaList;
    flcl = floor(sqrt(length(subplotList)));
    clcl = ceil(sqrt(length(subplotList)));
    if (flcl*clcl)<length(subplotList)
        flcl =ceil(sqrt(length(subplotList)));
    end
    
    for j = 1:length(plotCriteriaList)
        
        criteriastr = plotCriteriaList{j};
        [a,~] = regexp(criteriastr,'clone');
        if isempty(a)
            cloneName = criteriastr;
            criteriaName = criteriastr;
        else
            cloneName = criteriastr(a:end);
            criteriaName = criteriastr;
        end
        
            
        
        subplotstr = subplotMatchList{j};
        subidx =  strcmp(subplotList,subplotstr);
        subnum = find(subidx==1);  
        
        labelstr = char(labelMatchList{j});

            %%%
            timeMatrixNew = exportNucleiStruct(nidx(1)).timeMatrix;
            finalFrame = length(timeMatrixNew);
            tframeidx = true(1,finalFrame);
            %%%
        tMat = compiledTime(j,tframeidx);
        rMat = rMatrix(j,tframeidx);
        pMat = pMatrix(j,tframeidx);
        
        rMatSTD = rMatrixSTD(j,tframeidx);
        pMatSTD = pMatrixSTD(j,tframeidx);
        
        disp([num2str(j) '/' num2str(length(plotCriteriaList))])
        if sum(isnan(rMat))==length(rMat)
            rMat = zeros(size(rMat));
        end
%         rMat = smooth(rMat,3,'rlowess')';
%         pMat = smooth(pMat,10,'rlowess')';
       

%         f1 = figure(1);
%         sax2 = subplot(1,2,2);
%             e=errorbar(tMat,rMat,rMatSTD); hold on
% 
%                 e.Color = colorbasedonmap(j,:)./1.5;
%                 e.LineWidth = 1;
%                 e.Marker = 's';
%                 e.MarkerSize = 8;
%                 e.MarkerFaceColor = colorbasedonmap(j,:);
%                 e.DisplayName = criteriastr;
%             title('mCherry');
%             xlim(xlimits2)
% %             ylim(ylimits2)
%             ylabel('median reporter')
%             xlabel('minutes')
%             sax2.XGrid='on';
%             sax2.YGrid='on';
%             sax2.Color = [0.9 0.9 0.9];
        
%         sax4 = subplot(1,2,1);
%             e=errorbar(tMat,pMat,pMatSTD); hold on
%                 e.Color = colorbasedonmap(j,:)./1.5;
%                 e.LineWidth = 1;
%                 e.Marker = 's';
%                 e.MarkerSize = 8;
%                 e.MarkerFaceColor = colorbasedonmap(j,:);
%                 e.DisplayName = criteriastr;
%             title('Smad3');
% %             xlim(xlimits2)
% %             ylim(ylimits2)
%             ylabel('median reporter')
%             xlabel('minutes')
%             sax4.XGrid='on';
%             sax4.YGrid='on';
%             sax4.Color = [0.9 0.9 0.9];


        f1 = figure(1);
        sax2 = subplot(flcl,clcl,subnum);
            p=plot(tMat,rMat); hold on
%                 p.Color = colorbasedonmap(j,:)./1.5;
                p.Color = colorbasedonmap(j,:);
                p.LineWidth = 4;
                p.DisplayName = labelstr;
            title(subplotstr);
%             xlim(xlimits2)
%             ylim(ylimits2)
            ylabel('median reporter')
            xlabel('minutes')
            sax2.XGrid='on';
            sax2.YGrid='on';
            sax2.Color = [0.9 0.9 0.9];
            
%         sax4 = subplot(flcl,clcl,1);
%             p=plot(tMat,pMat); hold on
% %                 p.Color = colorbasedonmap(j,:)./1.5;
%                 p.Color = colorbasedonmap(j,:);
%                 p.LineWidth = 4;
%                 p.DisplayName = criteriastr;
%             title('Smad3');
% %             xlim(xlimits2)
% %             ylim(ylimits2)
%             ylabel('median reporter')
%             xlabel('minutes')
%             sax4.XGrid='on';
%             sax4.YGrid='on';
%             sax4.Color = [0.9 0.9 0.9];


    end
    
    

    

    for h = f1.Children'
        axis(h,'tight')
        ylimit = h.YLim;
        ylimit(1) = 0;
        ylimit(2) = 1;
        h.YLim = ylimit;
        h.XLim = [-100 1200];
        l = legend(h,'show');
        l.Location = 'eastoutside';
    end
    f1.Color='w';
    scrnsize = get(0,'screensize');
    figsize =scrnsize;
%     figsize(3) = figsize(3)./2;
    f1.Position = figsize;
    
    cd(mfiledir)
    dirnamestr = [mfiledir '/' 'compiled'];
    if ~isdir(dirnamestr)
        mkdir(dirnamestr)
    end
    olddir = pwd;
    cd(dirnamestr)
    savestr = includeStr;
    [a] = regexp(savestr,'('); savestr(a)=[]; a = regexp(savestr,'[|]'); savestr(a)='-';[a] = regexp(savestr,')'); savestr(a)=[];
    saveas(f1,[savestr ' ' normalizestr ' fig.fig']);
    saveas(f1,[savestr ' ' normalizestr ' image.png'],'png');
    cd(olddir)
end





function channelinputs = channelregexpmaker(channelstoinput)
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

