function script_plot_tracked_cells_hypothesis
close all


%set directory to location of code being used (generally external harddrive
%%%%
%determine path to .m file being executed
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b);
    exportdir = strcat(parentdir,'Export/');

%determine path to gparent folder
    [~,b ] = regexp(parentdir,'/');
        if isempty(b)
            [~,b] = regexp(parentdir,'\');
        end
        gparentdir = parentdir(1:b(end-1));

    cd(parentdir)
    cd ..


%%

% includeArray = {'(clone69|clone80)'};
    includeStr = '2';
%     excludeStr = '(cond|media|->)';
%     excludeStr = 'sb';
    excludeStr = '(exclude|excluder)';
    
%     datestr= '2017_07_08 plate exp1';
%     datestr= '2017_07_07 plate exp1';
%     datestr= '2017_07_01 plate exp2';
%     datestr= '2017_06_26 plate exp2';
    datestr= '2017_04_17 plate exp4';
%     datestr= '2017_07_01 plate exp1';
    
    
    
    FileName = [datestr '_tracking_export.mat'];
    
    %choose from the following (how to quantify nuclei)
    %total, median
    SmadTotalOrMedian = 'median'; % median, mean, total
    ReporterTotalOrMedian = 'total';
    
    %choose from the following
    %doseconddate, conddate, doseAndCondition, conditions, expdate,
    %dosestr, scene, wells, cellid
    groupBasedOn = {'scene','cellID'};
    labelBasedOn = {'scene','cellID'};
    subplotGroupBasedOn = {'scene','cellID'};
    colorBasedOn = {'scene','cellID'}; % can be 'individually' or {'cellID','scene'}
%     colorBasedOn = {'doseAndCondition','scene','cellID'};
    
    cmapchoice = 'lines';

    celltracestr = 'Reporter'; %Smad, Reporter
    smadnormalizestr = 'foldchange'; %none, abundance, difference, or foldchange, nuccytoXabundance, nuccyto, nuccytoXabundanceFC
    reporternormalizestr = 'difference';
    sortTime = 100; %minutes
    timeCutoff = 200;
    basalLength = 4;
%     determineExportNucleiStructCompiled(FileName,colorBasedOn)
    plot_tracked_cells_hypothesis(FileName,timeCutoff,sortTime,cmapchoice,groupBasedOn,colorBasedOn,subplotGroupBasedOn,labelBasedOn,basalLength,smadnormalizestr,reporternormalizestr,includeStr,excludeStr,SmadTotalOrMedian,ReporterTotalOrMedian,celltracestr)
    stopre=1;
end