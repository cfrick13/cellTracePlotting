function script_nucleiPlotterIter
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

includeArray = {'clone69','clone80','(clone69|clone80)','clone116','clone119','(clone116|clone119)','clone86','clone88','(clone86|clone88)'};
includeArray = {'(clone116|clone119)'};
    for i = 1:length(includeArray)
    includeStr = includeArray{i};
    excludeStr = '(cond|media|->|0.18|0.06|0.02)';
    % excludeStr = '(cond|media|->|sb)';
    FileName=[];
    %choose from the following
    %doseconddate, conddate, doseAndCondition, conditions, expdate, dosestr, scene, wells
    groupBasedOn = 'doseconddate';
    normalizeBasedOn = {'conditions','expdate'};
    labelBasedOn = 'dosestr';
    subplotGroupBasedOn = 'conddate';
    colorBasedOn = 'dosestr';

    toteORmedReporter = 'median'; %how to quantify nuclei
    normalizestr = 'difference'; %abundance, difference, or foldchange
    basalLength = 4;
%     determineExportNucleiStructCompiled(FileName,colorBasedOn)
    nucleiPlotterNew(FileName,groupBasedOn,colorBasedOn,subplotGroupBasedOn,normalizeBasedOn,labelBasedOn,basalLength,toteORmedReporter,normalizestr,includeStr,excludeStr)
    close all
    end
end