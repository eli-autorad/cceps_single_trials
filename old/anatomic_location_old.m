function ana = anatomic_location_old(pt_name,chLabels)

ana = cell(length(chLabels),1);


switch pt_name
    case 'HUP212_CCEP'
        for ich = 1:length(chLabels)
            label = chLabels{ich};
            
            % get the non numerical portion
            label_num_idx = regexp(label,'\d');
            if isempty(label_num_idx), continue; end
            
            label_non_num = label(1:label_num_idx-1);
            
            switch label_non_num
                case 'LC'
                    ana{ich} = 'HIPP';%'hippocampus';%'posterior\newlinehippocampus';
                case 'LB'
                    ana{ich} = 'HIPP';'hippocampus';%'anterior\newlinehippocampus';
                case 'LA'
                    ana{ich} = 'AM';%'amygdala';
                case 'LD'
                    ana{ich} = 'AI';%{'anterior','insula'};
                case 'LE'
                    ana{ich} = 'CS';%{'central','sulcus'};
                case 'LF'
                    ana{ich} = 'PI';%;{'posterior','insula'};
                case 'LG'
                    ana{ich} = 'PT';%{'pars','triangularis'};
                case 'LH'
                    ana{ich} = {'IFG'};
                case 'LI'
                    ana{ich} = 'AC';%{'anterior','cingulate'};
                case 'LJ'
                    ana{ich} = {'MFG'};
            end
            
        end 
    otherwise
        ana = [];
end


%{
switch pt_name
    case 'HUP212_CCEP'
        for ich = 1:length(chLabels)
            label = chLabels{ich};
            
            % get the non numerical portion
            label_num_idx = regexp(label,'\d');
            if isempty(label_num_idx), continue; end
            
            label_non_num = label(1:label_num_idx-1);
            
            switch label_non_num
                case 'LC'
                    ana{ich} = 'hippocampus';%'posterior\newlinehippocampus';
                case 'LB'
                    ana{ich} = 'hippocampus';%'anterior\newlinehippocampus';
                case 'LA'
                    ana{ich} = 'amygdala';
                case 'LD'
                    ana{ich} = 'anterior\newlineinsula';
                case 'LE'
                    ana{ich} = 'central\newlinesulcus';
                case 'LF'
                    ana{ich} = 'posterior\newlineinsula';
                case 'LG'
                    ana{ich} = 'pars\newlinetriangularis';
                case 'LH'
                    ana{ich} = 'inferior\newlinefrontal gyrus';
                case 'LI'
                    ana{ich} = 'anterior\newlinecingulate';
                case 'LJ'
                    ana{ich} = 'middle\newlinefrontal gyrus';
            end
            
        end 
    otherwise
        ana = [];
end
%}
end