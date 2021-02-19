function to_generate_DiLiPop_param_structure(name_matfile,path_save_param,nb_embryo)

% function to create the structure, containing the parameters needed for the
% DiLiPop analysis, to the condition you aim to study, filled with the
% proper values.

% general_param are parameters common to the whole embryo dataset
% param are parameters associated to each embryo

%% input args
% 1st arg: saving name of a xmlfile (mat file) for a given condition that contains the parameters associlated
% 2nd rag: saving path where the xmlfile (mat file) for a given condition that contains the parameters associlated
% 3rd arg: number of embryos associated to the studied dataset

if nargin < 1
    name_matfile = input_perso('Please give the name of the mat structure saving parameters of the given dataset (eg. PID4_TH65_none_DiLiPop_parameters.mat): ','s');
end

if nargin < 2
    path_save_param = uigetdir(pwd,'Please give a path where to save the mat structure saving parameters of the given dataset');
end

if nargin < 3
    nb_embryo = input('Please give the number of embryos in the given condition: ');
end

saveVarsMat = load('DiLiPop_model.mat'); % model of the structure used to create the structure of your condition
general_param = saveVarsMat.general_param;

%% edit general params here.

% wish to detect furrow ingression position along time and detect its onset: yes = 1, no = 0
furrow_detection_status = input_perso([' Do you wish to detect/time furrow ingression? (yes = 1, no = 0): '],1  );
general_param.furrow_detection.status_perform = furrow_detection_status;

% path to save...
general_param.work_path = '/Projects_Space/People/Helene/Working_dir';

% save new general_param
saveVarsMat_new.general_param = general_param;
clear general param

%% for each embryo

for k=1:nb_embryo
    
    param = saveVarsMat.params{1};
    
    %% edit param  here.
    
    if k == 1 % ask user only for the first embryo
        % below param that will be the same for all the embryos of the dataset
        format_image = input_perso([' Please, give the format of your images (mtif = tif stack): '],'mtif'  ); %param.format_image = 'mtif';
        param.format_image = format_image;
        frequency = input_perso([' Please, give the image acquisition rate (Hz): '],10  ); % image acquisition frequency in Hz
        param.sp6 = frequency;
        resolution = input_perso([' Please, give the spatial resolution (nm): '],139  );
        param.resol = resolution; % in nm
        decimate = input_perso([' Please, give image step between contour detection: '],10  );
        param.cortex_pass2.decimate = decimate; % embryo contour performed every 10 frames
        
        % path of your main data directory (directory for a given condition,
        % that will contain subfolder for each embryo to save informations for embryos of given condition
        basepath = uigetdir(pwd,'Please give a main path where to save data common to all embryos (param.basepath)');
        param.basepath = basepath;
        
        % duration between various events of the mitosis
        delta_furrowDetection_anaphaseOnset = input_perso([' Please, give the duration between anaphase onset and cytokinesis furrow detection (in s): '],101  ); % 101 true for TH65 at 23°C
        param.delta_furrowDetection_anaphaseOnset = delta_furrowDetection_anaphaseOnset;
        delta_furrowDetection_nebd = input_perso([' Please, give the duration between nebd and cytokinesis furrow detection (in s): '],295  );% 101 true for TH65 at 23°C
        param.delta_furrowDetection_nebd = delta_furrowDetection_nebd;
        NEBD_pseudoCleavage= input_perso([' Please, give the duration between pseudocleavage and nebd (in s): '],169  );% 101 true for TH65 at 23°C
        param.delta_NEBD_pseudoCleavage = NEBD_pseudoCleavage;
        delta_late_anaphase = input_perso([' Please, give the duration between anaphase onset and late anaphase (in s): '],200  );% 101 true for TH65 at 23°C
        param.delta_late_anaphase = delta_late_anaphase;
        delta_early_metaphase = input_perso([' Please, give the duration between nebd and anaphase onset (in s): '],193  );% 101 true for TH65 at 23°C
        param.delta_early_metaphase = delta_early_metaphase;
        
        % param below might be edited to user conditions
        param.sp7 = -1;
        param.channel_of_interest_AC = 1;
        param.mask = '';
        param.cortex_pass1.mask_image = '';
        param.cortex_pass2.mask_image = '';
        param.alpha = 0;
    else
        param.format_image = format_image;
        param.sp6 = frequency;
        param.resol = resolution;
        param.cortex_pass2.decimate = decimate;
        param.basepath = basepath;
        param.delta_furrowDetection_anaphaseOnset = delta_furrowDetection_anaphaseOnset;
        param.delta_furrowDetection_nebd = delta_furrowDetection_nebd;
        param.delta_NEBD_pseudoCleavage = NEBD_pseudoCleavage;
        param.delta_late_anaphase = delta_late_anaphase;
        param.delta_early_metaphase = delta_early_metaphase;
        % param below might be edited to user conditions
        param.sp7 = -1;
        param.channel_of_interest_AC = 1;
        param.mask = '';
        param.cortex_pass1.mask_image = '';
        param.cortex_pass2.mask_image = '';        
        param.alpha = 0;
    end
    
    %--------------------------------------------------------------------------------
    % ask questions to user to adjust the parameter to the studied embryo
    
    disp(['Embryo ' num2str(k)]);
    
    param.image_ref = input_perso([' Please, give the image reference index: '], 1  );
    param.stem_name = input_perso([' Please, give the name of given embryo (param.stem_name): (e.g. 2018-04-23_1) '], 's' );
    param.sp1 = input_perso([' Please, give the folder name associated to the studied embryo (param.sp1) (e.g. 2018-04-23_TH65_cortex_1): '], 's'  );
    param.sp2 = input_perso([' Please, give the first image index: '], 1  );
    param.sp3 = input_perso([' Please, give the last image index: '], 4000  );
    param.extra = input_perso([' Please, give a tag to embryo (empty = no tag): '], ''  );
    
    
    %--------------------------------------------------------------------------------------------------
%     % give the possibility to user to generate a ROI selection using polygonal selection tool: the mask will be saved
%     path_tiff_stack = strcat(param.basepath , '/' , param.sp1 , '/');
%     name_tiff_stack = param.stem_name;
%     name = fullfile(path_tiff_stack,name_tiff_stack);
%     
%     [end_reading,siz,fitsnom] = read_init(name,param.format_image,param.sp7,param.sp2,param.sp3,'');
%     [Image,~,error_reading] = read_with_trial(param.sp2,param.sp7,param.format_image,siz,fitsnom,'none',param.sp3,param.cortex_pass2.channel_interest_AC);
%     
%     % user will do a selection using roipoly tool
%     mask_BW_AC = roipoly(imadjust_perso(Image));
%     name_mask = strcat(param.stem_name,'_mask.tif');
%     name_path_mask = fullfile(path_tiff_stack,name_mask);
%     imwrite(mask_BW_AC,name_path_mask);
%     param.cortex_pass1.mask_image = name_path_mask;
%     param.cortex_pass2.mask_image = name_path_mask;
%     clear Image
    
    %----------------------
    % save the parameter in the new structure
    saveVarsMat_new.params{k} = param;
    clear param
    
end

% save the structure of the parameters of your dataset
save(fullfile(path_save_param,name_matfile), 'saveVarsMat_new');

clear saveVarsMat_new saveVarsMat param
clear all

end