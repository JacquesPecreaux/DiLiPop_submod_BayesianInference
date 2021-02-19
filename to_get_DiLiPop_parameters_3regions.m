function to_get_DiLiPop_parameters_3regions( xmlfile,save_stem,choiceModel,choiceCondition1,choiceCondition2,algo_choice, ...
    error_computation,time_reference_choice,folder_tag,limit_nb_tracks_for_fitting,minLength_tracks )
    

% This functions enables to do the DiLiPop statistical analysis, investigating the difference in the 3 cortical regions of the embryo
% From the durations of the tracks, it reveals the number of distinct dynamical behaviours and their parameters using
% Bayesian inference approach

% OUTPUT FILES of interest
% - mat file: tracks_duration_histo that contains the duration distributions of all the embryos for each period/region investigated
% - mat file: final_results-BayesianInference_sum_mle that contains the result of the fitting of the duration distributions 
% in each region/period couple
% - plots of the experimental duration distributions fitted by the different models
% - plots of the experimental duration distributions fitted with the best model
% - best model parameter values and errors saved in text file 'tracks_duration_histo_values-BayesianInference.txt' that can be opened by excel


%% input args

% 1st: xmlfile that contains information regarding :
% - the studied condition (general_param.),
% - each sample (e.g. embryo) composing this condition (param.)
% 2nd: path to save the results of the analysis
% 3rd: choice of the model that will be tested to adjust the experimental duration distributions 
% 4th: choice about the time period to investigate: eg. prophase, metaphase, anaphase, ...
% 5th: choice about the regions to investigate: whole embryo, region 1 /2/3.
% 6th: ask which algo to use to adjust the experimental duration distribution  with the model
% 7th: ask whether you wish to compute the errors associated to the model parameters
% 8th: reference time used to align the embryos : the two time reference possible are 1/ the onset of the furrow ingression 
% at about mid-anaphase, and 2/ the pseudo cleavage end happeening before pronuclei meeting
% 9th: tag of the folder: the data needed for this function are saved in the main folder "data_dir", and in the subfolder 
% corresponded to the project (PID), sub-subfolder correponddnt to the embyro (ID), and in the sub-sub-subfolder identified 
% with the given tag. This tag can be empty if no tag used.
% 10th: minimal number of tracks required to perform correctly the DiLiPop statistical analysis
% 11th: minimal duration of a track considered in the analysis


%% files that will be downloaded and required for the analysis

% These data  are saved in the main folder "data_dir", and in the subfolder corresponded to the project (PID), 
% sub-subfolder correpondant to the embyro (ID), and in the sub-sub-subfolder identified with the given tag.
% Files neded are:
% - region Area: area of the embryo along mitosis: whole, or divided alon AP axis
% - regionXlimit: coordinates of the embryo limits in the full frame along AP axis, 
% --> will be needed to assign the contact to the different regions of the embryo
% - regionXlength: width of the sub-regions of the embryo
% --> will be needed to assign the contact to the different regions of the embryo
% - furrow_position_convexity: gives the furrow ingression onset timing
% furrow_onset in sec = furrow_position_convexity.image_start_detection/param.sp6 -> in second from the start of the movie
% - dataTracks_rotated: gives the tracks: duration, coordinates of the position, first frame the track appears, 
% end frame of the tracks.


%%

clear global
global general_param;
global param;


%% process inputs and ask for them
if nargin<1
    [xmlfile, p] = uigetfile('*.mat','Please choose a job file to process');
    [~,xmlfile_bkp,~] = fileparts(xmlfile);
    xmlfile = fullfile(p, xmlfile);
end
if nargin<2
    [p,f] = fileparts(xmlfile);
    [save_stem, p] = uiputfile('*','Please provide a path and stem name for saving figures',fullfile(p,[f '_3regions_']));
else
    if ~exist(save_stem,'dir')
        [p,save_stem] = fileparts(save_stem);
        if isempty(p)
            p = pwd;
        end
    else
        p = save_stem;
        [~,save_stem] = fileparts(xmlfile);
    end
end
p0 = fullfile(p,[save_stem, date]);
p1 = unique_name(p0,'dir');
if ~exist(p1,'dir')
    mkdir(p1);
end
p = p1;
main_path = p;
save_stem = fullfile(p,save_stem);
if ~isempty(p)
    cd(p);
end
copyfile(xmlfile,p)
diary([save_stem '.log']);
diary on;


if nargin<3   
    choiceModel = input_perso(['Enter choice for model (a- MonoExpo, b- DoubleExpo, c- MonoExpo stretched, d- TripleExpo ): ' ...
        '\n (1) a+b, \n (2) a+b+c  \n (3) a+b+d  \n (4) a+b+c+d ' ], 1);   
    switch choiceModel
        case 1
            models = {'MonoExpo' 'DoubleExpo' };
        case 2
            models = {'MonoExpo' 'DoubleExpo' 'MonoExpo_stretched'};
        case 3
            models = {'MonoExpo' 'DoubleExpo' 'TripleExpo' };          
        case 4
            models = {'MonoExpo' 'DoubleExpo'  'MonoExpo_stretched' 'TripleExpo'};         
        otherwise
            disp('bad choice');
    end    
end

if nargin < 4   
    choiceCondition1 = input_perso(['Enter choice for temporal condition (a- entireRecording, b- prometaphase, c- anaphase -d late, e- prophase, f- metaphase ): ' ...
        '\n (1) a, \n (2) b  \n (3) c  \n (4) d, \n (5) a+b+c+d  \n (6) b+c \n (7) c+d  \n (8) b+c+d  \n (9) a+b+c  \n (10) a+e+f+c+d ' ...
        ' \n (11) a+e+f+c \n (12) e+f+c+d \n (13) a+f+c+d : '], 10);   
    switch choiceCondition1
        case 1
            conditions1 = {'entireRecording'};
        case 2
            conditions1 = {'metaphase'};
        case 3
            conditions1 = {'anaphase'};
        case 4
            conditions1 = {'late'};
        case 5
            conditions1 = {'entireRecording' 'prometaphase' 'anaphase' 'late'};
        case 6
            conditions1 = {'metaphase' 'anaphase'};
        case 7
            conditions1 = {'anaphase' 'late'};
        case 8
            conditions1 = {'metaphase' 'anaphase' 'late'};
        case 9
            conditions1 = {'entireRecording' 'metaphase' 'anaphase'};
        case 10
            conditions1 = {'entireRecording' 'prophase' 'metaphase' 'anaphase' 'late'};
        case 11
            conditions1 = {'entireRecording' 'prophase' 'metaphase' 'anaphase'};     
        case 12
            conditions1 = {'prophase' 'metaphase' 'anaphase' 'late'};    
        case 13
            conditions1 = {'entireRecording' 'metaphase' 'anaphase' 'late'};            
        otherwise
            disp('bad choice');
    end    
end

if nargin < 5    
    choiceCondition2 = input_perso(['Enter choice for spatial condition (a- entireEmbryo, b- regon1, c- region2, d- region3 ): ' ...
        '\n (1) a, \n (2) b  \n (3) c  \n (4) d  \n (5) b+c+d \n (6) a+b+c+d : '], 5);    
    switch choiceCondition2
        case 1
            conditions2 = {'entireEmbryo'};
        case 2
            conditions2 = {'region1'};
        case 3
            conditions2 = {'region2'};
        case 4
            conditions2 = {'region3'};
        case 5
            conditions2 = {'region1' 'region2' 'region3'};
        case 6
           conditions2 = {'entireEmbryo' 'region1' 'region2' 'region3'};           
        otherwise
            disp('bad choice');
    end
end

if nargin < 6
    algo_choice = input_perso(['Which algo of optimization?(fmincon =1 , patternsearch = 2) '], 1);
end

if nargin < 7
    error_computation = input_perso(['Do you wish to estimate model parameter errors? (yes = 1, no = 0) ?'], 1);
    if error_computation == 1
        choice_error_estimate = input_perso(['Which method do you wish to use to estimate the errors? (likelihood profile = 1, bootstrapping = 2, both = 3) '], 2);
    else
        choice_error_estimate = NaN;
    end
end

if nargin < 8
    time_reference_choice = input_perso(['Which reference time to use? (0= furrow ingression onset, 1 = pseudo-cleavage end '], 0);
end

if nargin < 9
    folder_tag = input_perso(['Set the tag used to get proper treatment at the cortex folder): '],'');
end


%% load parameters from the downloaded strcuture

load(xmlfile); % load structure containing general_param and param
general_param = saveVarsMat_new.general_param;% read general_param from structure

if nargin<10
    limit_nb_tracks_for_fitting = general_param.cortex_analysis.threshold_MTnumber_forMLEfit;
end

if nargin<11
minLength_tracks  = ...
    input_perso([sprintf('Please give minimal duration to be studied at cortex(default value = ')], general_param.cortex_analysis.minLength );
end

general_param.cortex_analysis.checking_tracks = 0;


%% treat individual embryo and get their data

nbEmbryo = 0;
nbEmbryo_prometaphase = 0;
index_embryo_prometaphase = [];
nbEmbryo_prophase = 0;
index_embryo_prophase = [];
nbEmbryo_metaphase = 0;
index_embryo_metaphase = [];
nbEmbryo_anaphase = 0;
index_embryo_anaphase = [];
nbEmbryo_late = 0;
index_embryo_late = [];

cell_name = {};
cell_name_prophase = {};
cell_name_prometaphase = {};
cell_name_metaphase = {};
cell_name_anaphase = {};
cell_name_late = {};

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp(param.sp1);
        
        late_fitting = 0;
        metaphase_fitting = 0;
        prometaphase_fitting = 0;
        prophase_fitting = 0;
        anaphase_fitting = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % DOWNLOAD PART
        
        %-------------------------
        % download Xlimit
        mainDirectory = 'contour detection';
        pathMainDirectory3 = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        clear mainDirectory
        nameData = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', short_name) '.mat']; % in squared pixels
        nameData_ = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat']; % in squared pixels
        % case that maskStack_rotated file saved to avoid to redo from_contour_to_mask, rotate_mask functions and others
        if exist(nameData_,'file') == 2
            nameData_ = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat'];
            regionArea = load(nameData_);  % in pixels**2
            nameData2_ = [sprintf('%s%s%s',pathMainDirectory3,'regionXlimit-', param.extra, short_name) '.mat'];
            regionXlimit = load(nameData2_);
            nameData3_ = [sprintf('%s%s%s',pathMainDirectory3,'regionXlength-', param.extra, short_name) '.mat'];
            regionXlength = load(nameData3_);
            clear nameData
            clear nameData2
            clear nameData3
        elseif exist(nameData,'file') == 2
            nameData = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', short_name) '.mat'];
            regionArea = load(nameData); % in pixels**2
            nameData2 = [sprintf('%s%s%s',pathMainDirectory3,'regionXlimit-', short_name) '.mat'];
            regionXlimit = load(nameData2);
            nameData3 = [sprintf('%s%s%s',pathMainDirectory3,'regionXlength-', short_name) '.mat'];
            regionXlength = load(nameData3);
            clear nameData
            clear nameData2
            clear nameData3
        else
            disp('Error: area region, length region and position region not found');
        end
        
        %-------------------------------
        % download furrow position
        mainDirectory = 'furrow_characterization';
        pathMainDirectory = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        filename = ['furrow_position_convexity-', short_name, '.mat'];
        filename2 = strcat('furrow_position_convexity-', short_name,  param.extra, '.mat');
        
        nameData = fullfile( pathMainDirectory, filename );
        nameData2 = fullfile( pathMainDirectory, filename2 );
        
        if exist(nameData,'file') == 2
            furrow_position = load(fullfile(pathMainDirectory,filename));
        elseif exist(nameData2,'file') == 2
            furrow_position = load(fullfile(pathMainDirectory,filename2));
        else
            disp('Error: furrow position file not found');
        end
        clear filename filename2 nameData nameData2 mainDirectory
        
        
        %------------------------------
        % download dataTracks_rotated
        mainDirectory = strcat('analysis at the cortex', folder_tag);
        pathMainDirectory = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        nameData0 = [sprintf('%s%s%s',pathMainDirectory,'dataTracks_rotated-', short_name) '.mat'];
        nameData1 = [sprintf('%s%s%s%s',pathMainDirectory,'dataTracks_rotated-', short_name) , param.extra,'.mat'];
        nameData3 = [sprintf('%s%s%s%s',pathMainDirectory,'dataTracks_rotated_part1-', short_name) , param.extra,'.mat'];
        if exist(nameData3,'file') == 2
            filename1 = ['dataTracks_rotated_part1-', short_name, param.extra, '.mat'];
            dataTracks_rotated_part1 = load(fullfile(pathMainDirectory,filename1));
            filename2 = ['dataTracks_rotated_part2-', short_name, param.extra, '.mat'];
            dataTracks_rotated_part2 = load(fullfile(pathMainDirectory,filename2));
            dataTracks_rotated.entireEmbryo.tracksX = dataTracks_rotated_part1.dataTracks_rotated_part1.entireEmbryo.tracksX;
            dataTracks_rotated.entireEmbryo.indexXStart = dataTracks_rotated_part1.dataTracks_rotated_part1.entireEmbryo.indexXStart;
            dataTracks_rotated.entireEmbryo.lengthTracks = dataTracks_rotated_part1.dataTracks_rotated_part1.entireEmbryo.lengthTracks;
            dataTracks_rotated.entireEmbryo.tracksY = dataTracks_rotated_part2.dataTracks_rotated_part2.entireEmbryo.tracksY;
            dataTracks_rotated.entireEmbryo.indexXEnd = dataTracks_rotated_part2.dataTracks_rotated_part2.entireEmbryo.indexXEnd;
            dataTracks_rotated.entireEmbryo.numTracks = dataTracks_rotated_part2.dataTracks_rotated_part2.entireEmbryo.numTracks;
            dataTracks_rotated.numTimePoints = dataTracks_rotated_part2.dataTracks_rotated_part2.numTimePoints;
            clear dataTracks_rotated_part1 dataTracks_rotated_part2
        elseif exist(nameData1,'file') == 2
            filename1 = ['dataTracks_rotated-', short_name, param.extra, '.mat'];
            dataTracks_rotated = load(fullfile(pathMainDirectory,filename1));
            FileInfo = dir(fullfile(pathMainDirectory,filename1));
            vector_datum = datevec(FileInfo.datenum);
            clear filename1
        elseif exist(nameData0,'file') == 2
            filename0 = ['dataTracks_rotated-', short_name, '.mat'];
            dataTracks_rotated = load(fullfile(pathMainDirectory,filename0));
            FileInfo = dir(fullfile(pathMainDirectory,filename0));
            vector_datum = datevec(FileInfo.datenum);
            clear filename0
        else
            disp('Error: dataTracks_rotated file not found');
        end
        
        
        %-----------------------------------
        % GET DATATRACKS STRUCTURE
        
        [ dataTracks_new,area_3regions,xStart_3regions,xEnd_3regions ] = to_get_dataTracks_for_global_residencyTime_analysis_3...
            ( dataTracks_rotated,furrow_position,regionXlength,regionXlimit,regionArea,conditions1,main_path,[],time_reference_choice );  % in pixels**2
        
        name = strcat('dataTracks_ID' , param.sp1, '.mat');
        save([save_stem name],  '-struct','dataTracks_new');
        
        nbEmbryo = nbEmbryo +1;
        cell_name{nbEmbryo} = short_name;
        
        
        %---------------------------------
        % identify which embryo can be studied in the different periods according to their number of tracks
        
        if ~isempty(find(contains('prophase',conditions1)))
            number_tracks_prophase1 = dataTracks_new.prophase.region1.numTracks;
            number_tracks_prophase2 = dataTracks_new.prophase.region2.numTracks;
            number_tracks_prophase3 = dataTracks_new.prophase.region3.numTracks;
            if ( number_tracks_prophase1 >= limit_nb_tracks_for_fitting )|| ( number_tracks_prophase2 >= limit_nb_tracks_for_fitting ) ||...
                    ( number_tracks_prophase3 >= limit_nb_tracks_for_fitting )
                nbEmbryo_prophase = nbEmbryo_prophase +1;
                index_embryo_prophase = [ [index_embryo_prophase] nbEmbryo ];
                cell_name_prophase{nbEmbryo_prophase} = short_name;
                prophase_fitting = 1;
            end
        end
        
        if ~isempty(find(contains('metaphase',conditions1))) && ~isempty(find(contains('prophase',conditions1)))
            number_tracks_metaphase1 = dataTracks_new.metaphase.region1.numTracks;
            number_tracks_metaphase2 = dataTracks_new.metaphase.region2.numTracks;
            number_tracks_metaphase3 = dataTracks_new.metaphase.region3.numTracks;
            if ( number_tracks_metaphase1 >= limit_nb_tracks_for_fitting )|| ( number_tracks_metaphase2 >= limit_nb_tracks_for_fitting ) ||...
                    ( number_tracks_metaphase3 >= limit_nb_tracks_for_fitting )
                nbEmbryo_metaphase = nbEmbryo_metaphase +1;
                index_embryo_metaphase = [ [index_embryo_metaphase] nbEmbryo ];
                cell_name_metaphase{nbEmbryo_metaphase} = short_name;
                metaphase_fitting = 1;
            end
        end
        
        if ~isempty(find(contains('anaphase',conditions1)))
            number_tracks_anaphase1 = dataTracks_new.anaphase.region1.numTracks;
            number_tracks_anaphase2 = dataTracks_new.anaphase.region2.numTracks;
            number_tracks_anaphase3 = dataTracks_new.anaphase.region3.numTracks;
            if ( number_tracks_anaphase1 >= limit_nb_tracks_for_fitting )|| ( number_tracks_anaphase2 >= limit_nb_tracks_for_fitting ) ||...
                    ( number_tracks_anaphase3 >= limit_nb_tracks_for_fitting )
                nbEmbryo_anaphase = nbEmbryo_anaphase +1;
                index_embryo_anaphase = [ [index_embryo_anaphase] nbEmbryo ];
                cell_name_anaphase{nbEmbryo_anaphase} = short_name;
                anaphase_fitting = 1;
            end
        end
        
        if ~isempty(find(contains('prometaphase',conditions1)))
            number_tracks_prometaphase1 = dataTracks_new.prometaphase.region1.numTracks;
            number_tracks_prometaphase2 = dataTracks_new.prometaphase.region2.numTracks;
            number_tracks_prometaphase3 = dataTracks_new.prometaphase.region3.numTracks;
            if ( number_tracks_prometaphase1 >= limit_nb_tracks_for_fitting )|| ( number_tracks_prometaphase2 >= limit_nb_tracks_for_fitting ) ||...
                    ( number_tracks_prometaphase3 >= limit_nb_tracks_for_fitting )
                nbEmbryo_prometaphase = nbEmbryo_prometaphase +1;
                index_embryo_prometaphase = [ [index_embryo_prometaphase] nbEmbryo ];
                cell_name_prometaphase{nbEmbryo_prometaphase} = short_name;
                prometaphase_fitting = 1;
            end
        end
        
        if ~isempty(find(contains('late',conditions1)))
            number_tracks_late1 = dataTracks_new.late.region1.numTracks;
            number_tracks_late2 = dataTracks_new.late.region2.numTracks;
            number_tracks_late3 = dataTracks_new.late.region3.numTracks;
            if ( number_tracks_late1 >= limit_nb_tracks_for_fitting )|| ( number_tracks_late2 >= limit_nb_tracks_for_fitting ) ||...
                    ( number_tracks_late3 >= limit_nb_tracks_for_fitting )
                nbEmbryo_late = nbEmbryo_late +1;
                index_embryo_late = [ [index_embryo_late] nbEmbryo ];
                cell_name_late{nbEmbryo_late} = short_name;
                late_fitting = 1;
            end
        end
        
        %---------------------------------
        % get experimental duration distributions for each embryo and in given period/region
        
        if ( nbEmbryo == 1 )
            [ tracks_duration_histo ] = to_get_tracks_duration_histo_finalAnalysis_3( dataTracks_new,nbEmbryo,...
                minLength_tracks,[],late_fitting,anaphase_fitting,metaphase_fitting,furrow_position,area_3regions,...
                conditions1,conditions2,main_path,prophase_fitting,prometaphase_fitting);
        else
            [ tracks_duration_histo ] = to_get_tracks_duration_histo_finalAnalysis_3( dataTracks_new,nbEmbryo,...
                minLength_tracks,tracks_duration_histo,late_fitting,anaphase_fitting,metaphase_fitting,furrow_position,area_3regions,...
                conditions1,conditions2,main_path,prophase_fitting,prometaphase_fitting);
        end
    end
    clear dataTracks
end

%% save structure with histo

name = strcat('tracks_duration_histo.mat');
save(fullfile(main_path,name), 'tracks_duration_histo');


%% fitting for each condition

% identify the conditions for which at least two embryos can be studied
if nbEmbryo_late < 2
    conditions1 = conditions1(~strcmp(conditions1,'late'));
end
if nbEmbryo_metaphase < 2
    conditions1 = conditions1(~strcmp(conditions1,'metaphase'));
end
if nbEmbryo_prometaphase < 2
    conditions1 = conditions1(~strcmp(conditions1,'prometaphase'));
end
if nbEmbryo_prophase < 2
    conditions1 = conditions1(~strcmp(conditions1,'prophase'));
end
if nbEmbryo_anaphase < 2
    conditions1 = conditions1(~strcmp(conditions1,'anaphase'));
end

n_condi1 = numel(conditions1);
n_condi2 = numel(conditions2);

first_condition_studied = 0;

% perform the fit of duration histogrem for each set of period/region
for iCondition=1:n_condi1
    
    name1 = conditions1{iCondition};
    
    if strcmp(name1,'metaphase')
        nbEmbryo_givenCondition = nbEmbryo_metaphase;
    elseif strcmp(name1,'late')
        nbEmbryo_givenCondition = nbEmbryo_late;
    elseif strcmp(name1,'prophase')
        nbEmbryo_givenCondition = nbEmbryo_prophase;
     elseif strcmp(name1,'prometaphase')
        nbEmbryo_givenCondition = nbEmbryo_prometaphase;   
     elseif strcmp(name1,'anaphase')
        nbEmbryo_givenCondition = nbEmbryo_anaphase;          
    else
        nbEmbryo_givenCondition = nbEmbryo;
    end
    
    
    for jCondition=1:n_condi2
        
        if nbEmbryo_givenCondition > 1
            
            name2 = conditions2{jCondition};
            
            first_condition_studied = first_condition_studied + 1;
            
            [fitting_results,model_choice] = tracks_duration_fitting_poisson_mle_all3...
                (models,name1,name2,nbEmbryo_givenCondition,tracks_duration_histo,index_embryo_metaphase,save_stem,...
                [],index_embryo_late,0,algo_choice,error_computation,choice_error_estimate,...
                index_embryo_prophase,index_embryo_prometaphase,index_embryo_anaphase);
            
            final_fitting_all.(name1).(name2).fitting = fitting_results;
            final_fitting_all.(name1).(name2).model = model_choice;
            
            clear fitting_results
            close all
            clear model_choice
            
        end
        
    end
end

%% save structure with results incorporated

name = strcat('final_results-BayesianInference_sum_mle.mat');
save(fullfile(main_path,name),'-struct', 'final_fitting_all');


%% save results in text file that could be opened with excel

save_fitting_results_inTextFile_mle_all( final_fitting_all,models,save_stem,1,nbEmbryo_metaphase,nbEmbryo,...
    [],minLength_tracks,nbEmbryo_late,conditions1,conditions2,folder_tag,0,[],[],algo_choice,error_computation,...
    choice_error_estimate,nbEmbryo_prophase,nbEmbryo_prometaphase,nbEmbryo_anaphase);


%% clear

clear all
close all


end

