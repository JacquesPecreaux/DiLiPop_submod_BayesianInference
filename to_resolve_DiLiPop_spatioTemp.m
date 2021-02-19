function to_resolve_DiLiPop_spatioTemp( xmlfile,save_stem,choiceModel,compute_wholeEmbryo,compute_3regions,error_compute,...
    time_reference_choice,additionnal_ref_nebd,folder_tag,limit_nb_tracks_for_fitting,minLength_tracks )

% This functions enables to do the DiLiPop statistical analysis, with a maximal resolution in time and space of the parameters

% OUTPUT FILES of interest
% - mat file: tracks_duration_histo that contains the duration distributions of all the embryos for each period/region investigated
% - mat file: final_results-BayesianInference_sum_mle that contains the result of the fitting of the duration distributions 
% in each region/period couple
% - plots of the experimental duration distributions fitted with the best for each period/region investigated
% - plots of the temporal evolutions in whole embro or in each region for the lifetimes, densities, frequencies and numbers. 
% - final structure of the temporal evolution of the model parameters: "results_temporalspatial_modulation"


%% input args

% 1st: xmlfile that contains information regarding :
% - the studied condition (general_param.),
% - each sample (e.g. embryo) composing this condition (param.)
% 2nd: path to save the results of the analysis
% 3rd: choice of the model that will be tested to adjust the experimental duration distributions 
% 4th: choice to resolve the whole embryo
% 5th: choice to resolve the 3 regions
% 6th: ask whether you wish to compute the errors associated to the model parameters
% 7th: reference time used to align the embryos : the two time reference possible are 1/ the onset of the furrow ingression 
% at about mid-anaphase, and 2/ the pseudo cleavage end happeening before pronuclei meeting
% 8th: if reference at cytokinesis furrow onset, ask for alignment to NEBD too
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
    [save_stem, p] = uiputfile('*','Please provide a path and stem name for saving figures',fullfile(p,[f '_spatioTemp_']));
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
    choiceModel = input_perso(['Enter choice for model (a- MonoExpo, b- DoubleExpo, c- MonoExpo stretched d- TripleExpo, ): ' ...
        '\n (0) b,\n (1) a+b, \n (2) a+b+c \n (3) a+b+d ' ], 1);    
    switch choiceModel
        case 0
            models = {'DoubleExpo' };
        case 1
            models = {'MonoExpo' 'DoubleExpo' };
        case 2
            models = {'MonoExpo' 'DoubleExpo' 'MonoExpo_stretched'};
        case 3
            models = {'MonoExpo' 'DoubleExpo' 'TripleExpo' };
        otherwise
            disp('bad choice');
    end    
end

if nargin<4
    compute_wholeEmbryo  = input_perso([sprintf('Do you wish to cmpute parameters in whole embryo? (yes = 1, No = 0 ')], 1 );
    if compute_wholeEmbryo == 1
        window_size_before_sec  = ...
            input_perso([sprintf('Please give window size during pro+metaphase/prophase in sec for whole embryo (default value = ')], 30 );
        
        window_size_after_sec = ...
            input_perso([sprintf('Please give window size during anaphase/meta+anaphase in sec  for whole embryo (default value = ')], 30 );
    end
end

if nargin<5
    compute_3regions  = input_perso([sprintf('Do you wish to compute parameters in 3 regions? (yes = 1, No = 0 ')], 1 );
    if compute_3regions == 1
        window_size_before_sec_3regions  = ...
            input_perso([sprintf('Please give window size during pro+metaphase/prophase in sec for 3 regions (default value = ')], 60 );
        
        window_size_after_sec_3regions = ...
            input_perso([sprintf('Please give window size during anaphase/meta+anaphase in sec  for 3 regions (default value = ')], 60 );
    end
end

if nargin<6
    error_compute = input_perso([sprintf('Do you wish to compute model parameter error (default value = ')], 0 );
end

if nargin < 7
    time_reference_choice = input_perso(['Which reference time to use? (0= furrow ingression onset, 1 = pseudo-cleavage end'], 0);
    if time_reference_choice == 0
        additionnal_ref_nebd = input_perso(['Do you wish to align data to NEBD? (0= no, 1 = yes) '], 0);
    else
        additionnal_ref_nebd = 0;
    end
end

if nargin < 9
    folder_tag = input_perso(['Set the tag used to get proper treatment at the cortex folder: '],'');
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


%% treat individual embryo and get their data
    
nbEmbryo = 0;
nbEmbryo_3regions = 0;
max_period_before = 0;
max_period_after = 0;
max_period_before_3regions = 0;
max_period_after_3regions = 0;
nbEmbryo_period_before = zeros(1,50);
nbEmbryo_period_after = zeros(1,50);
nbEmbryo_period_before_above250 = zeros(1,50);
nbEmbryo_period_after_above250 = zeros(1,50);
nbEmbryo_period_before_3regions = zeros(1,50);
nbEmbryo_period_after_3regions = zeros(1,50);
nbEmbryo_period_before_above250_3regions.region1 = zeros(1,50);
nbEmbryo_period_after_above250_3regions.region1 = zeros(1,50);
nbEmbryo_period_before_above250_3regions.region2 = zeros(1,50);
nbEmbryo_period_after_above250_3regions.region2 = zeros(1,50);
nbEmbryo_period_before_above250_3regions.region3 = zeros(1,50);
nbEmbryo_period_after_above250_3regions.region3 = zeros(1,50);

for i = 1 : 50
    if compute_wholeEmbryo == 1
        if time_reference_choice == 0
            extension1 = ['before_minus' num2str(i)];
            embryo_index_before.(extension1).entireEmbryo = [];
            extension2 = ['after_plus' num2str(i)];
            embryo_index_after.(extension2).entireEmbryo = [];
        elseif time_reference_choice == 1
            extension1 = ['NEBD_minus' num2str(i)];
            embryo_index_before.(extension1).entireEmbryo = [];
            extension2 = ['NEBD_plus' num2str(i)];
            embryo_index_after.(extension2).entireEmbryo = [];
        end
    end
    if compute_3regions == 1
        if time_reference_choice == 0
            extension1 = ['before_minus' num2str(i) '_3regions'];
            embryo_index_before.(extension1).region1 = [];
            embryo_index_before.(extension1).region2 = [];
            embryo_index_before.(extension1).region3 = [];
            extension2 = ['after_plus' num2str(i) '_3regions'];
            embryo_index_after.(extension2).region1 = [];
            embryo_index_after.(extension2).region2 = [];
            embryo_index_after.(extension2).region3 = [];
        elseif time_reference_choice == 1
            extension1 = ['NEBD_minus' num2str(i) '_3regions'];
            embryo_index_before.(extension1).region1 = [];
            embryo_index_before.(extension1).region2 = [];
            embryo_index_before.(extension1).region3 = [];
            extension2 = ['NEBD_plus' num2str(i) '_3regions'];
            embryo_index_after.(extension2).region1 = [];
            embryo_index_after.(extension2).region2 = [];
            embryo_index_after.(extension2).region3 = [];
        end
    end
end

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp(param.sp1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DOWNLOAD PART
        
        %-------------------------
        % download Xlimit, Xlength, Area
        mainDirectory = 'contour detection';
        pathMainDirectory3 = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        clear mainDirectory
        nameData = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', short_name) '.mat']; % in squared pixels
        nameData2 = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat']; % in squared pixels
        % case that maskStack_rotated file saved to avoid to redo from_contour_to_mask, rotate_mask functions and others
        if exist(nameData2,'file') == 2
            nameData = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat'];
            regionArea = load(nameData); % in pixels**2
            nameData2 = [sprintf('%s%s%s',pathMainDirectory3,'regionXlimit-', param.extra, short_name) '.mat'];
            regionXlimit = load(nameData2);
            nameData3 = [sprintf('%s%s%s',pathMainDirectory3,'regionXlength-', param.extra, short_name) '.mat'];
            regionXlength = load(nameData3);
            clear nameData
            clear nameData2
            clear nameData3
        elseif  exist(nameData,'file') == 2
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
            disp('Error: area region, length region and position region files not found');
        end
        
        
        %-------------------------------
        % download furrow position
        
        mainDirectory = 'furrow_characterization';
        pathMainDirectory = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        filename = ['furrow_position_convexity-', short_name, '.mat'];
        filename2 = strcat('furrow_position_convexity-', short_name,  param.extra, '.mat');
        nameData = fullfile( pathMainDirectory, filename );
        nameData2 = fullfile( pathMainDirectory, filename2 );
        if exist(nameData2,'file') == 2
            furrow_position = load(fullfile(pathMainDirectory,filename2));
        elseif exist(nameData,'file') == 2
            furrow_position = load(fullfile(pathMainDirectory,filename));
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
        elseif exist(nameData0,'file') == 2
            filename0 = ['dataTracks_rotated-', short_name, '.mat'];
            dataTracks_rotated = load(fullfile(pathMainDirectory,filename0));
        end
        % if saving not done as expected
        if  isfield(dataTracks_rotated,'dataTracks_rotated') == 1
            dataTracks_rotated_.entireEmbryo = dataTracks_rotated.dataTracks_rotated.entireEmbryo;
            dataTracks_rotated_.numTimePoints = dataTracks_rotated.dataTracks_rotated.numTimePoints;
            clear dataTracks_rotated
            dataTracks_rotated = dataTracks_rotated_;
            clear dataTracks_rotated_
        end
        %if dataTracks has been saved as dataTracks.entireEmbryo.entireRecording...
        if  isfield(dataTracks_rotated.entireEmbryo,'entireRecording') == 1
            dataTracks_rotated_.entireEmbryo = dataTracks_rotated.entireEmbryo.entireRecording;
            dataTracks_rotated_.numTimePoints = dataTracks_rotated.numTimePoints;
            clear dataTracks_rotated
            dataTracks_rotated = dataTracks_rotated_;
            clear dataTracks_rotated_
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get data split into 3 regions
        
        if compute_3regions == 1
            
            nbEmbryo_3regions = nbEmbryo_3regions + 1;
            %------------------------
            % 3 region area, limits
            
            limit1 = general_param.cortex_analysis.analysis_3SpatialRegions_limit1;
            limit2 = general_param.cortex_analysis.analysis_3SpatialRegions_limit2;
            [ area_3regions,xStart_3regions,xEnd_3regions ] = get_infos_3regions( regionXlimit,regionXlength,regionArea,limit1,limit2 ); % area in squared pixels
            
            % ------------------------
            % dataTracks split tinto 3 regions
            
            numTracks_region1 = 0;
            numTracks_region2 = 0;
            numTracks_region3 = 0;
            numTracks_outside = 0;
            
            sp2_bkp = param.sp2;
            for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
                if dataTracks_rotated.entireEmbryo.indexXStart(1,i) - param.sp2 <=0
                    param.sp2 = 0;
                end
                if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region1 = numTracks_region1 +1;
                    lengthTracks_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region1(:,numTracks_region1) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region1(:,numTracks_region1) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region2 = numTracks_region2 +1;
                    lengthTracks_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region2(:,numTracks_region2) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region2(:,numTracks_region2) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region3 = numTracks_region3 +1;
                    lengthTracks_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region3(:,numTracks_region3) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region3(:,numTracks_region3) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                else
                    numTracks_outside = numTracks_outside +1;
                    
                end
            end
            
            %save data from tracking in structure file
            dataTracks_rotated.region1.tracksX = tracksX_region1;
            dataTracks_rotated.region1.tracksY = tracksY_region1;
            dataTracks_rotated.region1.lengthTracks = lengthTracks_region1;
            dataTracks_rotated.region1.indexXStart = indexXStart_region1;
            dataTracks_rotated.region1.indexXEnd = indexXEnd_region1;
            dataTracks_rotated.region1.numTracks = numTracks_region1;
            
            dataTracks_rotated.region2.tracksX = tracksX_region2;
            dataTracks_rotated.region2.tracksY = tracksY_region2;
            dataTracks_rotated.region2.lengthTracks = lengthTracks_region2;
            dataTracks_rotated.region2.indexXStart = indexXStart_region2;
            dataTracks_rotated.region2.indexXEnd = indexXEnd_region2;
            dataTracks_rotated.region2.numTracks = numTracks_region2;
            
            dataTracks_rotated.region3.tracksX = tracksX_region3;
            dataTracks_rotated.region3.tracksY = tracksY_region3;
            dataTracks_rotated.region3.lengthTracks = lengthTracks_region3;
            dataTracks_rotated.region3.indexXStart = indexXStart_region3;
            dataTracks_rotated.region3.indexXEnd = indexXEnd_region3;
            dataTracks_rotated.region3.numTracks = numTracks_region3;
            
            dataTracks_rotated.numTracks_outside = numTracks_outside;
            
            clear tracksX_region3 tracksY_region3 lengthTracks_region3 indexXStart_region3 indexXEnd_region3
            clear tracksX_region2 tracksY_region2 lengthTracks_region2 indexXStart_region2 indexXEnd_region2
            clear tracksX_region1 tracksY_region1 lengthTracks_region1 indexXStart_region1 indexXEnd_region1
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform temporal separation of dataTracks
        
        % get temporal limits of the windows
        if time_reference_choice == 0
            image_reference = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6;
        elseif time_reference_choice == 1
            image_reference = ( param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage ) * param.sp6;
        end
        
        nbEmbryo = nbEmbryo + 1;
        
        % whole embryo windowing
        if compute_wholeEmbryo == 1
            window_size_before = window_size_before_sec * param.sp6; % in image nb
            window_size_after = window_size_after_sec * param.sp6; % in image nb
            
            % get reference time in metaphase or prophase
            if image_reference >= param.sp2
                image_reference_before = [ image_reference ];
                for ii = 1 : 50
                    tmp = image_reference - ii * window_size_before;
                    if tmp >  param.sp2
                        image_reference_before = [ [image_reference_before ] tmp ];
                    else
                        break
                    end
                end
                image_reference_before = [ [image_reference_before ] param.sp2 ];
            else
                image_reference_before = [];
                max_period_before = 0;
            end
            
            for jj = 1 : length(image_reference_before)
                nbEmbryo_period_before(1,jj) = nbEmbryo_period_before(1,jj) + 1;
            end
            if max_period_before < length(image_reference_before)
                max_period_before = length(image_reference_before);
            end
            
            % get reference time in anaphase
            if image_reference > param.sp2
                image_reference_after = [image_reference];
            else
                % image_reference_after = [param.sp2];
                image_reference_after = [];
            end
            for ii = 1 : 50
                tmp = image_reference + ii * window_size_after;
                if tmp > param.sp2
                    if tmp <  param.sp3
                        image_reference_after = [ [image_reference_after ] tmp ];
                    else
                        break
                    end
                end
            end
            image_reference_after = [ [image_reference_after ] param.sp3 ];
            
            for jj = 1 : length(image_reference_after)
                nbEmbryo_period_after(1,jj) = nbEmbryo_period_after(1,jj) + 1;
            end
            
            if max_period_after < length(image_reference_after)
                max_period_after = length(image_reference_after);
            end
        end
        
        % regions
        if compute_3regions == 1
            window_size_before_3regions = window_size_before_sec_3regions * param.sp6; % in image nb
            window_size_after_3regions = window_size_after_sec_3regions * param.sp6; % in image nb
            
            % get reference time in metaphase or prophase
            if image_reference >= param.sp2
                image_reference_before_3regions = [ image_reference ];
                for ii = 1 : 50
                    tmp = image_reference - ii * window_size_before_3regions;
                    if tmp >  param.sp2
                        image_reference_before_3regions = [ [image_reference_before_3regions ] tmp ];
                    else
                        break
                    end
                end
                image_reference_before_3regions = [ [image_reference_before_3regions ] param.sp2 ];
            else
                image_reference_before_3regions = [];
                max_period_before_3regions = 0;
            end
            
            for jj = 1 : length(image_reference_before_3regions)
                nbEmbryo_period_before_3regions(1,jj) = nbEmbryo_period_before_3regions(1,jj) + 1;
            end
            if max_period_before_3regions < length(image_reference_before_3regions)
                max_period_before_3regions = length(image_reference_before_3regions);
            end
            
            % get reference time in anaphase
            if image_reference > param.sp2
                image_reference_after_3regions = [image_reference];
            else
                image_reference_after_3regions = [param.sp2];
            end
            for ii = 1 : 50
                tmp = image_reference + ii * window_size_after_3regions;
                if tmp > param.sp2
                    if tmp <  param.sp3
                        image_reference_after_3regions = [ [image_reference_after_3regions ] tmp ];
                    else
                        break
                    end
                end
            end
            image_reference_after_3regions = [ [image_reference_after_3regions ] param.sp3 ];
            
            for jj = 1 : length(image_reference_after_3regions)
                nbEmbryo_period_after_3regions(1,jj) = nbEmbryo_period_after_3regions(1,jj) + 1;
            end
            
            if max_period_after_3regions < length(image_reference_after_3regions)
                max_period_after_3regions = length(image_reference_after_3regions);
            end
            
        end
        
        %--------------
        % separate tracks in temporal windows
        
        % whole embryo
        if compute_wholeEmbryo == 1
            
            for iPeriod_before = 1 : length(image_reference_before)-1
                if time_reference_choice == 0
                    extension = ['before_minus' num2str(iPeriod_before)];
                elseif time_reference_choice == 1
                    extension = ['NEBD_minus' num2str(iPeriod_before)];
                end
                entireEmbryo.(extension).numTracks = 0;
                entireEmbryo.(extension).lengthTracks = [];
            end
            
            for iPeriod_after = 1 : length(image_reference_after)-1
                if time_reference_choice == 0
                    extension = ['after_plus' num2str(iPeriod_after)];
                elseif time_reference_choice == 1
                    extension = ['NEBD_plus' num2str(iPeriod_after)];
                end
                entireEmbryo.(extension).numTracks = 0;
                entireEmbryo.(extension).lengthTracks = [];
            end
            
            for iTrack = 1 : dataTracks_rotated.entireEmbryo.numTracks
                for iPeriod_before = 1 : length(image_reference_before)-1
                    if time_reference_choice == 0
                        extension = ['before_minus' num2str(iPeriod_before)];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_minus' num2str(iPeriod_before)];
                    end
                    
                    if image_reference_before(iPeriod_before+1) <= dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) < image_reference_before(iPeriod_before)
                        
                        entireEmbryo.(extension).numTracks = entireEmbryo.(extension).numTracks +1;
                        entireEmbryo.(extension).lengthTracks(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.lengthTracks(iTrack);
                    end
                end
                for iPeriod_after = 1 : length(image_reference_after)-1
                    if time_reference_choice == 0
                        extension = ['after_plus' num2str(iPeriod_after)];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_plus' num2str(iPeriod_after)];
                    end
                    
                    if image_reference_after(iPeriod_after) <= dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) < image_reference_after(iPeriod_after+1)
                        
                        entireEmbryo.(extension).numTracks = entireEmbryo.(extension).numTracks +1;
                        entireEmbryo.(extension).lengthTracks(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.lengthTracks(iTrack);
                    end
                end
            end
        end
        
        % for 3 regions
        if compute_3regions == 1
            
            for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                if time_reference_choice == 0
                    extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
                elseif time_reference_choice == 1
                    extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
                end
                region1.(extension).numTracks = 0;
                region1.(extension).lengthTracks = [];
                region2.(extension).numTracks = 0;
                region2.(extension).lengthTracks = [];
                region3.(extension).numTracks = 0;
                region3.(extension).lengthTracks = [];
            end
            
            for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                if time_reference_choice == 0
                    extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
                elseif time_reference_choice == 1
                    extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
                end
                region1.(extension).numTracks = 0;
                region1.(extension).lengthTracks = [];
                region2.(extension).numTracks = 0;
                region2.(extension).lengthTracks = [];
                region3.(extension).numTracks = 0;
                region3.(extension).lengthTracks = [];
            end
            
            for iTrack = 1 : dataTracks_rotated.region1.numTracks
                for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                    if time_reference_choice == 0
                        extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
                    end
                    if image_reference_before_3regions(iPeriod_before+1) <= dataTracks_rotated.region1.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region1.indexXStart(1,iTrack) < image_reference_before_3regions(iPeriod_before)
                        region1.(extension).numTracks = region1.(extension).numTracks +1;
                        region1.(extension).lengthTracks(region1.(extension).numTracks) = dataTracks_rotated.region1.lengthTracks(iTrack);
                    end
                end
                for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                    if time_reference_choice == 0
                        extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
                    end
                    if image_reference_after_3regions(iPeriod_after) <= dataTracks_rotated.region1.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region1.indexXStart(1,iTrack) < image_reference_after_3regions(iPeriod_after+1)
                        region1.(extension).numTracks = region1.(extension).numTracks +1;
                        region1.(extension).lengthTracks(region1.(extension).numTracks) = dataTracks_rotated.region1.lengthTracks(iTrack);
                    end
                end
            end
            
            for iTrack = 1 : dataTracks_rotated.region2.numTracks
                for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                    if time_reference_choice == 0
                        extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
                    end
                    if image_reference_before_3regions(iPeriod_before+1) <= dataTracks_rotated.region2.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region2.indexXStart(1,iTrack) < image_reference_before_3regions(iPeriod_before)
                        region2.(extension).numTracks = region2.(extension).numTracks +1;
                        region2.(extension).lengthTracks(region2.(extension).numTracks) = dataTracks_rotated.region2.lengthTracks(iTrack);
                    end
                end
                for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                    if time_reference_choice == 0
                        extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
                    end
                    if image_reference_after_3regions(iPeriod_after) <= dataTracks_rotated.region2.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region2.indexXStart(1,iTrack) < image_reference_after_3regions(iPeriod_after+1)
                        region2.(extension).numTracks = region2.(extension).numTracks +1;
                        region2.(extension).lengthTracks(region2.(extension).numTracks) = dataTracks_rotated.region2.lengthTracks(iTrack);
                    end
                end
            end
            
            for iTrack = 1 : dataTracks_rotated.region3.numTracks
                for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                    if time_reference_choice == 0
                        extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
                    end
                    if image_reference_before_3regions(iPeriod_before+1) <= dataTracks_rotated.region3.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region3.indexXStart(1,iTrack) < image_reference_before_3regions(iPeriod_before)
                        region3.(extension).numTracks = region3.(extension).numTracks +1;
                        region3.(extension).lengthTracks(region3.(extension).numTracks) = dataTracks_rotated.region3.lengthTracks(iTrack);
                    end
                end
                for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                    if time_reference_choice == 0
                        extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
                    elseif time_reference_choice == 1
                        extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
                    end
                    if image_reference_after_3regions(iPeriod_after) <= dataTracks_rotated.region3.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region3.indexXStart(1,iTrack) < image_reference_after_3regions(iPeriod_after+1)
                        region3.(extension).numTracks = region3.(extension).numTracks +1;
                        region3.(extension).lengthTracks(region3.(extension).numTracks) = dataTracks_rotated.region3.lengthTracks(iTrack);
                    end
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get areas in each temporal window and regions
        
        % whole embryo
        if compute_wholeEmbryo == 1
            for iPeriod_before = 1 : length(image_reference_before)-1
                regionArea_period_before(iPeriod_before) = ...
                    mean( regionArea.entireEmbryo.nbR1(image_reference_before(iPeriod_before+1)-param.sp2+1: image_reference_before(iPeriod_before)-param.sp2+1) ).* ( (param.resol/1000)^2 );
            end
            
            for iPeriod_after = 1 : length(image_reference_after)-1
                regionArea_period_after(iPeriod_after) = ...
                    mean( regionArea.entireEmbryo.nbR1(image_reference_after(iPeriod_after)-param.sp2+1: image_reference_after(iPeriod_after+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
            end
        end
        
        % 3 regions
        if compute_3regions == 1
            for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                area_3regions.region1_period_before(iPeriod_before) = ...
                    mean( area_3regions.region1(image_reference_before_3regions(iPeriod_before+1)-param.sp2+1: image_reference_before_3regions(iPeriod_before)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                area_3regions.region2_period_before(iPeriod_before) = ...
                    mean( area_3regions.region2(image_reference_before_3regions(iPeriod_before+1)-param.sp2+1: image_reference_before_3regions(iPeriod_before)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                area_3regions.region3_period_before(iPeriod_before) = ...
                    mean( area_3regions.region3(image_reference_before_3regions(iPeriod_before+1)-param.sp2+1: image_reference_before_3regions(iPeriod_before)-param.sp2+1) ).* ( (param.resol/1000)^2 );
            end
            
            for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                area_3regions.region1_period_after(iPeriod_after) = ...
                    mean( area_3regions.region1(image_reference_after_3regions(iPeriod_after)-param.sp2+1: image_reference_after_3regions(iPeriod_after+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                area_3regions.region2_period_after(iPeriod_after) = ...
                    mean( area_3regions.region2(image_reference_after_3regions(iPeriod_after)-param.sp2+1: image_reference_after_3regions(iPeriod_after+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                area_3regions.region3_period_after(iPeriod_after) = ...
                    mean( area_3regions.region3(image_reference_after_3regions(iPeriod_after)-param.sp2+1: image_reference_after_3regions(iPeriod_after+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get duration histograms
        
        if compute_wholeEmbryo == 1 && ~isempty(image_reference_before)
            for iPeriod_before = 1 : length(image_reference_before)-1
                if time_reference_choice == 0
                    extension = ['before_minus' num2str(iPeriod_before)];
                elseif time_reference_choice == 1
                    extension = ['NEBD_minus' num2str(iPeriod_before)];
                end
                % entire embryo
                dataTracks_input = entireEmbryo.(extension).lengthTracks( entireEmbryo.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_before.(extension).entireEmbryo = [ [embryo_index_before.(extension).entireEmbryo] nbEmbryo ];
                    nbEmbryo_period_before_above250(1,iPeriod_before) = nbEmbryo_period_before_above250(1,iPeriod_before) + 1;
                end
                if nbEmbryo == 1 && iPeriod_before == 1
                    [ tracks_duration_histo_before ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo,[],window_size_before/param.sp6,regionArea_period_before(iPeriod_before));
                else
                    [ tracks_duration_histo_before ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo,tracks_duration_histo_before,window_size_before/param.sp6,regionArea_period_before(iPeriod_before));
                end
            end
        else
            tracks_duration_histo_before{1,nbEmbryo}.bincounts = [];
            tracks_duration_histo_before{1,nbEmbryo}.binranges =[]; % in sec
            tracks_duration_histo_before{1,nbEmbryo}.number_tracks = 0;
            tracks_duration_histo_before{1,nbEmbryo}.lengthTracks = [];
            tracks_duration_histo_before{1,nbEmbryo}.timeDuration_phase = 0;
            tracks_duration_histo_before{1,nbEmbryo}.area = 0;
        end
        
        if compute_3regions == 1 && ~isempty(image_reference_before_3regions)
            for iPeriod_before = 1 : length(image_reference_before_3regions)-1
                if time_reference_choice == 0
                    extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
                elseif time_reference_choice == 1
                    extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
                end
                % region 1
                dataTracks_input = region1.(extension).lengthTracks( region1.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_before.(extension).region1 = [ [embryo_index_before.(extension).region1] nbEmbryo_3regions ];
                    nbEmbryo_period_before_above250_3regions.region1(1,iPeriod_before) = nbEmbryo_period_before_above250_3regions.region1(1,iPeriod_before) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_before == 1
                    [ tracks_duration_histo_before_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,[],window_size_before_3regions/param.sp6,area_3regions.region1_period_before(iPeriod_before));
                else
                    [ tracks_duration_histo_before_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,tracks_duration_histo_before_region1,window_size_before_3regions/param.sp6,area_3regions.region1_period_before(iPeriod_before));
                end
                % region 2
                dataTracks_input = region2.(extension).lengthTracks( region2.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_before.(extension).region2 = [ [embryo_index_before.(extension).region2] nbEmbryo_3regions ];
                    nbEmbryo_period_before_above250_3regions.region2(1,iPeriod_before) = nbEmbryo_period_before_above250_3regions.region2(1,iPeriod_before) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_before == 1
                    [ tracks_duration_histo_before_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,[],window_size_before_3regions/param.sp6,area_3regions.region2_period_before(iPeriod_before));
                else
                    [ tracks_duration_histo_before_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,tracks_duration_histo_before_region2,window_size_before_3regions/param.sp6,area_3regions.region2_period_before(iPeriod_before));
                end
                % region 3
                dataTracks_input = region3.(extension).lengthTracks( region3.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_before.(extension).region3 = [ [embryo_index_before.(extension).region3] nbEmbryo_3regions ];
                    nbEmbryo_period_before_above250_3regions.region3(1,iPeriod_before) = nbEmbryo_period_before_above250_3regions.region3(1,iPeriod_before) + 1;
                end
                
                if nbEmbryo_3regions == 1 && iPeriod_before == 1
                    [ tracks_duration_histo_before_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,[],window_size_before_3regions/param.sp6,area_3regions.region3_period_before(iPeriod_before));
                else
                    [ tracks_duration_histo_before_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_before,nbEmbryo_3regions,tracks_duration_histo_before_region3,window_size_before_3regions/param.sp6,area_3regions.region3_period_before(iPeriod_before));
                end
                
            end
        else
            tracks_duration_histo_before_region1{1,nbEmbryo}.bincounts = [];
            tracks_duration_histo_before_region1{1,nbEmbryo}.binranges =[]; % in sec
            tracks_duration_histo_before_region1{1,nbEmbryo}.number_tracks = 0;
            tracks_duration_histo_before_region1{1,nbEmbryo}.lengthTracks = [];
            tracks_duration_histo_before_region1{1,nbEmbryo}.timeDuration_phase = 0;
            tracks_duration_histo_before_region1{1,nbEmbryo}.area = 0;
            tracks_duration_histo_before_region2{1,nbEmbryo}.bincounts = [];
            tracks_duration_histo_before_region2{1,nbEmbryo}.binranges =[]; % in sec
            tracks_duration_histo_before_region2{1,nbEmbryo}.number_tracks = 0;
            tracks_duration_histo_before_region2{1,nbEmbryo}.lengthTracks = [];
            tracks_duration_histo_before_region2{1,nbEmbryo}.timeDuration_phase = 0;
            tracks_duration_histo_before_region2{1,nbEmbryo}.area = 0;
            tracks_duration_histo_before_region3{1,nbEmbryo}.bincounts = [];
            tracks_duration_histo_before_region3{1,nbEmbryo}.binranges =[]; % in sec
            tracks_duration_histo_before_region3{1,nbEmbryo}.number_tracks = 0;
            tracks_duration_histo_before_region3{1,nbEmbryo}.lengthTracks = [];
            tracks_duration_histo_before_region3{1,nbEmbryo}.timeDuration_phase = 0;
            tracks_duration_histo_before_region3{1,nbEmbryo}.area = 0;
        end
        
        if compute_wholeEmbryo == 1 && ~isempty(image_reference_after)
            for iPeriod_after = 1 : length(image_reference_after)-1
                if time_reference_choice == 0
                    extension = ['after_plus' num2str(iPeriod_after)];
                elseif time_reference_choice == 1
                    extension = ['NEBD_plus' num2str(iPeriod_after)];
                end
                % entire embryo
                dataTracks_input = entireEmbryo.(extension).lengthTracks( entireEmbryo.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_after.(extension).entireEmbryo = [ [embryo_index_after.(extension).entireEmbryo] nbEmbryo ];
                    nbEmbryo_period_after_above250(1,iPeriod_after) = nbEmbryo_period_after_above250(1,iPeriod_after) + 1;
                end
                if nbEmbryo == 1 && iPeriod_after == 1
                    [ tracks_duration_histo_after ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo,[],window_size_after/param.sp6,regionArea_period_after(iPeriod_after));
                else
                    [ tracks_duration_histo_after ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo,tracks_duration_histo_after,window_size_after/param.sp6,regionArea_period_after(iPeriod_after));
                end
            end
        end
        
        if compute_3regions == 1 && ~isempty(image_reference_after_3regions)
            for iPeriod_after = 1 : length(image_reference_after_3regions)-1
                if time_reference_choice == 0
                    extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
                elseif time_reference_choice == 1
                    extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
                end
                % region 1
                dataTracks_input = region1.(extension).lengthTracks( region1.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_after.(extension).region1 = [ [embryo_index_after.(extension).region1] nbEmbryo_3regions ];
                    nbEmbryo_period_after_above250_3regions.region1(1,iPeriod_after) = nbEmbryo_period_after_above250_3regions.region1(1,iPeriod_after) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_after == 1
                    [ tracks_duration_histo_after_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,[],window_size_after_3regions/param.sp6,area_3regions.region1_period_after(iPeriod_after));
                else
                    [ tracks_duration_histo_after_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,tracks_duration_histo_after_region1,window_size_after_3regions/param.sp6,area_3regions.region1_period_after(iPeriod_after));
                end
                % region 2
                dataTracks_input = region2.(extension).lengthTracks( region2.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_after.(extension).region2 = [ [embryo_index_after.(extension).region2] nbEmbryo_3regions ];
                    nbEmbryo_period_after_above250_3regions.region2(1,iPeriod_after) = nbEmbryo_period_after_above250_3regions.region2(1,iPeriod_after) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_after == 1
                    [ tracks_duration_histo_after_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,[],window_size_after_3regions/param.sp6,area_3regions.region2_period_after(iPeriod_after));
                else
                    [ tracks_duration_histo_after_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,tracks_duration_histo_after_region2,window_size_after_3regions/param.sp6,area_3regions.region2_period_after(iPeriod_after));
                end
                % region 3
                dataTracks_input = region3.(extension).lengthTracks( region3.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_after.(extension).region3 = [ [embryo_index_after.(extension).region3] nbEmbryo_3regions ];
                    nbEmbryo_period_after_above250_3regions.region3(1,iPeriod_after) = nbEmbryo_period_after_above250_3regions.region3(1,iPeriod_after) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_after == 1
                    [ tracks_duration_histo_after_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,[],window_size_after_3regions/param.sp6,area_3regions.region3_period_after(iPeriod_after));
                else
                    [ tracks_duration_histo_after_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_after,nbEmbryo_3regions,tracks_duration_histo_after_region3,window_size_after_3regions/param.sp6,area_3regions.region3_period_after(iPeriod_after));
                end
            end
        end
        if compute_3regions == 1
            clear region1 region2 region3
            clear xStart_3regions xEnd_3regions area_3regions
        end
        if compute_wholeEmbryo == 1
            clear entireEmbryo
            clear regionArea_period_after regionArea_period_before area_blastomere
        end
        clear regionArea regionXlength regionXlimit dataTracks_input dataTracks_rotated furrow_position
    end
end

%% save structure with histo
if compute_wholeEmbryo == 1
    name = strcat('tracks_duration_histo_before_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_before');
    
    name = strcat('tracks_duration_histo_after_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_after');
end

if compute_3regions == 1
    name = strcat('tracks_duration_histo_before_region1_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_before_region1');
    
    name = strcat('tracks_duration_histo_after_region1_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_after_region1');
    
    name = strcat('tracks_duration_histo_before_region2_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_before_region2');
    
    name = strcat('tracks_duration_histo_after_region2_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_after_region2');
    
    name = strcat('tracks_duration_histo_before_region3_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_before_region3');
    
    name = strcat('tracks_duration_histo_after_region3_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_after_region3');
end


%% fitting for each condition

if compute_wholeEmbryo == 1
    for iPeriod_before = 1 : max_period_before
        
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before)];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before)];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250(iPeriod_before);
        
        if nbEmbryo_givenCondition > 1
            
            % entire embryo
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_before.(extension).entireEmbryo(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_before{iPeriod_before,iEmbryo}.binranges, ...
                    tracks_duration_histo_before{iPeriod_before,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_before{iPeriod_before,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_before{iPeriod_before,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_before{iPeriod_before,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'entireEmbryo',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'entireEmbryo',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'entireEmbryo',save_stem,2,0 );
            end
            final_fitting_all.(extension).entireEmbryo.fitting = fitting_results;
            final_fitting_all.(extension).entireEmbryo.model = model_choice;
            clear fitting_results model_choice input
            close all
        end
    end
end

if compute_3regions == 1
    for iPeriod_before = 1 : max_period_before_3regions
        
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region1(iPeriod_before);
        
        if nbEmbryo_givenCondition > 1
            
            % region1
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_before.(extension).region1(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_before_region1{iPeriod_before,iEmbryo}.binranges, ...
                    tracks_duration_histo_before_region1{iPeriod_before,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_before_region1{iPeriod_before,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_before_region1{iPeriod_before,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_before_region1{iPeriod_before,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region1',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region1',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region1',save_stem,2,0 );
            end
            final_fitting_all.(extension).region1.fitting = fitting_results;
            final_fitting_all.(extension).region1.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
            
              % extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region2(iPeriod_before);
        
        if nbEmbryo_givenCondition > 1
            % region2
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_before.(extension).region2(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_before_region2{iPeriod_before,iEmbryo}.binranges, ...
                    tracks_duration_histo_before_region2{iPeriod_before,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_before_region2{iPeriod_before,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_before_region2{iPeriod_before,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_before_region2{iPeriod_before,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region2',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region2',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region2',save_stem,2,0 );
            end
            final_fitting_all.(extension).region2.fitting = fitting_results;
            final_fitting_all.(extension).region2.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
        
              % extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region3(iPeriod_before);
        
        if nbEmbryo_givenCondition > 1
            
            % region3
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_before.(extension).region3(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_before_region3{iPeriod_before,iEmbryo}.binranges, ...
                    tracks_duration_histo_before_region3{iPeriod_before,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_before_region3{iPeriod_before,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_before_region3{iPeriod_before,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_before_region3{iPeriod_before,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region3',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region3',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region3',save_stem,2,0 );
            end
            final_fitting_all.(extension).region3.fitting = fitting_results;
            final_fitting_all.(extension).region3.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end

if compute_wholeEmbryo == 1
    for iPeriod_after = 1 : max_period_after
        
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after)];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after)];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250(iPeriod_after);
        
        if nbEmbryo_givenCondition > 1
            
            % entire embryo
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_after.(extension).entireEmbryo(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_after{iPeriod_after,iEmbryo}.binranges, ...
                    tracks_duration_histo_after{iPeriod_after,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_after{iPeriod_after,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_after{iPeriod_after,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_after{iPeriod_after,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'entireEmbryo',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'entireEmbryo',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'entireEmbryo',save_stem,2,0 );
            end
            final_fitting_all.(extension).entireEmbryo.fitting = fitting_results;
            final_fitting_all.(extension).entireEmbryo.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end

if compute_3regions == 1
    for iPeriod_after = 1 : max_period_after_3regions
        
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region1(iPeriod_after);
        
        if nbEmbryo_givenCondition > 1
            % region1
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_after.(extension).region1(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_after_region1{iPeriod_after,iEmbryo}.binranges, ...
                    tracks_duration_histo_after_region1{iPeriod_after,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_after_region1{iPeriod_after,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_after_region1{iPeriod_after,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_after_region1{iPeriod_after,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region1',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region1',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region1',save_stem,2,0 );
            end
            final_fitting_all.(extension).region1.fitting = fitting_results;
            final_fitting_all.(extension).region1.model = model_choice;
            clear fitting_results model_choice input
            close all
        end
        
              % extension = ['after_minus' num2str(iPeriod_after) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region2(iPeriod_after);
        
        if nbEmbryo_givenCondition > 1
            
            % region2
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_after.(extension).region2(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_after_region2{iPeriod_after,iEmbryo}.binranges, ...
                    tracks_duration_histo_after_region2{iPeriod_after,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_after_region2{iPeriod_after,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_after_region2{iPeriod_after,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_after_region2{iPeriod_after,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region2',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region2',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region2',save_stem,2,0 );
            end
            final_fitting_all.(extension).region2.fitting = fitting_results;
            final_fitting_all.(extension).region2.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
        
               %extension = ['after_minus' num2str(iPeriod_after) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region3(iPeriod_after);
        
        if nbEmbryo_givenCondition > 1
            
            % region3
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_after.(extension).region3(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_after_region3{iPeriod_after,iEmbryo}.binranges, ...
                    tracks_duration_histo_after_region3{iPeriod_after,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_after_region3{iPeriod_after,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_after_region3{iPeriod_after,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_after_region3{iPeriod_after,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            if error_compute == 0
                to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region3',save_stem,2 );
            elseif error_compute == 1
                [ fitting_results ] = to_calculate_errors_parameters_using_Boostrap_Mono_Dble_Expo...
                    ( input,models,fitting_results,best_model,extension,'region3',nbEmbryo_givenCondition,...
                    minLength_tracks,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,0 );
                to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,extension,'region3',save_stem,2,0 );
            end
            final_fitting_all.(extension).region3.fitting = fitting_results;
            final_fitting_all.(extension).region3.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end


%% save structure with results incorporated

name = strcat('final_results-BayesianInference_sum_mle' , xmlfile_bkp, '.mat');
save(fullfile(main_path,name), '-struct','final_fitting_all');


%% display results as temporal modulation plot

if compute_wholeEmbryo == 1
    time_range_before = [];
    count_period_before = 0; % to count only period for which nb embryo above 1
    for iPeriod_before = max_period_before : -1 : 1
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before)];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before)];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250(iPeriod_before);
        if nbEmbryo_givenCondition > 1
            count_period_before = count_period_before +1;
            time_range_before =  [ [time_range_before] -window_size_before * iPeriod_before ];
        end
    end
    time_range_after = [];
    count_period_after = 0;
    for iPeriod_after = 1: 1: max_period_after
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after)];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after)];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250(iPeriod_after);
        if nbEmbryo_givenCondition > 1
            count_period_after = count_period_after +1;
            time_range_after =  [ [time_range_after] window_size_after * iPeriod_after ];
        end
    end
    all_time_range = [ time_range_before ,0, time_range_after ]./param.sp6;
    if additionnal_ref_nebd == 1
        all_time_range_nebd = all_time_range +( param.delta_furrowDetection_nebd - param.delta_furrowDetection_anaphaseOnset );
    end
end

if compute_3regions == 1
    time_range_before_3regions.region1 = [];
    count_period_before_3regions.region1 = 0; % to count only period for which nb embryo above 1
    time_range_before_3regions.region2 = [];
    count_period_before_3regions.region2 = 0; % to count only period for which nb embryo above 1
    time_range_before_3regions.region3 = [];
    count_period_before_3regions.region3 = 0; % to count only period for which nb embryo above 1
    for iPeriod_before = max_period_before_3regions : -1 : 1
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region1(iPeriod_before);
        if nbEmbryo_givenCondition > 1
            count_period_before_3regions.region1 = count_period_before_3regions.region1 +1;
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region2(iPeriod_before);
        if nbEmbryo_givenCondition > 1
            count_period_before_3regions.region2 = count_period_before_3regions.region2 +1;
        end
        nbEmbryo_givenCondition = nbEmbryo_period_before_above250_3regions.region3(iPeriod_before);
        if nbEmbryo_givenCondition > 1
            count_period_before_3regions.region3 = count_period_before_3regions.region3 +1;
        end        
    end
    % get max nb of metaphase period wuthin the three regions
    max_meta_period = max([ count_period_before_3regions.region1 count_period_before_3regions.region2 count_period_before_3regions.region3 ]);
     for iPeriod_before = max_meta_period : -1 : 1
        time_range_before_3regions.region1 =  [ [time_range_before_3regions.region1] -window_size_before_3regions * iPeriod_before ];
        time_range_before_3regions.region2 =  [ [time_range_before_3regions.region2] -window_size_before_3regions * iPeriod_before ];
        time_range_before_3regions.region3 =  [ [time_range_before_3regions.region3] -window_size_before_3regions * iPeriod_before ];      
     end
     
     count_period_before_3regions.region1 = max_meta_period;
     count_period_before_3regions.region2 = max_meta_period;
     count_period_before_3regions.region3 = max_meta_period;
     
    time_range_after_3regions.region1 = [];
    count_period_after_3regions.region1 = 0;
    time_range_after_3regions.region2 = [];
    count_period_after_3regions.region2 = 0;
    time_range_after_3regions.region3 = [];
    count_period_after_3regions.region3 = 0;    
    for iPeriod_after = 1: 1: max_period_after_3regions
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
        end
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region1(iPeriod_after);
        if nbEmbryo_givenCondition > 1
            count_period_after_3regions.region1 = count_period_after_3regions.region1 +1;
        end
         nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region2(iPeriod_after);
        if nbEmbryo_givenCondition > 1
            count_period_after_3regions.region2 = count_period_after_3regions.region2 +1;
        end
        nbEmbryo_givenCondition = nbEmbryo_period_after_above250_3regions.region3(iPeriod_after);
        if nbEmbryo_givenCondition > 1
            count_period_after_3regions.region3 = count_period_after_3regions.region3 +1;
        end        
    end
    % get max nb of metaphase period wuthin the three regions
    max_ana_period = max([ count_period_after_3regions.region1 count_period_after_3regions.region2 count_period_after_3regions.region3 ]);
     for iPeriod_after = 1 : 1 : max_ana_period
        time_range_after_3regions.region1 =  [ [time_range_after_3regions.region1] window_size_after_3regions * iPeriod_after ];
        time_range_after_3regions.region2 =  [ [time_range_after_3regions.region2] window_size_after_3regions * iPeriod_after ];
        time_range_after_3regions.region3 =  [ [time_range_after_3regions.region3] window_size_after_3regions * iPeriod_after ];    
     end
    count_period_after_3regions.region1 = max_ana_period;
    count_period_after_3regions.region2 = max_ana_period;
    count_period_after_3regions.region3 = max_ana_period;
     
    all_time_range_3regions.region1 = [ time_range_before_3regions.region1 ,0, time_range_after_3regions.region1 ]./param.sp6;
    all_time_range_3regions.region2 = [ time_range_before_3regions.region2 ,0, time_range_after_3regions.region2 ]./param.sp6;
    all_time_range_3regions.region3 = [ time_range_before_3regions.region3 ,0, time_range_after_3regions.region3 ]./param.sp6;
    
    if additionnal_ref_nebd == 1
         all_time_range_3regions_nebd.region1 =  all_time_range_3regions.region1 +( param.delta_furrowDetection_nebd - param.delta_furrowDetection_anaphaseOnset );
         all_time_range_3regions_nebd.region2 =  all_time_range_3regions.region2 +( param.delta_furrowDetection_nebd - param.delta_furrowDetection_anaphaseOnset );
         all_time_range_3regions_nebd.region3 =  all_time_range_3regions.region3 +( param.delta_furrowDetection_nebd - param.delta_furrowDetection_anaphaseOnset );
    end
    
end



% plot to display temporal regulation
if compute_wholeEmbryo == 1
    vector_T1 = nan( length(all_time_range)-1,1 );
    vector_T2 = nan( length(all_time_range)-1,1 );
    vector_P1 = nan( length(all_time_range)-1,1 );
    vector_P2 = nan( length(all_time_range)-1,1 );
    vector_N1 = nan( length(all_time_range)-1,1);
    vector_N2 = nan( length(all_time_range)-1,1 );
    vector_d1 = nan( length(all_time_range)-1,1 );
    vector_d2 = nan( length(all_time_range)-1,1 );
    vector_f1 = nan( length(all_time_range)-1,1 );
    vector_f2 = nan( length(all_time_range)-1,1 );
    if error_compute == 1
        vector_T1_std = nan( length(all_time_range)-1,1 );
        vector_T2_std = nan( length(all_time_range)-1,1 );
        vector_P1_std = nan( length(all_time_range)-1,1 );
        vector_P2_std = nan( length(all_time_range)-1,1 );
        vector_N1_std = nan( length(all_time_range)-1,1);
        vector_N2_std = nan( length(all_time_range)-1,1 );
        vector_d1_std = nan( length(all_time_range)-1,1);
        vector_d2_std = nan( length(all_time_range)-1,1 );
        vector_f1_std = nan( length(all_time_range)-1,1);
        vector_f2_std = nan( length(all_time_range)-1,1 );     
    end
end

if compute_3regions == 1
    vector_T1_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_T1_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_T1_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_T2_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_T2_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_T2_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_P1_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_P1_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_P1_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_P2_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_P2_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_P2_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_N1_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_N1_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_N1_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_N2_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_N2_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_N2_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_d1_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_d1_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_d1_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_d2_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_d2_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_d2_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_f1_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_f1_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_f1_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    vector_f2_region1 = nan( length(all_time_range_3regions.region1)-1,1);
    vector_f2_region2 = nan( length(all_time_range_3regions.region2)-1,1);
    vector_f2_region3 = nan( length(all_time_range_3regions.region3)-1,1);
    
    if error_compute == 1
        vector_T1_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_T1_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_T1_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_T2_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_T2_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_T2_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_P1_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_P1_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_P1_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_P2_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_P2_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_P2_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_N1_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_N1_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_N1_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_N2_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_N2_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_N2_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_d1_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_d1_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_d1_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_d2_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_d2_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_d2_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_f1_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_f1_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_f1_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
        vector_f2_region1_std = nan( length(all_time_range_3regions.region1)-1,1);
        vector_f2_region2_std = nan( length(all_time_range_3regions.region2)-1,1);
        vector_f2_region3_std = nan( length(all_time_range_3regions.region3)-1,1);
    end
end

if compute_wholeEmbryo == 1
    for iPeriod_before = count_period_before : -1 : 1
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before)];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before)];
        end
        
        if isfield(final_fitting_all,extension) == 1
            
        vector_T1(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T1;
        vector_T2(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T2;
        if error_compute == 1
            vector_T1_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T1_se_bootstrap;
            vector_T2_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T2_se_bootstrap;          
        end
        
        vector_P1(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1;
        vector_P2(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2;
        if error_compute == 1
            vector_P1_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1_se_bootstrap;
            vector_P2_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2_se_bootstrap;             
        end
        
        if window_size_after ~= window_size_before
            multipicative_factor = window_size_after / window_size_before;
            vector_N1(count_period_before-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total1 ...
                ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
            vector_N2(count_period_before-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total2 ...
                ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
        elseif window_size_after == window_size_before
            vector_N1(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total1 ...
                ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
            vector_N2(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total2 ...
                ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
        end
        if error_compute == 1
            vector_N1_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1_se_bootstrap .* ...
                final_fitting_all.(extension).entireEmbryo.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
            vector_N2_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2_se_bootstrap .* ...
                final_fitting_all.(extension).entireEmbryo.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
        end
        
        vector_d1(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
        vector_d2(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
        
        if error_compute == 1
            vector_d1_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.d1_se_bootstrap; 
            vector_d2_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.d2_se_bootstrap;
        end
        
        vector_f1(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_time.mean1 ;
        vector_f2(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_time.mean2;
        if error_compute == 1
            vector_f1_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.f1_se_bootstrap; 
            vector_f2_std(count_period_before-iPeriod_before+1,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.f2_se_bootstrap; 
        end
        
        else
            
            vector_T1(count_period_before-iPeriod_before+1,1) = NaN;
            vector_T2(count_period_before-iPeriod_before+1,1) = NaN;
            vector_P1(count_period_before-iPeriod_before+1,1) = NaN;
            vector_P2(count_period_before-iPeriod_before+1,1) = NaN;
            vector_N1(count_period_before-iPeriod_before+1,1) = NaN;
            vector_N2(count_period_before-iPeriod_before+1,1) = NaN;
            vector_d1(count_period_before-iPeriod_before+1,1) = NaN;
            vector_d2(count_period_before-iPeriod_before+1,1) = NaN;
            vector_f1(count_period_before-iPeriod_before+1,1) = NaN;
            vector_f2(count_period_before-iPeriod_before+1,1) = NaN;
            if error_compute == 1
                vector_T1_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_T2_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_P1_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_P2_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_N1_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_N2_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_d1_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_d2_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_f1_std(count_period_before-iPeriod_before+1,1) = NaN;
                vector_f2_std(count_period_before-iPeriod_before+1,1) = NaN;
            end
        end
    end
end

if compute_3regions == 1
    for iPeriod_before = count_period_before_3regions.region1 : -1 : 1
        if time_reference_choice == 0
            extension = ['before_minus' num2str(iPeriod_before) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_minus' num2str(iPeriod_before) '_3regions'];
        end
        if isfield(final_fitting_all.(extension),'region1')
            vector_T1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T1;
            vector_T2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T2;
            vector_P1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1;
            vector_P2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
            end
            vector_d1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_T2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_P1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_P2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_N1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_N2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_d1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_d2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_f1_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            vector_f2_region1(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
        end

        if isfield(final_fitting_all.(extension),'region2')
            vector_T1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T1;
            vector_T2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T2;
            vector_P1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P1;
            vector_P2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
            end
            vector_d1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_T2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_P1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_P2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_N1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_N2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_d1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_d2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_f1_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
            vector_f2_region2(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
        end

        if isfield(final_fitting_all.(extension),'region3')
            vector_T1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T1;
            vector_T2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T2;
            vector_P1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1;
            vector_P2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = multipicative_factor * final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
            end
            vector_d1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_T2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_P1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_P2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_N1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_N2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_d1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_d2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_f1_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            vector_f2_region3(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
        end
        
        
        if error_compute == 1

            if isfield(final_fitting_all.(extension),'region1')
                vector_T1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2_se_bootstrap;             
                vector_N1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region1.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region1.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_d1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.d1_se_bootstrap; 
                vector_d2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.d2_se_bootstrap;
                vector_f1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.f2_se_bootstrap;                 
           else
                vector_T1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_T2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_P1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_P2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_N1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_N2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_d1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_d2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_f1_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
                vector_f2_region1_std(count_period_before_3regions.region1-iPeriod_before+1,1) = NaN;
            end

            if isfield(final_fitting_all.(extension),'region2')
                vector_T1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2_se_bootstrap;
                vector_N1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region2.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region2.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_d1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.d1_se_bootstrap; 
                vector_d2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.d2_se_bootstrap; 
                vector_f1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.f2_se_bootstrap; 
          else
                vector_T1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_T2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_P1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_P2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_N1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_N2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_d1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_d2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_f1_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;
                vector_f2_region2_std(count_period_before_3regions.region2-iPeriod_before+1,1) = NaN;      
            
            end

            if isfield(final_fitting_all.(extension),'region3')
                vector_T1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2_se_bootstrap;
                vector_N1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region3.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region3.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_d1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.d1_se_bootstrap;
                vector_d2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.d2_se_bootstrap; 
                vector_f1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.f2_se_bootstrap; 
            else
                vector_T1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_T2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_P1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_P2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_N1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_N2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_d1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_d2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_f1_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
                vector_f2_region3_std(count_period_before_3regions.region3-iPeriod_before+1,1) = NaN;
            end
        end
        
    end
end

if compute_wholeEmbryo == 1
    for iPeriod_after = 1 : count_period_after
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after)];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after)];
        end
        
        if isfield(final_fitting_all,extension) == 1
            
        vector_T1(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T1;
        vector_T2(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T2;
        if error_compute == 1
            vector_T1_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T1_se_bootstrap;
            vector_T2_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.T2_se_bootstrap;        
        end
        
        vector_P1(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1;
        vector_P2(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2;
        if error_compute == 1
            vector_P1_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1_se_bootstrap;
            vector_P2_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2_se_bootstrap;
        end
        
        vector_N1(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total1 ...
            ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
        vector_N2(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population.total2 ...
            ./ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
        
        if error_compute == 1
            vector_N1_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P1_se_bootstrap .* ...
                final_fitting_all.(extension).entireEmbryo.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);
            vector_N2_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.P2_se_bootstrap .* ...
                final_fitting_all.(extension).entireEmbryo.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).entireEmbryo.fitting.size_population)) -1);           
        end
        
        vector_d1(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
        vector_d2(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
        
        if error_compute == 1
            vector_d1_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.d1_se_bootstrap;
            vector_d2_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.d2_se_bootstrap;            
        end
        
        vector_f1(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_time.mean1 ;
        vector_f2(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.size_population_normalized_time.mean2;
        
        if error_compute == 1
            vector_f1_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.f1_se_bootstrap ;
            vector_f2_std(count_period_before+iPeriod_after,1) = final_fitting_all.(extension).entireEmbryo.fitting.DoubleExpo.f2_se_bootstrap;            
        end
        
        else
            vector_T1(count_period_before+iPeriod_after,1) = NaN;
            vector_T2(count_period_before+iPeriod_after,1) = NaN;
            vector_P1(count_period_before+iPeriod_after,1) = NaN;
            vector_P2(count_period_before+iPeriod_after,1) = NaN;
            vector_N1(count_period_before+iPeriod_after,1) = NaN;
            vector_N2(count_period_before+iPeriod_after,1) = NaN;
            vector_d1(count_period_before+iPeriod_after,1) = NaN;
            vector_d2(count_period_before+iPeriod_after,1) = NaN;
            vector_f1(count_period_before+iPeriod_after,1) = NaN;
            vector_f2(count_period_before+iPeriod_after,1) = NaN;
            if error_compute == 1
                vector_T1_std(count_period_before+iPeriod_after,1) = NaN;
                vector_T2_std(count_period_before+iPeriod_after,1) = NaN;
                vector_P1_std(count_period_before+iPeriod_after,1) = NaN;
                vector_P2_std(count_period_before+iPeriod_after,1) = NaN;
                vector_N1_std(count_period_before+iPeriod_after,1) = NaN;
                vector_N2_std(count_period_before+iPeriod_after,1) = NaN;
                vector_d1_std(count_period_before+iPeriod_after,1) = NaN;
                vector_d2_std(count_period_before+iPeriod_after,1) = NaN;
                vector_f1_std(count_period_before+iPeriod_after,1) = NaN;
                vector_f2_std(count_period_before+iPeriod_after,1) = NaN;
            end
        end
    end
end


if compute_3regions == 1
    for iPeriod_after = 1 : count_period_after_3regions.region1
        if time_reference_choice == 0
            extension = ['after_plus' num2str(iPeriod_after) '_3regions'];
        elseif time_reference_choice == 1
            extension = ['NEBD_plus' num2str(iPeriod_after) '_3regions'];
        end

        if isfield(final_fitting_all.(extension),'region1')
            vector_T1_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T1;
            vector_T2_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T2;
            vector_P1_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1;
            vector_P2_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region1(count_period_before_3regions.region1+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1(count_period_before_3regions.region1+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
            end
            vector_d1_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region1(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_T2_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_P1_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_P2_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_N1_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_N2_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_d1_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_d2_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_f1_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
            vector_f2_region1(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
        end

        if isfield(final_fitting_all.(extension),'region2')
            vector_T1_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T1;
            vector_T2_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T2;
            vector_P1_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P1;
            vector_P2_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region2(count_period_before_3regions.region2+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2(count_period_before_3regions.region2+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
            end
            vector_d1_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region2(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_T2_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_P1_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_P2_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_N1_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_N2_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_d1_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_d2_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_f1_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
            vector_f2_region2(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
        end

        if isfield(final_fitting_all.(extension),'region3')
            vector_T1_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T1;
            vector_T2_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T2;
            vector_P1_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1;
            vector_P2_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2;
            if window_size_after_3regions ~= window_size_before_3regions
                multipicative_factor = window_size_after_3regions / window_size_before_3regions;
                vector_N1_region3(count_period_before_3regions.region3+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3(count_period_before_3regions.region3+iPeriod_after,1) = multipicative_factor * final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
            elseif window_size_after_3regions == window_size_before_3regions
                vector_N1_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total1 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population.total2 ...
                    ./ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
            end
            vector_d1_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean1;
            vector_d2_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_timeAndArea.mean2;
            vector_f1_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_time.mean1;
            vector_f2_region3(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.size_population_normalized_time.mean2;
        else
            vector_T1_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_T2_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_P1_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_P2_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_N1_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_N2_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_d1_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_d2_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_f1_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
            vector_f2_region3(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
        end
        
        
        if error_compute == 1

            if isfield(final_fitting_all.(extension),'region1')
                vector_T1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2_se_bootstrap;             
                vector_N1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region1.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_N2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region1.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region1.fitting.size_population)) -1);
                vector_d1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.d1_se_bootstrap;
                vector_d2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.d2_se_bootstrap; 
                vector_f1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = final_fitting_all.(extension).region1.fitting.DoubleExpo.f2_se_bootstrap;               
            else
                vector_T1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_T2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_P1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_P2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_N1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_N2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_d1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_d2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_f1_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;
                vector_f2_region1_std(count_period_before_3regions.region1+iPeriod_after,1) = NaN;             
            end

            if isfield(final_fitting_all.(extension),'region2')
                vector_T1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2_se_bootstrap;
                vector_N1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region2.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_N2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region2.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region2.fitting.size_population)) -1);
                vector_d1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.d1_se_bootstrap; 
                vector_d2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.d2_se_bootstrap; 
                vector_f1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = final_fitting_all.(extension).region2.fitting.DoubleExpo.f2_se_bootstrap;                
            else
                vector_T1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_T2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_P1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_P2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_N1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_N2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_d1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_d2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_f1_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;
                vector_f2_region2_std(count_period_before_3regions.region2+iPeriod_after,1) = NaN;             
            end

            if isfield(final_fitting_all.(extension),'region3')
                vector_T1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T1_se_bootstrap;
                vector_T2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.T2_se_bootstrap;
                vector_P1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1_se_bootstrap;
                vector_P2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2_se_bootstrap;
                vector_N1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P1_se_bootstrap .* ...
                    final_fitting_all.(extension).region3.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_N2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.P2_se_bootstrap .* ...
                    final_fitting_all.(extension).region3.fitting.size_population.total.raw/ ( numel(fieldnames(final_fitting_all.(extension).region3.fitting.size_population)) -1);
                vector_d1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.d1_se_bootstrap; 
                vector_d2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.d2_se_bootstrap; 
                vector_f1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.f1_se_bootstrap; 
                vector_f2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = final_fitting_all.(extension).region3.fitting.DoubleExpo.f2_se_bootstrap;                
            else
                vector_T1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_T2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_P1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_P2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_N1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_N2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_d1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_d2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_f1_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;
                vector_f2_region3_std(count_period_before_3regions.region3+iPeriod_after,1) = NaN;                
            end
        end
        
    end
end


if compute_wholeEmbryo == 1
    for i = 1 : length(all_time_range)-1
        all_time_range_bis(i) = (all_time_range(i)+all_time_range(i+1))/2;
    end
    if additionnal_ref_nebd == 1
        for i = 1 : length(all_time_range_nebd)-1
            all_time_range_nebd_bis(i) = (all_time_range_nebd(i)+all_time_range_nebd(i+1))/2;
        end
    end
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis,vector_T1,'-or');
        hold all
        plot(all_time_range_bis,vector_T2,'-ob');
    elseif error_compute == 1
        errorbar(all_time_range_bis,vector_T1,vector_T1_std,'-or');
        hold all
        errorbar(all_time_range_bis,vector_T2,vector_T2_std,'-ob');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('Lifetime (s)');
    legend('T1','T2')
    saveas(gcf,[save_stem 'Temporal_modulation_lifetime.fig']);
    
  
    %------------------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis,vector_P1.*100,'-xr');
        hold all
        plot(all_time_range_bis,vector_P2.*100,'-xb');
    elseif error_compute == 1
        errorbar(all_time_range_bis,vector_P1.*100,vector_P1_std.*100,'-xr');
        hold all
        errorbar(all_time_range_bis,vector_P2.*100,vector_P2_std.*100,'-xb');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end    
    ylim([0 100])
    ylabel('Proportion (%)');
    legend('P1','P2')
    saveas(gcf,[save_stem 'Temporal_modulation_proportion.fig']);
    
    
    %-----------------------------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis,vector_N1,'-or');
        hold all
        plot(all_time_range_bis,vector_N2,'-ob');
    elseif error_compute == 1
        errorbar(all_time_range_bis,vector_N1,vector_N1_std,'-or');
        hold all
        errorbar(all_time_range_bis,vector_N2,vector_N2_std,'-ob');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT count (a.u.)');
    legend('N1','N2')
    saveas(gcf,[save_stem 'Temporal_modulation_count_entireEmbryo.fig']);
       
        %-------------------------------
        
    figure,
    if error_compute == 0
        plot(all_time_range_bis,vector_d1,'-or');
        hold all
        plot(all_time_range_bis,vector_d2,'-ob');
    elseif error_compute == 1
        errorbar(all_time_range_bis,vector_d1,vector_d1_std,'-or');
        hold all
        errorbar(all_time_range_bis,vector_d2,vector_d2_std,'-ob');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT density (/min/um2)');
    legend('d1','d2')
    saveas(gcf,[save_stem 'Temporal_modulation_density_entireEmbryo.fig']);
    
    %----------------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis,vector_f1,'-or');
        hold all
        plot(all_time_range_bis,vector_f2,'-ob');
    elseif error_compute == 1
        errorbar(all_time_range_bis,vector_f1,vector_f1_std,'-or');
        hold all
        errorbar(all_time_range_bis,vector_f2,vector_f2_std,'-ob');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT frequency (/s)');
    legend('f1','f2')
    saveas(gcf,[save_stem 'Temporal_modulation_frequency_entireEmbryo.fig']);
        
    %-------------------------
    % display also data with ref time = nebd (still using onset of furrow
    % ingression as reference at the cortex)
    
    if additionnal_ref_nebd == 1
        
        figure,
        if error_compute == 0
            plot(all_time_range_nebd_bis,vector_T1,'-or');
            hold all
            plot(all_time_range_nebd_bis,vector_T2,'-ob');
        elseif error_compute == 1
            errorbar(all_time_range_nebd_bis,vector_T1,vector_T1_std,'-or');
            hold all
            errorbar(all_time_range_nebd_bis,vector_T2,vector_T2_std,'-ob');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('Lifetime (s)');
        legend('T1','T2')
        saveas(gcf,[save_stem 'Temporal_modulation_lifetime_ref-nebd.fig']);
                
        figure,
        if error_compute == 0
            plot(all_time_range_nebd_bis,vector_P1.*100,'-xr');
            hold all
            plot(all_time_range_nebd_bis,vector_P2.*100,'-xb');
        elseif error_compute == 1
            errorbar(all_time_range_nebd_bis,vector_P1.*100,vector_P1_std.*100,'-xr');
            hold all
            errorbar(all_time_range_nebd_bis,vector_P2.*100,vector_P2_std.*100,'-xb');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylim([0 100])
        ylabel('Proportion (%)');
        legend('P1','P2')
        saveas(gcf,[save_stem 'Temporal_modulation_proportion_nebd-ref.fig']);
             
                  
        figure,
        if error_compute == 0
            plot(all_time_range_nebd_bis,vector_N1,'-or');
            hold all
            plot(all_time_range_nebd_bis,vector_N2,'-ob');
        elseif error_compute == 1
            errorbar(all_time_range_nebd_bis,vector_N1,vector_N1_std,'-or');
            hold all
            errorbar(all_time_range_nebd_bis,vector_N2,vector_N2_std,'-ob');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT count (a.u.)');
        legend('N1','N2')
        saveas(gcf,[save_stem 'Temporal_modulation_count_entireEmbryo_nebd-ref.fig']);
        
       
        figure,
        if error_compute == 0
            plot(all_time_range_nebd_bis,vector_d1,'-or');
            hold all
            plot(all_time_range_nebd_bis,vector_d2,'-ob');
        elseif error_compute == 1
            errorbar(all_time_range_nebd_bis,vector_d1,vector_d1_std,'-or');
            hold all
            errorbar(all_time_range_nebd_bis,vector_d2,vector_d2_std,'-ob');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT density (/min/um2)');
        legend('d1','d2')
        saveas(gcf,[save_stem 'Temporal_modulation_density_entireEmbryo_nebd-ref.fig']);
        
        %---------------------
        
        figure,
        if error_compute == 0
            plot(all_time_range_nebd_bis,vector_f1,'-or');
            hold all
            plot(all_time_range_nebd_bis,vector_f2,'-ob');
        elseif error_compute == 1
            errorbar(all_time_range_nebd_bis,vector_f1,vector_f1_std,'-or');
            hold all
            errorbar(all_time_range_nebd_bis,vector_f2,vector_f2_std,'-ob');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT frequency (/s)');
        legend('f1','f2')
        saveas(gcf,[save_stem 'Temporal_modulation_frequency_entireEmbryo_nebd-ref.fig']);      
    end
    
end


if compute_3regions == 1
    
    for i = 1 : length(all_time_range_3regions.region1)-1
        all_time_range_bis_3regions.region1(i) = (all_time_range_3regions.region1(i)+all_time_range_3regions.region1(i+1))/2;
    end
    for i = 1 : length(all_time_range_3regions.region2)-1
        all_time_range_bis_3regions.region2(i) = (all_time_range_3regions.region2(i)+all_time_range_3regions.region2(i+1))/2;
    end
    for i = 1 : length(all_time_range_3regions.region3)-1
        all_time_range_bis_3regions.region3(i) = (all_time_range_3regions.region3(i)+all_time_range_3regions.region3(i+1))/2;
    end
    
    if additionnal_ref_nebd == 1
        for i = 1 : length(all_time_range_3regions_nebd.region1)-1
            all_time_range_bis_3regions_nebd.region1(i) = (all_time_range_3regions_nebd.region1(i)+all_time_range_3regions_nebd.region1(i+1))/2;
        end
        for i = 1 : length(all_time_range_3regions_nebd.region2)-1
            all_time_range_bis_3regions_nebd.region2(i) = (all_time_range_3regions_nebd.region2(i)+all_time_range_3regions_nebd.region2(i+1))/2;
        end
        for i = 1 : length(all_time_range_3regions_nebd.region3)-1
            all_time_range_bis_3regions_nebd.region3(i) = (all_time_range_3regions_nebd.region3(i)+all_time_range_3regions_nebd.region3(i+1))/2;
        end        
    end
    
    %--------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis_3regions.region1,vector_N1_region1,'-or');
        hold all
        plot(all_time_range_bis_3regions.region1,vector_N2_region1,'-ob');
        plot(all_time_range_bis_3regions.region2,vector_N1_region2,'-xr');
        plot(all_time_range_bis_3regions.region2,vector_N2_region2,'-xb');
        plot(all_time_range_bis_3regions.region3,vector_N1_region3,'-sr');
        plot(all_time_range_bis_3regions.region3,vector_N2_region3,'-sb');
    elseif error_compute == 1
        errorbar(all_time_range_bis_3regions.region1,vector_N1_region1,vector_N1_region1_std,'-or');
        hold all
        errorbar(all_time_range_bis_3regions.region1,vector_N2_region1,vector_N2_region1_std,'-ob');
        errorbar(all_time_range_bis_3regions.region2,vector_N1_region2,vector_N1_region2_std,'-xr');
        errorbar(all_time_range_bis_3regions.region2,vector_N2_region2,vector_N2_region2_std,'-xb');
        errorbar(all_time_range_bis_3regions.region3,vector_N1_region3,vector_N1_region3_std,'-sr');
        errorbar(all_time_range_bis_3regions.region3,vector_N2_region3,vector_N2_region3_std,'-sb');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        ylabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT count (a.u.)');
    legend('N1 0-45%','N2 0-45%','N1 45-70%','N2 45-70%','N1 70-100%','N2 70-100%')
    saveas(gcf,[save_stem 'Temporal_modulation_count_3regions.fig']);
    
    %----------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis_3regions.region1,vector_T1_region1,'-or');
        hold all
        plot(all_time_range_bis_3regions.region1,vector_T2_region1,'-ob');
        plot(all_time_range_bis_3regions.region2,vector_T1_region2,'-xr');
        plot(all_time_range_bis_3regions.region2,vector_T2_region2,'-xb');
        plot(all_time_range_bis_3regions.region3,vector_T1_region3,'-sr');
        plot(all_time_range_bis_3regions.region3,vector_T2_region3,'-sb');
    elseif error_compute == 1
        errorbar(all_time_range_bis_3regions.region1,vector_T1_region1,vector_T1_region1_std,'-or');
        hold all
        errorbar(all_time_range_bis_3regions.region1,vector_T2_region1,vector_T2_region1_std,'-ob');
        errorbar(all_time_range_bis_3regions.region2,vector_T1_region2,vector_T1_region2_std,'-xr');
        errorbar(all_time_range_bis_3regions.region2,vector_T2_region2,vector_T2_region2_std,'-xb');
        errorbar(all_time_range_bis_3regions.region3,vector_T1_region3,vector_T1_region3_std,'-sr');
        errorbar(all_time_range_bis_3regions.region3,vector_T2_region3,vector_T2_region3_std,'-sb');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT lifetime (s)');
    legend('T1 0-45%','T2 0-45%','T1 45-70%','T2 45-70%','T1 70-100%','T2 70-100%')
    saveas(gcf,[save_stem 'Temporal_modulation_lifetime_3regions.fig']);
    
    %-----------------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis_3regions.region1,vector_P1_region1,'-or');
        hold all
        plot(all_time_range_bis_3regions.region2,vector_P1_region2,'-xr');
        plot(all_time_range_bis_3regions.region3,vector_P1_region3,'-sr');
    elseif error_compute == 1
        errorbar(all_time_range_bis_3regions.region1,vector_P1_region1,vector_P1_region1_std,'-or');
        hold all
        errorbar(all_time_range_bis_3regions.region2,vector_P1_region2,vector_P1_region2_std,'-xr');
        errorbar(all_time_range_bis_3regions.region3,vector_P1_region3,vector_P1_region3_std,'-sr');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('Proportion (%)');
    legend('P1 0-45%','P1 45-70%','P1 70-100%')
    saveas(gcf,[save_stem 'Temporal_modulation_proportion_3regions.fig']);
    
    %---------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis_3regions.region1,vector_f1_region1,'-or');
        hold all
        plot(all_time_range_bis_3regions.region1,vector_f2_region1,'-ob');
        plot(all_time_range_bis_3regions.region2,vector_f1_region2,'-xr');
        plot(all_time_range_bis_3regions.region2,vector_f2_region2,'-xb');
        plot(all_time_range_bis_3regions.region3,vector_f1_region3,'-sr');
        plot(all_time_range_bis_3regions.region3,vector_f2_region3,'-sb');
    elseif error_compute == 1
        errorbar(all_time_range_bis_3regions.region1,vector_f1_region1,vector_f1_region1_std,'-or');
        hold all
        errorbar(all_time_range_bis_3regions.region1,vector_f2_region1,vector_f2_region1_std,'-ob');
        errorbar(all_time_range_bis_3regions.region2,vector_f1_region2,vector_f1_region2_std,'-xr');
        errorbar(all_time_range_bis_3regions.region2,vector_f2_region2,vector_f2_region2_std,'-xb');
        errorbar(all_time_range_bis_3regions.region3,vector_f1_region3,vector_f1_region3_std,'-sr');
        errorbar(all_time_range_bis_3regions.region3,vector_f2_region3,vector_f2_region3_std,'-sb');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT frequency (/s)');
    legend('f1 0-45%','f2 0-45%','f1 45-70%','f2 45-70%','f1 70-100%','f2 70-100%')
    saveas(gcf,[save_stem 'Temporal_modulation_frequency_3regions.fig']);
    
    %---------------------------
    
    figure,
    if error_compute == 0
        plot(all_time_range_bis_3regions.region1,vector_d1_region1,'-or');
        hold all
        plot(all_time_range_bis_3regions.region1,vector_d2_region1,'-ob');
        plot(all_time_range_bis_3regions.region2,vector_d1_region2,'-xr');
        plot(all_time_range_bis_3regions.region2,vector_d2_region2,'-xb');
        plot(all_time_range_bis_3regions.region3,vector_d1_region3,'-sr');
        plot(all_time_range_bis_3regions.region3,vector_d2_region3,'-sb');
    elseif error_compute == 1
        errorbar(all_time_range_bis_3regions.region1,vector_d1_region1,vector_d1_region1_std,'-or');
        hold all
        errorbar(all_time_range_bis_3regions.region1,vector_d2_region1,vector_d2_region1_std,'-ob');
        errorbar(all_time_range_bis_3regions.region2,vector_d1_region2,vector_d1_region2_std,'-xr');
        errorbar(all_time_range_bis_3regions.region2,vector_d2_region2,vector_d2_region2_std,'-xb');
        errorbar(all_time_range_bis_3regions.region3,vector_d1_region3,vector_d1_region3_std,'-sr');
        errorbar(all_time_range_bis_3regions.region3,vector_d2_region3,vector_d2_region3_std,'-sb');
    end
    if time_reference_choice == 0
        xlabel('Time from anaphase onset (s)');
    elseif time_reference_choice == 1
        xlabel(sprintf('Time from NEBD (s)'));
    end
    ylabel('MT density (/min/um2)');
    legend('d1 0-45%','d2 0-45%','d1 45-70%','d2 45-70%','d1 70-100%','d2 70-100%')
    saveas(gcf,[save_stem 'Temporal_modulation_density_3regions.fig']);
    
    %-------------------------
    % display also data with ref time = nebd (still using onset of furrow
    % ingression as reference at the cortex)
    
    if additionnal_ref_nebd == 1
        
        figure,
        if error_compute == 0
            plot(all_time_range_bis_3regions_nebd.region1,vector_N1_region1,'-or');
            hold all
            plot(all_time_range_bis_3regions_nebd.region1,vector_N2_region1,'-ob');
            plot(all_time_range_bis_3regions_nebd.region2,vector_N1_region2,'-xr');
            plot(all_time_range_bis_3regions_nebd.region2,vector_N2_region2,'-xb');
            plot(all_time_range_bis_3regions_nebd.region3,vector_N1_region3,'-sr');
            plot(all_time_range_bis_3regions_nebd.region3,vector_N2_region3,'-sb');
        elseif error_compute == 1
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_N1_region1,vector_N1_region1_std,'-or');
            hold all
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_N2_region1,vector_N2_region1_std,'-ob');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_N1_region2,vector_N1_region2_std,'-xr');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_N2_region2,vector_N2_region2_std,'-xb');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_N1_region3,vector_N1_region3_std,'-sr');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_N2_region3,vector_N2_region3_std,'-sb');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT count (a.u.)');
        legend('N1 0-45%','N2 0-45%','N1 45-70%','N2 45-70%','N1 70-100%','N2 70-100%')
        saveas(gcf,[save_stem 'Temporal_modulation_count_3regions_nebd.fig']);
        
        %-------------------
        
        figure,
        if error_compute == 0
            plot(all_time_range_bis_3regions_nebd.region1,vector_T1_region1,'-or');
            hold all
            plot(all_time_range_bis_3regions_nebd.region1,vector_T2_region1,'-ob');
            plot(all_time_range_bis_3regions_nebd.region2,vector_T1_region2,'-xr');
            plot(all_time_range_bis_3regions_nebd.region2,vector_T2_region2,'-xb');
            plot(all_time_range_bis_3regions_nebd.region3,vector_T1_region3,'-sr');
            plot(all_time_range_bis_3regions_nebd.region3,vector_T2_region3,'-sb');
        elseif error_compute == 1
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_T1_region1,vector_T1_region1_std,'-or');
            hold all
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_T2_region1,vector_T2_region1_std,'-ob');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_T1_region2,vector_T1_region2_std,'-xr');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_T2_region2,vector_T2_region2_std,'-xb');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_T1_region3,vector_T1_region3_std,'-sr');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_T2_region3,vector_T2_region3_std,'-sb');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT lifetime (s)');
        legend('T1 0-45%','T2 0-45%','T1 45-70%','T2 45-70%','T1 70-100%','T2 70-100%')
        saveas(gcf,[save_stem 'Temporal_modulation_lifetime_3regions_nebd.fig']);
        
        %------------------------
        
        figure,
        if error_compute == 0
            plot(all_time_range_bis_3regions_nebd.region1,vector_P1_region1,'-or');
            hold all
            plot(all_time_range_bis_3regions_nebd.region2,vector_P1_region2,'-xr');
            plot(all_time_range_bis_3regions_nebd.region3,vector_P1_region3,'-sr');
        elseif error_compute == 1
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_P1_region1,vector_P1_region1_std,'-or');
            hold all
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_P1_region2,vector_P1_region2_std,'-xr');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_P1_region3,vector_P1_region3_std,'-sr');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('Proportion (%)');
        legend('P1 0-45%','P1 45-70%','P1 70-100%')
        saveas(gcf,[save_stem 'Temporal_modulation_proportion_3regions_nebd.fig']);
        
        %------------------------
        
        figure,
        if error_compute == 0
            plot(all_time_range_bis_3regions_nebd.region1,vector_f1_region1,'-or');
            hold all
            plot(all_time_range_bis_3regions_nebd.region1,vector_f2_region1,'-ob');
            plot(all_time_range_bis_3regions_nebd.region2,vector_f1_region2,'-xr');
            plot(all_time_range_bis_3regions_nebd.region2,vector_f2_region2,'-xb');
            plot(all_time_range_bis_3regions_nebd.region3,vector_f1_region3,'-sr');
            plot(all_time_range_bis_3regions_nebd.region3,vector_f2_region3,'-sb');
        elseif error_compute == 1
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_f1_region1,vector_f1_region1_std,'-or');
            hold all
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_f2_region1,vector_f2_region1_std,'-ob');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_f1_region2,vector_f1_region2_std,'-xr');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_f2_region2,vector_f2_region2_std,'-xb');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_f1_region3,vector_f1_region3_std,'-sr');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_f2_region3,vector_f2_region3_std,'-sb');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT frequency (/s)');
        legend('f1 0-45%','f2 0-45%','f1 45-70%','f2 45-70%','f1 70-100%','f2 70-100%')
        saveas(gcf,[save_stem 'Temporal_modulation_frequency_3regions_nebd.fig']);
        
        %-------------------------
        
        figure,
        if error_compute == 0
            plot(all_time_range_bis_3regions_nebd.region1,vector_d1_region1,'-or');
            hold all
            plot(all_time_range_bis_3regions_nebd.region1,vector_d2_region1,'-ob');
            plot(all_time_range_bis_3regions_nebd.region2,vector_d1_region2,'-xr');
            plot(all_time_range_bis_3regions_nebd.region2,vector_d2_region2,'-xb');
            plot(all_time_range_bis_3regions_nebd.region3,vector_d1_region3,'-sr');
            plot(all_time_range_bis_3regions_nebd.region3,vector_d2_region3,'-sb');
        elseif error_compute == 1
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_d1_region1,vector_d1_region1_std,'-or');
            hold all
            errorbar(all_time_range_bis_3regions_nebd.region1,vector_d2_region1,vector_d2_region1_std,'-ob');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_d1_region2,vector_d1_region2_std,'-xr');
            errorbar(all_time_range_bis_3regions_nebd.region2,vector_d2_region2,vector_d2_region2_std,'-xb');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_d1_region3,vector_d1_region3_std,'-sr');
            errorbar(all_time_range_bis_3regions_nebd.region3,vector_d2_region3,vector_d2_region3_std,'-sb');
        end
        xlabel(sprintf('Time from NEBD (s)'));
        ylabel('MT density (a.u.)');
        legend('d1 0-45%','d2 0-45%','d1 45-70%','d2 45-70%','d1 70-100%','d2 70-100%')
        saveas(gcf,[save_stem 'Temporal_modulation_density_3regions_nebd.fig']);
        
    end
    
end

%% final results

if compute_wholeEmbryo == 1
    results_temporalspatial_modulation.entireEmbryo.lifetime.T1 = vector_T1;
    results_temporalspatial_modulation.entireEmbryo.lifetime.T2 = vector_T2;
    results_temporalspatial_modulation.entireEmbryo.proportion.P1 = vector_P1;
    results_temporalspatial_modulation.entireEmbryo.proportion.P2 = vector_P2;
    results_temporalspatial_modulation.entireEmbryo.proportion.N1 = vector_N1;
    results_temporalspatial_modulation.entireEmbryo.proportion.N2 = vector_N2;
    results_temporalspatial_modulation.entireEmbryo.proportion.f1 = vector_f1;
    results_temporalspatial_modulation.entireEmbryo.proportion.f2 = vector_f2;
    results_temporalspatial_modulation.entireEmbryo.proportion.d1 = vector_d1;
    results_temporalspatial_modulation.entireEmbryo.proportion.d2 = vector_d2;
    
    if error_compute == 1
        results_temporalspatial_modulation.entireEmbryo.lifetime.T1_std = vector_T1_std;
        results_temporalspatial_modulation.entireEmbryo.lifetime.T2_std = vector_T2_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.P1_std = vector_P1_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.P2_std = vector_P2_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.N1_std = vector_N1_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.N2_std = vector_N2_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.d1_std = vector_d1_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.d2_std = vector_d2_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.f1_std = vector_f1_std;
        results_temporalspatial_modulation.entireEmbryo.proportion.f2_std = vector_f2_std;             
    end
end

if compute_3regions == 1
    results_temporalspatial_modulation.region1.lifetime.T1 = vector_T1_region1;
    results_temporalspatial_modulation.region2.lifetime.T1 = vector_T1_region2;
    results_temporalspatial_modulation.region3.lifetime.T1 = vector_T1_region3;
    results_temporalspatial_modulation.region1.lifetime.T2 = vector_T2_region1;
    results_temporalspatial_modulation.region2.lifetime.T2 = vector_T2_region2;
    results_temporalspatial_modulation.region3.lifetime.T2 = vector_T2_region3;
    
    results_temporalspatial_modulation.region1.proportion.P1 = vector_P1_region1;
    results_temporalspatial_modulation.region2.proportion.P1 = vector_P1_region2;
    results_temporalspatial_modulation.region3.proportion.P1 = vector_P1_region3;
    results_temporalspatial_modulation.region1.proportion.P2 = vector_P2_region1;
    results_temporalspatial_modulation.region2.proportion.P2 = vector_P2_region2;
    results_temporalspatial_modulation.region3.proportion.P2 = vector_P2_region3;
    
    results_temporalspatial_modulation.region1.proportion.N1 = vector_N1_region1;
    results_temporalspatial_modulation.region2.proportion.N1 = vector_N1_region2;
    results_temporalspatial_modulation.region3.proportion.N1 = vector_N1_region3;
    results_temporalspatial_modulation.region1.proportion.N2 = vector_N2_region1;
    results_temporalspatial_modulation.region2.proportion.N2 = vector_N2_region2;
    results_temporalspatial_modulation.region3.proportion.N2 = vector_N2_region3;
    
    results_temporalspatial_modulation.region1.proportion.f1 = vector_f1_region1;
    results_temporalspatial_modulation.region2.proportion.f1 = vector_f1_region2;
    results_temporalspatial_modulation.region3.proportion.f1 = vector_f1_region3;
    results_temporalspatial_modulation.region1.proportion.f2 = vector_f2_region1;
    results_temporalspatial_modulation.region2.proportion.f2 = vector_f2_region2;
    results_temporalspatial_modulation.region3.proportion.f2 = vector_f2_region3;
    
    results_temporalspatial_modulation.region1.proportion.d1 = vector_d1_region1;
    results_temporalspatial_modulation.region2.proportion.d1 = vector_d1_region2;
    results_temporalspatial_modulation.region3.proportion.d1 = vector_d1_region3;
    results_temporalspatial_modulation.region1.proportion.d2 = vector_d2_region1;
    results_temporalspatial_modulation.region2.proportion.d2 = vector_d2_region2;
    results_temporalspatial_modulation.region3.proportion.d2 = vector_d2_region3;
    
    if error_compute == 1
        results_temporalspatial_modulation.region1.lifetime.T1_std = vector_T1_region1_std;
        results_temporalspatial_modulation.region2.lifetime.T1_std = vector_T1_region2_std;
        results_temporalspatial_modulation.region3.lifetime.T1_std = vector_T1_region3_std;
        results_temporalspatial_modulation.region1.lifetime.T2_std = vector_T2_region1_std;
        results_temporalspatial_modulation.region2.lifetime.T2_std = vector_T2_region2_std;
        results_temporalspatial_modulation.region3.lifetime.T2_std = vector_T2_region3_std;
        
        results_temporalspatial_modulation.region1.proportion.P1_std = vector_P1_region1_std;
        results_temporalspatial_modulation.region2.proportion.P1_std = vector_P1_region2_std;
        results_temporalspatial_modulation.region3.proportion.P1_std = vector_P1_region3_std;
        results_temporalspatial_modulation.region1.proportion.P2_std = vector_P2_region1_std;
        results_temporalspatial_modulation.region2.proportion.P2_std = vector_P2_region2_std;
        results_temporalspatial_modulation.region3.proportion.P2_std = vector_P2_region3_std;
        
        results_temporalspatial_modulation.region1.proportion.N1_std = vector_N1_region1_std;
        results_temporalspatial_modulation.region2.proportion.N1_std = vector_N1_region2_std;
        results_temporalspatial_modulation.region3.proportion.N1_std = vector_N1_region3_std;
        results_temporalspatial_modulation.region1.proportion.N2_std = vector_N2_region1_std;
        results_temporalspatial_modulation.region2.proportion.N2_std = vector_N2_region2_std;
        results_temporalspatial_modulation.region3.proportion.N2_std = vector_N2_region3_std;
        
        results_temporalspatial_modulation.region1.proportion.d1_std = vector_d1_region1_std;
        results_temporalspatial_modulation.region2.proportion.d1_std = vector_d1_region2_std;
        results_temporalspatial_modulation.region3.proportion.d1_std = vector_d1_region3_std;
        results_temporalspatial_modulation.region1.proportion.d2_std = vector_d2_region1_std;
        results_temporalspatial_modulation.region2.proportion.d2_std = vector_d2_region2_std;
        results_temporalspatial_modulation.region3.proportion.d2_std = vector_d2_region3_std;
        
        results_temporalspatial_modulation.region1.proportion.f1_std = vector_f1_region1_std;
        results_temporalspatial_modulation.region2.proportion.f1_std = vector_f1_region2_std;
        results_temporalspatial_modulation.region3.proportion.f1_std = vector_f1_region3_std;
        results_temporalspatial_modulation.region1.proportion.f2_std = vector_f2_region1_std;
        results_temporalspatial_modulation.region2.proportion.f2_std = vector_f2_region2_std;
        results_temporalspatial_modulation.region3.proportion.f2_std = vector_f2_region3_std;
    end
end

name = strcat('Results_temporalspatial_modulation_' , xmlfile_bkp, '.mat');
save(fullfile(main_path,name), '-struct','results_temporalspatial_modulation');


%% clear

clear all
close all


end

