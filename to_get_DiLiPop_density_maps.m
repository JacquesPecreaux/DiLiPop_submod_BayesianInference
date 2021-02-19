function to_get_DiLiPop_density_maps( xmlfile,save_stem,folder_tag,compute_wholeEmbryo,compute_3regions,time_reference_choice,...
    minLength_tracks,limit_nb_tracks_for_fitting)

% this function enables to generate :
% 1/ the number of contacts per frame that are assigned to short-lived or
% long-lived populations
% 2/ the DiLiPop density maps, i.e. the density of MT contacts assigned to short-lived and long-lived populations
% for each embryo of the studied xmlfile, and the averaged density maps for the given condition

% The assignment can be done , 
% - either considering the whole embryo to measure the two distinct dynamical behaviours, 
% - or distinguishing the three cortical regions (0-40% of AP axis, 40-70% and 70-100%) of C. elegans embryo

% OUTPUT FILES of interest
% - mat file: several tracks_duration_histo... that contains the duration distributions of all the embryos 
% for each period (metaphase/anaphase)/region investigated (whole, region1, region2, region3
% - mat file: final_results-BayesianInference_sum_mle that contains the result of the fitting of the duration distributions 
% in each region/period couple
% - mat file: familyAssignement that has eachtrack assigned to short-lived or long-lived populations
% - mat file: embryoSelection that contains the counters per embryo for each assigned population
% - mat file: averageCountersWholeEmbryo which contains the data of the averaged count of assigned contcat along mitosis.
% - mat file: averagecounters_nbR10Regions-area which contains the data of the density DiLipop maps
% - averaged plots of the number of contacts along time in each assigned
% population and also without assignment (with individual data surperosed or not)
% - plot of Averaged DiLiPop density map in 3D, or in 2D with density encoded in color.
% - plot of individual DiLipop density in 2D
% NB: at the end of the script, script will ask user to give maximal
% density level and time period to be displayed for the density maps

clear global

global general_param;
global param;


%% input args

% 1st: xmlfile that contains information regarding :
% - the studied condition (general_param.),
% - each sample (e.g. embryo) composing this condition (param.)
% 2nd: path to save the results of the analysis
% 3rd: tag of the folder: the data needed for this function are saved in the main folder "data_dir", and in the subfolder 
% corresponded to the project (PID), sub-subfolder correponddnt to the embyro (ID), and in the sub-sub-subfolder identified 
% with the given tag. This tag can be empty if no tag used.
% 4th: to choose whether the study of two distinct dynamical behaviors is done considering the whole embryo. 
% If yes, choose the period duration to resolve temporally the parameters. It requires a minimal number of tracks to 
% perform the statistical analysis, this will the limiting parameter that sets the window size.
% 5th: to choose whether the study of two distinct dynamical behaviors is done considering the three cortical regions. 
% If yes, choose the period duration to resolve temporally the parameters. It requires a minimal number of tracks to 
% perform the statistical analysis, this will the limiting parameter that sets the window size.
% 6th: reference time used to align the embryos : the two time reference possible are 1/ the onset of the furrow ingression 
% at about mid-anaphase, and 2/ the pseudo cleavage end happeening before pronuclei meeting
% 7th: minimal duration of a track considered in the analysis
% 8th: minimal number of tracks required to perform correctly the DiLiPop statistical analysis


%% files that will be downloaded and required for the analysis

% These data  are saved in the main folder "data_dir", and in the subfolder corresponded to the project (PID), 
% sub-subfolder correpondant to the embyro (ID), and in the sub-sub-subfolder identified with the given tag.
% Files neded are:
% - region Area: area of the embryo along mitosis: whole, or divided alon AP axis, will be needed to compute the density
% - regionXlimit: coordinates of the embryo limits in the full frame along AP axis, 
% --> will be needed to assign the contact to the different regions of the embryo
% - regionXlength: width of the sub-regions of the embryo: = full embryo length along AP axis divided by 
% the number of regions (here, 10 regiosn to generate the denisty maps).
% --> will be needed to assign the contact to the different regions of the embryo
% - furrow_position_convexity: gives the furrow ingression onset timing
% furrow_onset in sec = furrow_position_convexity.image_start_detection/param.sp6 -> in second from the start of the movie
% - dataTracks_rotated: gives the tracks: duration, coordinates of the position, first frame the track appears, 
% end frame of the tracks.



%% process inputs and ask for them

if nargin<1
    [xmlfile, p] = uigetfile('*.mat','Please choose a job file to process');
    [~,xmlfile_bkp,~] = fileparts(xmlfile);
    xmlfile = fullfile(p, xmlfile);
end
if nargin<2
    [p,f] = fileparts(xmlfile);
    [save_stem, p] = uiputfile('*','Please provide a path and stem name for saving figures',fullfile(p,[f 'DiLiPop_density_maps_']));
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

if nargin < 3
    folder_tag = input_perso(['Set the tag used to get proper treatment at the cortex folder: '],'');
end

visualize_individual_embryo = 1;
models = {'MonoExpo' 'DoubleExpo' };

if nargin < 4
    compute_wholeEmbryo  = input_perso([sprintf('Do you wish to cmpute parameters in whole embryo? (yes = 1, No = 0 ')], 1 );
    if compute_wholeEmbryo == 1
        window_size_metaphase_sec  = ...
            input_perso([sprintf('Please give window size during metaphase in sec for whole embryo (default value = ')], 30 );
        
        window_size_anaphase_sec = ...
            input_perso([sprintf('Please give window size during anaphase in sec  for whole embryo (default value = ')], 30 );
    end
end

if nargin < 5
    compute_3regions  = input_perso([sprintf('Do you wish to cmpute parameters in 3 regions? (yes = 1, No = 0 ')], 1 );
    if compute_3regions == 1
        window_size_metaphase_sec_3regions  = ...
            input_perso([sprintf('Please give window size during metaphase in sec for 3 regions (default value = ')], 60 );
        
        window_size_anaphase_sec_3regions = ...
            input_perso([sprintf('Please give window size during anaphase in sec  for 3 regions (default value = ')], 60 );
    end
end

if nargin < 6
    time_reference_choice = input_perso(['Which reference time to use? (0= furrow ingression onset, 1 = pseudo-cleavage end'], 0);
end


%% load parameters from the downloaded strcuture

load(xmlfile); % load structure containing general_param and param
general_param = saveVarsMat_new.general_param;% read general_param from structure

if nargin < 7
    minLength_tracks = general_param.cortex_analysis.minLength;
end

if nargin < 8
    limit_nb_tracks_for_fitting = general_param.cortex_analysis.threshold_MTnumber_forMLEfit;
end


%% treat individual embryo and get their data

nbEmbryo = 0;
nbEmbryo_3regions = 0;
max_period_metaphase = 0;
max_period_anaphase = 0;
max_period_metaphase_3regions = 0;
max_period_anaphase_3regions = 0;
nbEmbryo_period_metaphase = zeros(1,50);
nbEmbryo_period_anaphase = zeros(1,50);
nbEmbryo_period_metaphase_above250 = zeros(1,50);
nbEmbryo_period_anaphase_above250 = zeros(1,50);
nbEmbryo_period_metaphase_3regions = zeros(1,50);
nbEmbryo_period_anaphase_3regions = zeros(1,50);
nbEmbryo_period_metaphase_above250_3regions.region1 = zeros(1,50);
nbEmbryo_period_anaphase_above250_3regions.region1 = zeros(1,50);
nbEmbryo_period_metaphase_above250_3regions.region2 = zeros(1,50);
nbEmbryo_period_anaphase_above250_3regions.region2 = zeros(1,50);
nbEmbryo_period_metaphase_above250_3regions.region3 = zeros(1,50);
nbEmbryo_period_anaphase_above250_3regions.region3 = zeros(1,50);

for i = 1 : 50
    if compute_wholeEmbryo == 1
        extension1 = ['meta_minus' num2str(i)];
        embryo_index_metaphase.(extension1).entireEmbryo = [];
        extension2 = ['ana_plus' num2str(i)];
        embryo_index_anaphase.(extension2).entireEmbryo = [];
    end
    if compute_3regions == 1
        extension1 = ['meta_minus' num2str(i) '_3regions'];
        embryo_index_metaphase.(extension1).region1 = [];
        embryo_index_metaphase.(extension1).region2 = [];
        embryo_index_metaphase.(extension1).region3 = [];
        extension2 = ['ana_plus' num2str(i) '_3regions'];
        embryo_index_anaphase.(extension2).region1 = [];
        embryo_index_anaphase.(extension2).region2 = [];
        embryo_index_anaphase.(extension2).region3 = [];
    end
end

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp(param.sp1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DOWNLOAD PART
        % NB: data saved in two ways: with or without the use of param.extra
        
        %-------------------------
        % download Xlimit, Xlength, Area
        mainDirectory = 'contour detection';
        pathMainDirectory3 = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        clear mainDirectory
        nameData2 = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat'];
        nameData = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', short_name) '.mat'];
        if exist(nameData2,'file') == 2
            nameData_ = [sprintf('%s%s%s',pathMainDirectory3,'regionArea-', param.extra, short_name) '.mat'];
            regionArea = load(nameData_); % in pixels**2
            nameData2_ = [sprintf('%s%s%s',pathMainDirectory3,'regionXlimit-', param.extra, short_name) '.mat'];
            regionXlimit = load(nameData2_);
            nameData3_ = [sprintf('%s%s%s',pathMainDirectory3,'regionXlength-', param.extra, short_name) '.mat'];
            regionXlength = load(nameData3_);
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
        pathMainDirectory2 = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        clear mainDirectory
        filename0 = ['furrow_position_convexity-', short_name, '.mat'];
        filename1 = ['furrow_position_convexity-', short_name, param.extra, '.mat'];
        if exist(fullfile(pathMainDirectory2,filename1),'file') == 2
            furrow_position = load(fullfile(pathMainDirectory2,filename1));
            clear filename1
        elseif exist(fullfile(pathMainDirectory2,filename0),'file') == 2
            furrow_position = load(fullfile(pathMainDirectory2,filename0));
            clear filename0
        else
            disp('Error: furrow position file not found');
        end
        
        %------------------------------
        % download dataTracks_rotated
        mainDirectory = strcat('analysis at the cortex', folder_tag);
        pathMainDirectory = strcat(param.basepath , '/' , param.sp1 , '/', mainDirectory, '/');
        nameData1 = [sprintf('%s%s%s%s',pathMainDirectory,'dataTracks_rotated-', short_name) , param.extra,'.mat'];
        nameData2 = [sprintf('%s%s%s%s',pathMainDirectory,'dataTracks_rotated-', short_name) , '.mat'];
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
            clear filename1
        elseif exist(nameData2,'file') == 2
            filename2 = ['dataTracks_rotated-', short_name, '.mat'];
            dataTracks_rotated = load(fullfile(pathMainDirectory,filename2));
            clear filename2
        else
            disp('Error: datatracks file not found');
        end
        % to correct when saving of dataTracks done as dataTracks_rotated.dataTracks_rotated.entireEmbryo      
        if  isfield(dataTracks_rotated,'dataTracks_rotated') == 1 % account for data not save properly.
            dataTracks_rotated_.entireEmbryo = dataTracks_rotated.dataTracks_rotated.entireEmbryo;
            dataTracks_rotated_.numTimePoints = dataTracks_rotated.dataTracks_rotated.numTimePoints;
            clear dataTracks_rotated
            dataTracks_rotated = dataTracks_rotated_;
            clear dataTracks_rotated_
        end
        % to correct when saving of dataTracks done as
        % entireEmbryo.entireRecording.tracksX, while expected as
        % entireEmbryo.tarcksX
        if  isfield(dataTracks_rotated.entireEmbryo,'entireRecording') == 1 % account for data not save properly.
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
            
            limit1 = general_param.cortex_analysis.analysis_3SpatialRegions_limit1; % 40% of AP axis
            limit2 = general_param.cortex_analysis.analysis_3SpatialRegions_limit2; % 70% of AP axis
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
            image_reference_m2a = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6;
            image_reference_nebd = furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset + param.delta_early_metaphase)*param.sp6;
        elseif time_reference_choice == 1
            image_reference_m2a =  (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) * param.sp6;
            image_reference_nebd =  (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage ) * param.sp6;
        end
        nbEmbryo = nbEmbryo + 1;
        
        % whole embryo windowing
        if compute_wholeEmbryo == 1
            
            window_size_metaphase = window_size_metaphase_sec * param.sp6; % in image nb
            window_size_anaphase = window_size_anaphase_sec * param.sp6; % in image nb
            % get reference time in metaphase
            image_reference_metaphase = [ image_reference_m2a ];
            for ii = 1 : 50
                tmp = image_reference_m2a - ii * window_size_metaphase;
                if tmp >  0
                    image_reference_metaphase = [ [image_reference_metaphase ] tmp ];
                else
                    break
                end
            end
            image_reference_metaphase = [ [image_reference_metaphase ] 1 ]; % param.sp2
            
            for jj = 1 : length(image_reference_metaphase)
                nbEmbryo_period_metaphase(1,jj) = nbEmbryo_period_metaphase(1,jj) + 1;
            end
            if max_period_metaphase < length(image_reference_metaphase)
                max_period_metaphase = length(image_reference_metaphase);
            end
            
            % get reference time in anaphase
            image_reference_anaphase = [image_reference_m2a];
            for ii = 1 : 50
                tmp = image_reference_m2a + ii * window_size_anaphase;
                if tmp <  param.sp3
                    image_reference_anaphase = [ [image_reference_anaphase ] tmp ];
                else
                    break
                end
            end
            image_reference_anaphase = [ [image_reference_anaphase ] param.sp3 ];
            
            for jj = 1 : length(image_reference_anaphase)
                nbEmbryo_period_anaphase(1,jj) = nbEmbryo_period_anaphase(1,jj) + 1;
            end
            
            if max_period_anaphase < length(image_reference_anaphase)
                max_period_anaphase = length(image_reference_anaphase);
            end
        end
        
        % regions
        if compute_3regions == 1
            
            window_size_metaphase_3regions = window_size_metaphase_sec_3regions * param.sp6; % in image nb
            window_size_anaphase_3regions = window_size_anaphase_sec_3regions * param.sp6; % in image nb
            
            % get reference time in metaphase
            image_reference_metaphase_3regions = [ image_reference_m2a ];
            for ii = 1 : 50
                tmp = image_reference_m2a - ii * window_size_metaphase_3regions;
                if tmp >  0
                    image_reference_metaphase_3regions = [ [image_reference_metaphase_3regions ] tmp ];
                else
                    break
                end
            end
            image_reference_metaphase_3regions = [ [image_reference_metaphase_3regions ] 1 ]; % param.sp2
            
            for jj = 1 : length(image_reference_metaphase_3regions)
                nbEmbryo_period_metaphase_3regions(1,jj) = nbEmbryo_period_metaphase_3regions(1,jj) + 1;
            end
            if max_period_metaphase_3regions < length(image_reference_metaphase_3regions)
                max_period_metaphase_3regions = length(image_reference_metaphase_3regions);
            end
            
            % get reference time in anaphase
            image_reference_anaphase_3regions = [image_reference_m2a];
            for ii = 1 : 50
                tmp = image_reference_m2a + ii * window_size_anaphase_3regions;
                if tmp <  param.sp3
                    image_reference_anaphase_3regions = [ [image_reference_anaphase_3regions ] tmp ];
                else
                    break
                end
            end
            image_reference_anaphase_3regions = [ [image_reference_anaphase_3regions ] param.sp3 ];
            
            for jj = 1 : length(image_reference_anaphase_3regions)
                nbEmbryo_period_anaphase_3regions(1,jj) = nbEmbryo_period_anaphase_3regions(1,jj) + 1;
            end
            
            if max_period_anaphase_3regions < length(image_reference_anaphase_3regions)
                max_period_anaphase_3regions = length(image_reference_anaphase_3regions);
            end
            
        end
        
        %--------------
        % separate tracks in temporal windows
        
        % whole embryo
        if compute_wholeEmbryo == 1
            
            for iPeriod_metaphase = 1 : length(image_reference_metaphase)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase)];
                entireEmbryo.(extension).numTracks = 0;
                entireEmbryo.(extension).lengthTracks = [];
                entireEmbryo.(extension).index = [];
                entireEmbryo.(extension).indexXStart = [];
                entireEmbryo.(extension).indexXEnd = [];
            end
            
            for iPeriod_anaphase = 1 : length(image_reference_anaphase)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase)];
                entireEmbryo.(extension).numTracks = 0;
                entireEmbryo.(extension).lengthTracks = [];
                entireEmbryo.(extension).index = [];
                entireEmbryo.(extension).indexXStart = [];
                entireEmbryo.(extension).indexXEnd = [];
            end
            
            for iTrack = 1 : dataTracks_rotated.entireEmbryo.numTracks
                for iPeriod_metaphase = 1 : length(image_reference_metaphase)-1
                    extension = ['meta_minus' num2str(iPeriod_metaphase)];
                    
                    if image_reference_metaphase(iPeriod_metaphase+1) <= dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) < image_reference_metaphase(iPeriod_metaphase)
                        entireEmbryo.(extension).numTracks = entireEmbryo.(extension).numTracks +1;
                        entireEmbryo.(extension).lengthTracks(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.lengthTracks(iTrack);
                        entireEmbryo.(extension).index(entireEmbryo.(extension).numTracks) = iTrack;
                        entireEmbryo.(extension).indexXStart(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.indexXStart(iTrack);
                        entireEmbryo.(extension).indexXEnd(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.indexXEnd(iTrack);
                        entireEmbryo.(extension).tracksX(:,entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.tracksX(:,iTrack);
                        entireEmbryo.(extension).tracksY(:,entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.tracksY(:,iTrack);
                    end
                end
                for iPeriod_anaphase = 1 : length(image_reference_anaphase)-1
                    extension = ['ana_plus' num2str(iPeriod_anaphase)];
                    
                    if image_reference_anaphase(iPeriod_anaphase) <= dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,iTrack) < image_reference_anaphase(iPeriod_anaphase+1)
                        entireEmbryo.(extension).numTracks = entireEmbryo.(extension).numTracks +1;
                        entireEmbryo.(extension).lengthTracks(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.lengthTracks(iTrack);
                        entireEmbryo.(extension).index(entireEmbryo.(extension).numTracks) = iTrack;
                        entireEmbryo.(extension).indexXStart(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.indexXStart(iTrack);
                        entireEmbryo.(extension).indexXEnd(entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.indexXEnd(iTrack);
                        entireEmbryo.(extension).tracksX(:,entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.tracksX(:,iTrack);
                        entireEmbryo.(extension).tracksY(:,entireEmbryo.(extension).numTracks) = dataTracks_rotated.entireEmbryo.tracksY(:,iTrack);
                    end
                end
            end
        end
        
        % for 3 regions
        if compute_3regions == 1
            
            for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                region1.(extension).numTracks = 0;
                region1.(extension).lengthTracks = [];
                region1.(extension).index = [];
                region1.(extension).indexXStart = [];
                region1.(extension).indexXEnd = [];
                region2.(extension).numTracks = 0;
                region2.(extension).lengthTracks = [];
                region2.(extension).index = [];
                region2.(extension).indexXStart = [];
                region2.(extension).indexXEnd = [];
                region3.(extension).numTracks = 0;
                region3.(extension).lengthTracks = [];
                region3.(extension).index = [];
                region3.(extension).indexXStart = [];
                region3.(extension).indexXEnd = [];
            end
            
            for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                region1.(extension).numTracks = 0;
                region1.(extension).lengthTracks = [];
                region1.(extension).index = [];
                region1.(extension).indexXStart = [];
                region1.(extension).indexXEnd = [];
                region2.(extension).numTracks = 0;
                region2.(extension).lengthTracks = [];
                region2.(extension).index = [];
                region2.(extension).indexXStart = [];
                region2.(extension).indexXEnd = [];
                region3.(extension).numTracks = 0;
                region3.(extension).lengthTracks = [];
                region3.(extension).index = [];
                region3.(extension).indexXStart = [];
                region3.(extension).indexXEnd = [];
            end
            
            for iTrack = 1 : dataTracks_rotated.region1.numTracks
                for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                    extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                    if image_reference_metaphase_3regions(iPeriod_metaphase+1) <= dataTracks_rotated.region1.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region1.indexXStart(1,iTrack) < image_reference_metaphase_3regions(iPeriod_metaphase)
                        region1.(extension).numTracks = region1.(extension).numTracks +1;
                        region1.(extension).lengthTracks(region1.(extension).numTracks) = dataTracks_rotated.region1.lengthTracks(iTrack);
                        region1.(extension).index(region1.(extension).numTracks) = iTrack;
                        region1.(extension).indexXStart(region1.(extension).numTracks) = dataTracks_rotated.region1.indexXStart(iTrack);
                        region1.(extension).indexXEnd(region1.(extension).numTracks) = dataTracks_rotated.region1.indexXEnd(iTrack);
                        region1.(extension).tracksX(:,region1.(extension).numTracks) = dataTracks_rotated.region1.tracksX(:,iTrack);
                        region1.(extension).tracksY(:,region1.(extension).numTracks) = dataTracks_rotated.region1.tracksY(:,iTrack);
                    end
                end
                for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                    extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                    if image_reference_anaphase_3regions(iPeriod_anaphase) <= dataTracks_rotated.region1.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region1.indexXStart(1,iTrack) < image_reference_anaphase_3regions(iPeriod_anaphase+1)
                        region1.(extension).numTracks = region1.(extension).numTracks +1;
                        region1.(extension).lengthTracks(region1.(extension).numTracks) = dataTracks_rotated.region1.lengthTracks(iTrack);
                        region1.(extension).index(region1.(extension).numTracks) = iTrack;
                        region1.(extension).indexXStart(region1.(extension).numTracks) = dataTracks_rotated.region1.indexXStart(iTrack);
                        region1.(extension).indexXEnd(region1.(extension).numTracks) = dataTracks_rotated.region1.indexXEnd(iTrack);
                        region1.(extension).tracksX(:,region1.(extension).numTracks) = dataTracks_rotated.region1.tracksX(:,iTrack);
                        region1.(extension).tracksY(:,region1.(extension).numTracks) = dataTracks_rotated.region1.tracksY(:,iTrack);
                    end
                end
            end
            
            for iTrack = 1 : dataTracks_rotated.region2.numTracks
                for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                    extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                    if image_reference_metaphase_3regions(iPeriod_metaphase+1) <= dataTracks_rotated.region2.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region2.indexXStart(1,iTrack) < image_reference_metaphase_3regions(iPeriod_metaphase)
                        region2.(extension).numTracks = region2.(extension).numTracks +1;
                        region2.(extension).lengthTracks(region2.(extension).numTracks) = dataTracks_rotated.region2.lengthTracks(iTrack);
                        region2.(extension).index(region2.(extension).numTracks) = iTrack;
                        region2.(extension).indexXStart(region2.(extension).numTracks) = dataTracks_rotated.region2.indexXStart(iTrack);
                        region2.(extension).indexXEnd(region2.(extension).numTracks) = dataTracks_rotated.region2.indexXEnd(iTrack);
                        region2.(extension).tracksX(:,region2.(extension).numTracks) = dataTracks_rotated.region2.tracksX(:,iTrack);
                        region2.(extension).tracksY(:,region2.(extension).numTracks) = dataTracks_rotated.region2.tracksY(:,iTrack);
                    end
                end
                for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                    extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                    if image_reference_anaphase_3regions(iPeriod_anaphase) <= dataTracks_rotated.region2.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region2.indexXStart(1,iTrack) < image_reference_anaphase_3regions(iPeriod_anaphase+1)
                        region2.(extension).numTracks = region2.(extension).numTracks +1;
                        region2.(extension).lengthTracks(region2.(extension).numTracks) = dataTracks_rotated.region2.lengthTracks(iTrack);
                        region2.(extension).index(region2.(extension).numTracks) = iTrack;
                        region2.(extension).indexXStart(region2.(extension).numTracks) = dataTracks_rotated.region2.indexXStart(iTrack);
                        region2.(extension).indexXEnd(region2.(extension).numTracks) = dataTracks_rotated.region2.indexXEnd(iTrack);
                        region2.(extension).tracksX(:,region2.(extension).numTracks) = dataTracks_rotated.region2.tracksX(:,iTrack);
                        region2.(extension).tracksY(:,region2.(extension).numTracks) = dataTracks_rotated.region2.tracksY(:,iTrack);
                    end
                end
            end
            
            for iTrack = 1 : dataTracks_rotated.region3.numTracks
                for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                    extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                    if image_reference_metaphase_3regions(iPeriod_metaphase+1) <= dataTracks_rotated.region3.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region3.indexXStart(1,iTrack) < image_reference_metaphase_3regions(iPeriod_metaphase)
                        region3.(extension).numTracks = region3.(extension).numTracks +1;
                        region3.(extension).lengthTracks(region3.(extension).numTracks) = dataTracks_rotated.region3.lengthTracks(iTrack);
                        region3.(extension).index(region3.(extension).numTracks) = iTrack;
                        region3.(extension).indexXStart(region3.(extension).numTracks) = dataTracks_rotated.region3.indexXStart(iTrack);
                        region3.(extension).indexXEnd(region3.(extension).numTracks) = dataTracks_rotated.region3.indexXEnd(iTrack);
                        region3.(extension).tracksX(:,region3.(extension).numTracks) = dataTracks_rotated.region3.tracksX(:,iTrack);
                        region3.(extension).tracksY(:,region3.(extension).numTracks) = dataTracks_rotated.region3.tracksY(:,iTrack);
                    end
                end
                for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                    extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                    if image_reference_anaphase_3regions(iPeriod_anaphase) <= dataTracks_rotated.region3.indexXStart(1,iTrack) ...
                            && dataTracks_rotated.region3.indexXStart(1,iTrack) < image_reference_anaphase_3regions(iPeriod_anaphase+1)
                        region3.(extension).numTracks = region3.(extension).numTracks +1;
                        region3.(extension).lengthTracks(region3.(extension).numTracks) = dataTracks_rotated.region3.lengthTracks(iTrack);
                        region3.(extension).index(region3.(extension).numTracks) = iTrack;
                        region3.(extension).indexXStart(region3.(extension).numTracks) = dataTracks_rotated.region3.indexXStart(iTrack);
                        region3.(extension).indexXEnd(region3.(extension).numTracks) = dataTracks_rotated.region3.indexXEnd(iTrack);
                        region3.(extension).tracksX(:,region3.(extension).numTracks) = dataTracks_rotated.region3.tracksX(:,iTrack);
                        region3.(extension).tracksY(:,region3.(extension).numTracks) = dataTracks_rotated.region3.tracksY(:,iTrack);
                    end
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get areas in each temporal window and regions
        
        
        area_withRepartition.raw.entireEmbryo = regionArea.entireEmbryo.nbR1;
        xStart_withRepartition.raw.entireEmbryo = regionXlimit.entireEmbryo.nbR1;
        deltaX_withRepartition.raw.entireEmbryo = regionXlength.entireEmbryo.nbR1;
        area_withRepartition.raw.nbR10 = regionArea.entireEmbryo.nbR10;
        xStart_withRepartition.raw.nbR10 = regionXlimit.entireEmbryo.nbR10;
        deltaX_withRepartition.raw.nbR10 = regionXlength.entireEmbryo.nbR10;
        
        if compute_3regions == 1
            area_withRepartition.raw.region1 =  area_3regions.region1;
            xStart_withRepartition.raw.region1 = xStart_3regions.region1;
            deltaX_withRepartition.raw.region1 = xEnd_3regions.region1;
            area_withRepartition.raw.region2 =  area_3regions.region2;
            xStart_withRepartition.raw.region2 = xStart_3regions.region2;
            deltaX_withRepartition.raw.region2 = xEnd_3regions.region2;
            area_withRepartition.raw.region3 =  area_3regions.region3;
            xStart_withRepartition.raw.region3 = xStart_3regions.region3;
            deltaX_withRepartition.raw.region3 = xEnd_3regions.region3;
        end
        
        % save timing furrow onset in this strcuture
        area_withRepartition.furrow_onset = (furrow_position.image_start_detection - param.sp2)/param.sp6; % in sec with image_detection refrence which is at start recording
        
        % whole embryo
        if compute_wholeEmbryo == 1
            for iPeriod_metaphase = 1 : length(image_reference_metaphase)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase)];
                area_withRepartition.entireEmbryo.(extension) = ...
                    mean( regionArea.entireEmbryo.nbR1(image_reference_metaphase(iPeriod_metaphase+1): image_reference_metaphase(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
            end
            
            for iPeriod_anaphase = 1 : length(image_reference_anaphase)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase)];
                if ( image_reference_anaphase(iPeriod_anaphase)< 0 || image_reference_anaphase(iPeriod_anaphase)> 0  ) && image_reference_anaphase(iPeriod_anaphase+1)> 0
                    if image_reference_anaphase(iPeriod_anaphase) < 0
                        image_reference_anaphase(iPeriod_anaphase) = 1;
                    end
                    area_withRepartition.entireEmbryo.(extension) = ...
                        mean( regionArea.entireEmbryo.nbR1(image_reference_anaphase(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                end
            end
        end
        
        % 3 regions
        if compute_3regions == 1
            for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                area_withRepartition.region1.(extension)= ...
                    mean( area_3regions.region1(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
                area_withRepartition.region2.(extension) = ...
                    mean( area_3regions.region2(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
                area_withRepartition.region3.(extension) = ...
                    mean( area_3regions.region3(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
            end
            
            for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                if ( image_reference_anaphase_3regions(iPeriod_anaphase)< 0 || image_reference_anaphase_3regions(iPeriod_anaphase)> 0  ) && image_reference_anaphase_3regions(iPeriod_anaphase+1)> 0
                    if image_reference_anaphase_3regions(iPeriod_anaphase) < 0
                        image_reference_anaphase_3regions(iPeriod_anaphase) = 1;
                    end
                    area_withRepartition.region1.(extension) = ...
                        mean( area_3regions.region1(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                    area_withRepartition.region2.(extension)= ...
                        mean( area_3regions.region2(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                    area_withRepartition.region3.(extension) = ...
                        mean( area_3regions.region3(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                end
            end
        end
        
        % get area in 10 regions - used later for density map
        if compute_3regions == 1
            for iRegion = 1 : 10
                name = ['nbR' num2str(iRegion)];
                given_area = regionArea.entireEmbryo.(name);
                given_xStart = regionXlimit.entireEmbryo.(name);
                given_deltaX = regionXlength.entireEmbryo.(name);
                for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                    extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                    area_withRepartition.(name).(extension)= ...
                        mean( given_area(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
                    xStart_withRepartition.(name).(extension)= ...
                        mean( given_xStart(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
                    deltaX_withRepartition.(name).(extension)= ...
                        mean( given_deltaX(image_reference_metaphase_3regions(iPeriod_metaphase+1): image_reference_metaphase_3regions(iPeriod_metaphase)) ).* ( (param.resol/1000)^2 );
                end
                
                for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                    extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                    if ( image_reference_anaphase_3regions(iPeriod_anaphase)< 0 || image_reference_anaphase_3regions(iPeriod_anaphase)> 0  ) && image_reference_anaphase_3regions(iPeriod_anaphase+1)> 0
                        if image_reference_anaphase_3regions(iPeriod_anaphase) < 0
                            image_reference_anaphase_3regions(iPeriod_anaphase) = 1;
                        end
                        area_withRepartition.(name).(extension) = ...
                            mean( given_area(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                        xStart_withRepartition.(name).(extension) = ...
                            mean( given_xStart(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                        deltaX_withRepartition.(name).(extension) = ...
                            mean( given_deltaX(image_reference_anaphase_3regions(iPeriod_anaphase)-param.sp2+1: image_reference_anaphase_3regions(iPeriod_anaphase+1)-param.sp2+1) ).* ( (param.resol/1000)^2 );
                        
                    end
                end
            end
        end
        filename = strcat('area_ID' , param.sp1, '.mat');
        save(fullfile(main_path,filename),  '-struct', 'area_withRepartition');
        filename = strcat('xStart_ID' , param.sp1, '.mat');
        save(fullfile(main_path,filename),  '-struct', 'xStart_withRepartition');
        filename = strcat('deltaX_ID' , param.sp1, '.mat');
        save(fullfile(main_path,filename),  '-struct', 'deltaX_withRepartition');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get duration histograms
        
        if compute_wholeEmbryo == 1
            for iPeriod_metaphase = 1 : length(image_reference_metaphase)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase)];
                given_area = area_withRepartition.entireEmbryo.(extension);
                % entire embryo
                dataTracks_input = entireEmbryo.(extension).lengthTracks( entireEmbryo.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_metaphase.(extension).entireEmbryo = [ [embryo_index_metaphase.(extension).entireEmbryo] nbEmbryo ];
                    nbEmbryo_period_metaphase_above250(1,iPeriod_metaphase) = nbEmbryo_period_metaphase_above250(1,iPeriod_metaphase) + 1;
                end
                if nbEmbryo == 1 && iPeriod_metaphase == 1
                    [ tracks_duration_histo_metaphase ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo,[],window_size_metaphase/param.sp6,given_area);
                else
                    [ tracks_duration_histo_metaphase ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo,tracks_duration_histo_metaphase,window_size_metaphase/param.sp6,given_area);
                end
            end
        end
        
        if compute_3regions == 1
            for iPeriod_metaphase = 1 : length(image_reference_metaphase_3regions)-1
                extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                given_area = area_withRepartition.region1.(extension);
                % region 1
                dataTracks_input = region1.(extension).lengthTracks( region1.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_metaphase.(extension).region1 = [ [embryo_index_metaphase.(extension).region1] nbEmbryo_3regions ];
                    nbEmbryo_period_metaphase_above250_3regions.region1(1,iPeriod_metaphase) = nbEmbryo_period_metaphase_above250_3regions.region1(1,iPeriod_metaphase) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_metaphase == 1
                    [ tracks_duration_histo_metaphase_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,[],window_size_metaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_metaphase_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,tracks_duration_histo_metaphase_region1,window_size_metaphase_3regions/param.sp6,given_area);
                end
                % region 2
                given_area = area_withRepartition.region2.(extension);
                dataTracks_input = region2.(extension).lengthTracks( region2.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_metaphase.(extension).region2 = [ [embryo_index_metaphase.(extension).region2] nbEmbryo_3regions ];
                    nbEmbryo_period_metaphase_above250_3regions.region2(1,iPeriod_metaphase) = nbEmbryo_period_metaphase_above250_3regions.region2(1,iPeriod_metaphase) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_metaphase == 1
                    [ tracks_duration_histo_metaphase_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,[],window_size_metaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_metaphase_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,tracks_duration_histo_metaphase_region2,window_size_metaphase_3regions/param.sp6,given_area);
                end
                % region 3
                given_area = area_withRepartition.region3.(extension);
                dataTracks_input = region3.(extension).lengthTracks( region3.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_metaphase.(extension).region3 = [ [embryo_index_metaphase.(extension).region3] nbEmbryo_3regions ];
                    nbEmbryo_period_metaphase_above250_3regions.region3(1,iPeriod_metaphase) = nbEmbryo_period_metaphase_above250_3regions.region3(1,iPeriod_metaphase) + 1;
                end
                
                if nbEmbryo_3regions == 1 && iPeriod_metaphase == 1
                    [ tracks_duration_histo_metaphase_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,[],window_size_metaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_metaphase_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_metaphase,nbEmbryo_3regions,tracks_duration_histo_metaphase_region3,window_size_metaphase_3regions/param.sp6,given_area);
                end
                
            end
        end
        
        if compute_wholeEmbryo == 1
            for iPeriod_anaphase = 1 : length(image_reference_anaphase)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase)];
                given_area = area_withRepartition.entireEmbryo.(extension);
                % entire embryo
                dataTracks_input = entireEmbryo.(extension).lengthTracks( entireEmbryo.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_anaphase.(extension).entireEmbryo = [ [embryo_index_anaphase.(extension).entireEmbryo] nbEmbryo ];
                    nbEmbryo_period_anaphase_above250(1,iPeriod_anaphase) = nbEmbryo_period_anaphase_above250(1,iPeriod_anaphase) + 1;
                end
                if nbEmbryo == 1 && iPeriod_anaphase == 1
                    [ tracks_duration_histo_anaphase ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo,[],window_size_anaphase/param.sp6,given_area);
                else
                    [ tracks_duration_histo_anaphase ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo,tracks_duration_histo_anaphase,window_size_anaphase/param.sp6,given_area);
                end
            end
        end
        
        if compute_3regions == 1
            for iPeriod_anaphase = 1 : length(image_reference_anaphase_3regions)-1
                extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                % region 1
                given_area = area_withRepartition.region1.(extension);
                dataTracks_input = region1.(extension).lengthTracks( region1.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_anaphase.(extension).region1 = [ [embryo_index_anaphase.(extension).region1] nbEmbryo_3regions ];
                    nbEmbryo_period_anaphase_above250_3regions.region1(1,iPeriod_anaphase) = nbEmbryo_period_anaphase_above250_3regions.region1(1,iPeriod_anaphase) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_anaphase == 1
                    [ tracks_duration_histo_anaphase_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,[],window_size_anaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_anaphase_region1 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,tracks_duration_histo_anaphase_region1,window_size_anaphase_3regions/param.sp6,given_area);
                end
                % region 2
                given_area = area_withRepartition.region2.(extension);
                dataTracks_input = region2.(extension).lengthTracks( region2.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_anaphase.(extension).region2 = [ [embryo_index_anaphase.(extension).region2] nbEmbryo_3regions ];
                    nbEmbryo_period_anaphase_above250_3regions.region2(1,iPeriod_anaphase) = nbEmbryo_period_anaphase_above250_3regions.region2(1,iPeriod_anaphase) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_anaphase == 1
                    [ tracks_duration_histo_anaphase_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,[],window_size_anaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_anaphase_region2 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,tracks_duration_histo_anaphase_region2,window_size_anaphase_3regions/param.sp6,given_area);
                end
                % region 3
                given_area = area_withRepartition.region3.(extension);
                dataTracks_input = region3.(extension).lengthTracks( region3.(extension).lengthTracks >= minLength_tracks * param.sp6);
                if length(dataTracks_input) > limit_nb_tracks_for_fitting
                    embryo_index_anaphase.(extension).region3 = [ [embryo_index_anaphase.(extension).region3] nbEmbryo_3regions ];
                    nbEmbryo_period_anaphase_above250_3regions.region3(1,iPeriod_anaphase) = nbEmbryo_period_anaphase_above250_3regions.region3(1,iPeriod_anaphase) + 1;
                end
                if nbEmbryo_3regions == 1 && iPeriod_anaphase == 1
                    [ tracks_duration_histo_anaphase_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,[],window_size_anaphase_3regions/param.sp6,given_area);
                else
                    [ tracks_duration_histo_anaphase_region3 ] = to_get_tracks_duration_histo_simple...
                        ( dataTracks_input,iPeriod_anaphase,nbEmbryo_3regions,tracks_duration_histo_anaphase_region3,window_size_anaphase_3regions/param.sp6,given_area);
                end
            end
        end
        
        if compute_3regions == 1
            dataTracks_repartition.region1 = region1;
            dataTracks_repartition.region2 = region2;
            dataTracks_repartition.region3 = region3;
            clear region1 region2 region3
            clear xStart_3regions xEnd_3regions
        end
        if compute_wholeEmbryo == 1
            dataTracks_repartition.entireEmbryo = entireEmbryo;
            clear entireEmbryo
            clear regionArea_period_anaphase regionArea_period_metaphase area_blastomere
        end
        filename = strcat('dataTracks_withRepartition_ID' , param.sp1, '.mat');
        save(fullfile(main_path,filename), '-struct', 'dataTracks_repartition');
        clear regionArea regionXlength regionXlimit dataTracks_input dataTracks_rotated furrow_position dataTracks_repartition
        
    end
end

% save structure with histo
if compute_wholeEmbryo == 1
    name = strcat('tracks_duration_histo_metaphase_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_metaphase');
    
    name = strcat('tracks_duration_histo_anaphase_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_anaphase');
end

if compute_3regions == 1
    name = strcat('tracks_duration_histo_metaphase_region1_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_metaphase_region1');
    
    name = strcat('tracks_duration_histo_anaphase_region1_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_anaphase_region1');
    
    name = strcat('tracks_duration_histo_metaphase_region2_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name),'tracks_duration_histo_metaphase_region2');
    
    name = strcat('tracks_duration_histo_anaphase_region2_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_anaphase_region2');
    
    name = strcat('tracks_duration_histo_metaphase_region3_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_metaphase_region3');
    
    name = strcat('tracks_duration_histo_anaphase_region3_' , xmlfile_bkp, '.mat');
    save(fullfile(main_path,name), 'tracks_duration_histo_anaphase_region3');
end


%% to perform fitting in each given timewindow and region

if compute_wholeEmbryo == 1
    for iPeriod_metaphase = 1 : max_period_metaphase
        
        extension = ['meta_minus' num2str(iPeriod_metaphase)];
        nbEmbryo_givenCondition = nbEmbryo_period_metaphase_above250(iPeriod_metaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % entire embryo
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_metaphase.(extension).entireEmbryo(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_metaphase{iPeriod_metaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_metaphase{iPeriod_metaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_metaphase{iPeriod_metaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_metaphase{iPeriod_metaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_metaphase{iPeriod_metaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).entireEmbryo.fitting = fitting_results;
            final_fitting_all.(extension).entireEmbryo.model = model_choice;
            clear fitting_results model_choice input
            close all
        end
    end
end

if compute_3regions == 1
    for iPeriod_metaphase = 1 : max_period_metaphase_3regions
        
        extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_metaphase_above250_3regions.region1(iPeriod_metaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % region1
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_metaphase.(extension).region1(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_metaphase_region1{iPeriod_metaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_metaphase_region1{iPeriod_metaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_metaphase_region1{iPeriod_metaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_metaphase_region1{iPeriod_metaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_metaphase_region1{iPeriod_metaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region1.fitting = fitting_results;
            final_fitting_all.(extension).region1.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
        
        % extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_metaphase_above250_3regions.region2(iPeriod_metaphase);
        
        if nbEmbryo_givenCondition > 1
            % region2
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_metaphase.(extension).region2(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_metaphase_region2{iPeriod_metaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_metaphase_region2{iPeriod_metaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_metaphase_region2{iPeriod_metaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_metaphase_region2{iPeriod_metaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_metaphase_region2{iPeriod_metaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region2.fitting = fitting_results;
            final_fitting_all.(extension).region2.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
        
        % extension = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_metaphase_above250_3regions.region3(iPeriod_metaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % region3
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_metaphase.(extension).region3(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_metaphase_region3{iPeriod_metaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_metaphase_region3{iPeriod_metaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_metaphase_region3{iPeriod_metaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_metaphase_region3{iPeriod_metaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_metaphase_region3{iPeriod_metaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region3.fitting = fitting_results;
            final_fitting_all.(extension).region3.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end

if compute_wholeEmbryo == 1
    for iPeriod_anaphase = 1 : max_period_anaphase
        
        extension = ['ana_plus' num2str(iPeriod_anaphase)];
        nbEmbryo_givenCondition = nbEmbryo_period_anaphase_above250(iPeriod_anaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % entire embryo
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_anaphase.(extension).entireEmbryo(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_anaphase{iPeriod_anaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_anaphase{iPeriod_anaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_anaphase{iPeriod_anaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_anaphase{iPeriod_anaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_anaphase{iPeriod_anaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).entireEmbryo.fitting = fitting_results;
            final_fitting_all.(extension).entireEmbryo.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end

if compute_3regions == 1
    for iPeriod_anaphase = 1 : max_period_anaphase_3regions
        
        extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_anaphase_above250_3regions.region1(iPeriod_anaphase);
        
        if nbEmbryo_givenCondition > 1
            % region1
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_anaphase.(extension).region1(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_anaphase_region1{iPeriod_anaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_anaphase_region1{iPeriod_anaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_anaphase_region1{iPeriod_anaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_anaphase_region1{iPeriod_anaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_anaphase_region1{iPeriod_anaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region1.fitting = fitting_results;
            final_fitting_all.(extension).region1.model = model_choice;
            clear fitting_results model_choice input
            close all
        end
        
        % extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_anaphase_above250_3regions.region2(iPeriod_anaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % region2
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_anaphase.(extension).region2(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_anaphase_region2{iPeriod_anaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_anaphase_region2{iPeriod_anaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_anaphase_region2{iPeriod_anaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_anaphase_region2{iPeriod_anaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_anaphase_region2{iPeriod_anaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region2.fitting = fitting_results;
            final_fitting_all.(extension).region2.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
        
        %extension = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
        nbEmbryo_givenCondition = nbEmbryo_period_anaphase_above250_3regions.region3(iPeriod_anaphase);
        
        if nbEmbryo_givenCondition > 1
            
            % region3
            for iiEmbryo = 1:nbEmbryo_givenCondition
                iEmbryo = embryo_index_anaphase.(extension).region3(iiEmbryo);
                name_embryo = ['embryo' num2str(iiEmbryo)];
                input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo_anaphase_region3{iPeriod_anaphase,iEmbryo}.binranges, ...
                    tracks_duration_histo_anaphase_region3{iPeriod_anaphase,iEmbryo}.bincounts ) ) );
                input.(name_embryo).duration_phase = tracks_duration_histo_anaphase_region3{iPeriod_anaphase,iEmbryo}.timeDuration_phase;
                input.(name_embryo).area = tracks_duration_histo_anaphase_region3{iPeriod_anaphase,iEmbryo}.area; % area in squared um
                input.(name_embryo).raw_data = tracks_duration_histo_anaphase_region3{iPeriod_anaphase,iEmbryo}.lengthTracks;
            end
            [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
                ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
            [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );
            final_fitting_all.(extension).region3.fitting = fitting_results;
            final_fitting_all.(extension).region3.model = model_choice;
            clear fitting_results model_choice input
            close all
            
        end
    end
end

% save structure with results incorporated

name = strcat('final_results-BayesianInference_sum_mle' , xmlfile_bkp, '.mat');
save(fullfile(main_path,name), '-struct','final_fitting_all');


%% to perform the assignment

%-------------------------------------------------------------------------------------------
% assignment for each population & get corresponding dataTarcks
nb_embryo = 0;
ID_embryo = {};


for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp(param.sp1);
            
            nb_embryo = nb_embryo +1;
            name0 = ['embryo', num2str(nb_embryo)];
            ID_embryo{nb_embryo} = param.sp1;
            
            name_dataTracks = strcat('dataTracks_withRepartition_ID' , param.sp1, '.mat');
            dataTracks_repartition = load(fullfile(main_path,name_dataTracks));
            
            if compute_3regions == 1 && compute_wholeEmbryo == 1
                conditions1 =  {'entireEmbryo' 'region1' 'region2' 'region3' };
            end
            if compute_3regions == 0 && compute_wholeEmbryo == 1
                conditions1 =  {'entireEmbryo'};
            end
            if compute_3regions == 1 && compute_wholeEmbryo == 0
                conditions1 =  {'region1' 'region2' 'region3' };
                max_period_metaphase = max_period_metaphase_3regions;
                max_period_anaphase = max_period_anaphase_3regions;
            end
            n_condi1 = length(conditions1);
            
            n_condi2 = 0;
            conditions2 = {};
            for iPeriod_metaphase = 1 : max_period_metaphase
                n_condi2 = n_condi2 + 1;
                if compute_3regions == 1 % maybe add the both sets together
                    conditions2{n_condi2} = ['meta_minus' num2str(iPeriod_metaphase) '_3regions'];
                else
                    conditions2{n_condi2} = ['meta_minus' num2str(iPeriod_metaphase)];
                end
            end
            
            for iPeriod_anaphase = 1 : max_period_anaphase
                n_condi2 = n_condi2 + 1;
                if compute_3regions == 1
                    conditions2{n_condi2} = ['ana_plus' num2str(iPeriod_anaphase) '_3regions'];
                else
                    conditions2{n_condi2} = ['ana_plus' num2str(iPeriod_anaphase)];
                end
            end
                       
            for jCondition=1:n_condi2
                name2 = conditions2{jCondition};
                
                if ~isempty(strfind(name2,'_3regions'))
                    conditions1_ =  {'region1' 'region2' 'region3' };
                else
                    conditions1_ =  {'entireEmbryo'};
                end
                n_condi1_ = numel(conditions1_);
                
                for iCondition=1:n_condi1_
                    name1 = conditions1_{iCondition};
                    
                    embryo_present = 0;
                    if  findstr(name2,'meta')
                        embryo_present = ismember(nb_embryo,embryo_index_metaphase.(name2).(name1));
                    end
                    if  findstr(name2,'ana')
                        embryo_present = ismember(nb_embryo,embryo_index_anaphase.(name2).(name1));
                    end
                    
                    if isfield(final_fitting_all,name2) && isfield(final_fitting_all.(name2),name1) && embryo_present
                        
                        fitting_results = final_fitting_all.(name2).(name1).fitting;
                        best_model = 'DoubleExpo';  %final_fitting_all.(name1).(name2).model.best_model_BIC;
                        
                        proba_threshold = NaN;
                        if  findstr(name2,'meta')
                            given_number = find(embryo_index_metaphase.(name2).(name1)==nb_embryo);
                        end
                        if  findstr(name2,'ana')
                            given_number = find(embryo_index_anaphase.(name2).(name1)==nb_embryo);
                        end
                        
                        [ familyAssignement,dataTracks_repartition ] = to_assign_MTs_to_subPopulations_shortOrLong...
                            ( dataTracks_repartition,best_model,fitting_results,name1,name2,[],proba_threshold,[],nb_embryo,main_path,param.sp6,given_number );
                        
                        all_familyAssignment.(name0).(name1).(name2) = familyAssignement;
                        
                    end
                end
            end
            name = strcat('dataTracks_withAssignment_ID' , param.sp1, '.mat');
            save(fullfile(main_path, name), '-struct','dataTracks_repartition');
            clear dataTracks_repartition
    end
end

name = strcat('familyAssignement_' , xmlfile_bkp, '.mat');
save(fullfile(main_path,name),'-struct','all_familyAssignment');


%% to get MT contact count evolution along mitosis and in the 3 regions of each subpopulations


nb_embryo = 0; % count of all embryo for which param.sp9 = param. status in xmlfile

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp(param.sp1);
            
            nb_embryo = nb_embryo +1;
            
            disp(param.sp1);
            
            % download datatracks after assignment
            name_dataTracks = strcat('dataTracks_withAssignment_ID' , ID_embryo{nb_embryo}, '.mat');
            dataTracks_afterAssignment = load(fullfile(main_path,name_dataTracks));
            % download regionArea
            name_dataTracks = strcat('area_ID' , ID_embryo{nb_embryo}, '.mat');
            load(fullfile(main_path,name_dataTracks));
            % download xLimit
            name_dataTracks = strcat('xStart_ID' , ID_embryo{nb_embryo}, '.mat');
            load(fullfile(main_path,name_dataTracks));
            % download xStart
            name_dataTracks = strcat('deltaX_ID' , ID_embryo{nb_embryo}, '.mat');
            load(fullfile(main_path,name_dataTracks));
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get counters in the three regions (not spatial count
            % within regions)
            
            local_conditions = {'short_duration' 'long_duration' 'notAssigned'}; % notice that notAssign are all the tracks (short+long)
            
            for jCondition=1:n_condi1  % for entireEmbryo/region1/region2/region3
                name1 = conditions1{jCondition};
                
                %disp(name1);
                
                for iSet = 1 : numel(local_conditions)
                    
                    given_set = local_conditions{iSet};
                    dataTracks_afterAssignment.numTimePoints = param.sp3-param.sp2 +1;
                    
                    %-----------------------------------------
                    % allocate variables
                    counterEachImageBrut = zeros(dataTracks_afterAssignment.numTimePoints,1);
                    counterEachImagePercent = zeros(dataTracks_afterAssignment.numTimePoints,1);
                    counterEachImageAreaNormalized = zeros(dataTracks_afterAssignment.numTimePoints,1);                  
                    
                    for iCondition=1:n_condi2  % for meta_minus../ana_plus...
                        name2 = conditions2{iCondition};
                        
                        embryo_present = 0;
                        if  findstr(name2,'meta')
                            embryo_present = ismember(nb_embryo,embryo_index_metaphase.(name2).(name1));
                        end
                        if  findstr(name2,'ana')
                            embryo_present = ismember(nb_embryo,embryo_index_anaphase.(name2).(name1));
                        end
                        
                        if isfield(dataTracks_afterAssignment,name1) && isfield(dataTracks_afterAssignment.(name1),name2) && ...
                                isfield(dataTracks_afterAssignment.(name1).(name2),given_set) && embryo_present
                            
                            for yourNumberTrack = 1 : dataTracks_afterAssignment.(name1).(name2).(given_set).numTracks
                                if dataTracks_afterAssignment.(name1).(name2).(given_set).lengthTracks(yourNumberTrack) >= general_param.cortex_analysis.minLength * param.sp6
                                    for iTime =  1 : dataTracks_afterAssignment.numTimePoints
                                        if  ( dataTracks_afterAssignment.(name1).(name2).(given_set).indexXStart(yourNumberTrack) <= param.sp2+iTime ) && ...
                                                ( param.sp2+iTime < dataTracks_afterAssignment.(name1).(name2).(given_set).indexXEnd(yourNumberTrack) )
                                            counterEachImageBrut(iTime,1) = counterEachImageBrut(iTime,1) + 1;
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                    
                    % normalized with embryo area
                    for iTime = 1 : dataTracks_afterAssignment.numTimePoints
                        counterEachImageAreaNormalized(iTime,1) = counterEachImageBrut(iTime,1)/area_withRepartition.raw.(name1)(iTime); % aire en pixel**2
                    end
                    
                    
                    %------------------------------------
                    % save data in a structure for a given embryo
                    countersEachImage.(given_set).brut.(name1)= counterEachImageBrut;
                    countersEachImage.(given_set).area.(name1) = counterEachImageAreaNormalized;
                    clear counterEachImageBrut counterEachImageAreaNormalized
                    
                end
                
                ext = ['number' num2str(nb_embryo)];
                legend_array{nb_embryo} = ['ID ' ID_embryo{nb_embryo}];
                embryo.(ext).(name1).countersArea.notAssigned = countersEachImage.notAssigned.area.(name1);
                embryo.(ext).(name1).countersBrut.notAssigned = countersEachImage.notAssigned.brut.(name1)(:,1);
                embryo.(ext).(name1).countersSmoothed.notAssigned = smooth(countersEachImage.notAssigned.brut.(name1)(:,1),20); %20 frames = 2 sec
                embryo.(ext).(name1).countersArea.short_duration = countersEachImage.short_duration.area.(name1);
                embryo.(ext).(name1).countersBrut.short_duration = countersEachImage.short_duration.brut.(name1)(:,1);
                embryo.(ext).(name1).countersSmoothed.short_duration = smooth(countersEachImage.short_duration.brut.(name1)(:,1),20); %20 frames = 2 sec
                embryo.(ext).(name1).countersArea.long_duration = countersEachImage.long_duration.area.(name1);
                embryo.(ext).(name1).countersBrut.long_duration = countersEachImage.long_duration.brut.(name1)(:,1);
                embryo.(ext).(name1).countersSmoothed.long_duration = smooth(countersEachImage.long_duration.brut.(name1)(:,1),20); %20 frames = 2 sec
                nameData = strcat('countersEachImageCortex_ID', ID_embryo{nb_embryo}, '_min-', num2str(general_param.cortex_analysis.minLength*10) , '.mat');
                save(fullfile(pathMainDirectory,nameData), '-struct', 'countersEachImage');
                clear countersEachImage
            end
            
            %--------------------
            %number of images and time vector
            numTimePoints = size(embryo.(ext).(name1).countersArea.short_duration,1);
            time = 1/param.sp6*[1:numTimePoints];
            
            
            %----------------
            % organize data of each embryo in a structure
            embryo.(ext).name = short_name;
            embryo.(ext).time = time;
            embryo.(ext).onsetCyto = area_withRepartition.furrow_onset ;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get density map within 10 regions along AP axis
            % (spatial resolution)
            
            iZone = 10;
            name = ['nbR' num2str(iZone)];
            area = area_withRepartition.raw.(name);
            xStart = xStart_withRepartition.raw.(name);
            deltaX = deltaX_withRepartition.raw.(name);
            
            if compute_wholeEmbryo == 1
                
                local_conditions = {'short_duration' 'long_duration' 'notAssigned'};
                for iSet = 1 : numel(local_conditions)
                    given_set = local_conditions{iSet};
                    
                    % allocate local counters
                    % remporal & spatial evolutions
                    
                    counterEachImageBrut = zeros(dataTracks_afterAssignment.numTimePoints,iZone+1);
                    counterEachImagePercent = zeros(dataTracks_afterAssignment.numTimePoints,iZone);
                    counterEachImageAreaNormalized = zeros(dataTracks_afterAssignment.numTimePoints,iZone);
                    
                    for iCondition=1:n_condi2  % for meta_minus../ana_plus...
                        name2 = conditions2{iCondition};
                        
                        embryo_present = 0;
                        if  findstr(name2,'meta')
                            embryo_present = ismember(nb_embryo,embryo_index_metaphase.(name2).entireEmbryo);
                        end
                        if  findstr(name2,'ana')
                            embryo_present = ismember(nb_embryo,embryo_index_anaphase.(name2).entireEmbryo);
                        end
                        
                        if isfield(dataTracks_afterAssignment.entireEmbryo,name2) &&  isfield(dataTracks_afterAssignment.entireEmbryo.(name2),given_set) ...
                                && embryo_present
                            
                            dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).numTracks = ...
                                length(dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).lengthTracks);
                            for yourNumberTrack = 1 : dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).numTracks
                                if dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).lengthTracks(yourNumberTrack) >= general_param.cortex_analysis.minLength
                                    for iTime =  1 : dataTracks_afterAssignment.numTimePoints
                                        if  ( dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).indexXStart(yourNumberTrack) <= param.sp2+iTime ) && ...
                                                ( param.sp2+iTime < dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).indexXEnd(yourNumberTrack) )
                                            counterEachImageBrut(iTime,iZone+1) = counterEachImageBrut(iTime,iZone+1) + 1;
                                            for iiZone = 1 : iZone
                                                if ( xStart(iTime,iiZone) <= dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).tracksX(dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).indexXStart(yourNumberTrack)-param.sp2,yourNumberTrack) ) ...
                                                        && ( dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).tracksX(dataTracks_afterAssignment.entireEmbryo.(name2).(given_set).indexXStart(yourNumberTrack)-param.sp2,yourNumberTrack) < xStart(iTime,iiZone)+deltaX(iTime) )
                                                    counterEachImageBrut(iTime,iiZone) = counterEachImageBrut(iTime,iiZone) + 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    for iTime = 1 : dataTracks_afterAssignment.numTimePoints
                        for iiZone = 1 : iZone
                            counterEachImagePercent(iTime,iiZone) = 100*counterEachImageBrut(iTime,iiZone)/counterEachImageBrut(iTime,end);
                            counterEachImageAreaNormalized(iTime,iiZone) = counterEachImageBrut(iTime,iiZone)/area(iTime,iiZone); % aire en pixel**2
                        end
                    end
                    
                    
                    % save different counters in structure
                    countersEachImage.(given_set).brut.(name)= counterEachImageBrut;  % name = nbR10
                    countersEachImage.(given_set).percent.(name)= counterEachImagePercent;
                    countersEachImage.(given_set).area.(name) = counterEachImageAreaNormalized;
                    
                    ext = ['number' num2str(nb_embryo)];
                    embryo.(ext).counters_nbR10.(given_set) = countersEachImage.(given_set).area.(name);
                    
                end
            end
            
            
            if compute_3regions == 1
                
                local_conditions = {'short_duration' 'long_duration' 'notAssigned'};
                for iSet = 1 : numel(local_conditions)
                    given_set = local_conditions{iSet};
                    
                    % allocate local counters
                    % remporal & spatial evolutions
                    
                    counterEachImageBrut = zeros(dataTracks_afterAssignment.numTimePoints,iZone+1);
                    counterEachImagePercent = zeros(dataTracks_afterAssignment.numTimePoints,iZone);
                    counterEachImageAreaNormalized = zeros(dataTracks_afterAssignment.numTimePoints,iZone);
                    
                    for iCondition=1:n_condi2  % for meta_minus../ana_plus...
                        name2 = conditions2{iCondition};
                        
                        conditions1_ =  {'region1' 'region2' 'region3' };
                        for jCondition=1:n_condi1  % for meta_minus../ana_plus...
                            name1 = conditions1_{jCondition};
                            
                            embryo_present = 0;
                            if  findstr(name2,'meta')
                                embryo_present = ismember(nb_embryo,embryo_index_metaphase.(name2).(name1));
                            end
                            if  findstr(name2,'ana')
                                embryo_present = ismember(nb_embryo,embryo_index_anaphase.(name2).(name1));
                            end
                            
                            if isfield(dataTracks_afterAssignment.(name1),name2) &&  isfield(dataTracks_afterAssignment.(name1).(name2),given_set) ...
                                    && embryo_present
                                
                                dataTracks_afterAssignment.(name1).(name2).(given_set).numTracks = ...
                                    length(dataTracks_afterAssignment.(name1).(name2).(given_set).lengthTracks);
                                for yourNumberTrack = 1 : dataTracks_afterAssignment.(name1).(name2).(given_set).numTracks
                                    if dataTracks_afterAssignment.(name1).(name2).(given_set).lengthTracks(yourNumberTrack) >= general_param.cortex_analysis.minLength
                                        for iTime =  1 : dataTracks_afterAssignment.numTimePoints
                                            if  ( dataTracks_afterAssignment.(name1).(name2).(given_set).indexXStart(yourNumberTrack) <= param.sp2+iTime ) && ...
                                                    ( param.sp2+iTime < dataTracks_afterAssignment.(name1).(name2).(given_set).indexXEnd(yourNumberTrack) )
                                                counterEachImageBrut(iTime,iZone+1) = counterEachImageBrut(iTime,iZone+1) + 1;
                                                for iiZone = 1 : iZone
                                                    if ( xStart(iTime,iiZone) <= dataTracks_afterAssignment.(name1).(name2).(given_set).tracksX(dataTracks_afterAssignment.(name1).(name2).(given_set).indexXStart(yourNumberTrack)-param.sp2,yourNumberTrack) ) ...
                                                            && ( dataTracks_afterAssignment.(name1).(name2).(given_set).tracksX(dataTracks_afterAssignment.(name1).(name2).(given_set).indexXStart(yourNumberTrack)-param.sp2,yourNumberTrack) < xStart(iTime,iiZone)+deltaX(iTime) )
                                                        counterEachImageBrut(iTime,iiZone) = counterEachImageBrut(iTime,iiZone) + 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        
                        for iTime = 1 : dataTracks_afterAssignment.numTimePoints
                            for iiZone = 1 : iZone
                                counterEachImagePercent(iTime,iiZone) = 100*counterEachImageBrut(iTime,iiZone)/counterEachImageBrut(iTime,end);
                                counterEachImageAreaNormalized(iTime,iiZone) = counterEachImageBrut(iTime,iiZone)/area(iTime,iiZone); % aire en pixel**2
                            end
                        end
                        
                        
                        % save different counters in structure
                        countersEachImage.(given_set).brut.(name)= counterEachImageBrut;  % name = nbR10
                        countersEachImage.(given_set).percent.(name)= counterEachImagePercent;
                        countersEachImage.(given_set).area.(name) = counterEachImageAreaNormalized;
                        
                        ext = ['number' num2str(nb_embryo)];
                        embryo.(ext).counters_nbR10_3regions.(given_set) = countersEachImage.(given_set).area.(name);
                        
                    end
                end
                
            end
    end  
end
% save embryo selection in structure
name = ['embryoSelection' param.extra 'nb' num2str(nb_embryo) '.mat'];
save(fullfile(main_path,name), '-struct', 'embryo');


%% AVERAGING OF THE COUNTERS

%% averaging of counters in the three regions (no spatial count within
% region)

steps = 10; % in sec
half_width = 1000; % in time unit

%---------------------------------------------------------------------------
%-----------------------------------------------------------------------
% averaging of area normalized counters over different embryos
nb_col = nb_embryo +1;
col=jet(nb_col);

for jCondition=1:n_condi1
    name1 = conditions1{jCondition}; %name1 = /entireEmbryoregion 1/2/3
    
    
    local_conditions = {'short_duration' 'long_duration' 'notAssigned'};
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        figure
        histo = zeros(2*ceil(half_width/steps)+1,3);
        for embryo_idx=1:nb_embryo
            ext2 = ['number' num2str(embryo_idx)];
            time = transpose(embryo.(ext2).time); %transpose( 1/param.sp6*[1: length(embryo.(ext2).entireEmbryo.countersArea.short_duration)] ); %
            data = embryo.(ext2).(name1).countersArea.(given_set);
            time_reference = embryo.(ext2).onsetCyto; % in time unit
            [histo,maxamp]=add_histo_helper_ln(cat(2,time,data),histo,steps,half_width,time_reference); %to average
            amp(:,embryo_idx) = maxamp.data(:,1);
            plot( (time -time_reference) , data,'LineWidth',0.5,'color',col(embryo_idx,:));
            hold all
        end
        legend(legend_array{:})
        ylabel('Number of MTs per \pix^2 ');
        xlabel(sprintf('Time from furrow ingression onset (s)'));
        title(given_set)
        namePlot3 = ['Superposition_cortexCountersAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '-' name1 '-' given_set '.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        % finalize
        amp(:,nb_embryo+1) = histo(:,1);
        amp(amp(:,nb_embryo+1)==0,:) = [];
        amp(:,nb_embryo+2) = nanmean(amp(:,1:nb_embryo),2);
        amp(:,nb_embryo+3) = 0;
        for j = 1 : size(amp,1)
            for i = 1 : nb_embryo
                if ~isnan(amp(j,i))
                    amp(j,nb_embryo+3) = amp(j,nb_embryo+3) + (amp(j,i) - amp(j,nb_embryo+2)).^2;
                end
            end
        end
        for j = 1 : size(amp,1)
            amp(j,nb_embryo+3) = sqrt(amp(j,nb_embryo+3)/amp(j,nb_embryo+1));
        end
        histo_f = finalise_histo_helper_ln(histo,half_width,steps);
        
        % time / average / sem
        final_values(:,1) = histo_f(:,1);
        final_values(:,2) = 1/((param.resol/1000)^2) * amp(:,nb_embryo+2);
        final_values(:,3) = 1/((param.resol/1000)^2) * amp(:,nb_embryo+3); % 1/((param.resol/1000)^2) = 51.757 to put in um2 unit since 1pixel = 139nm
        ext3 = ['threeRegions'];
        averageRegions.(ext3).(name1).area_normalized.(given_set) = final_values;
        
        clear amp histo_f maxamp histo final_values
        
    end
    
    %save data
    name_ = ['averageCountersWholeEmbryo-area' param.extra 'nb' num2str(nb_embryo) '-' name1];
    save(fullfile(main_path,name_), '-struct', 'averageRegions');
    
    
    %-----------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------
    % averaging of instantaneous counters over different embryos
    
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        figure
        histo = zeros(2*ceil(half_width/steps)+1,3);
        for embryo_idx=1:nb_embryo
            ext2 = ['number' num2str(embryo_idx)];
            time = transpose(embryo.(ext2).time);
            data = embryo.(ext2).(name1).countersBrut.(given_set);
            time_reference = embryo.(ext2).onsetCyto; % in time unit
            [histo,maxamp]=add_histo_helper_ln(cat(2,time,data),histo,steps,half_width,time_reference); %to average
            amp(:,embryo_idx) = maxamp.data(:,1);
            plot( (time -time_reference) , data,'LineWidth',0.5,'color',col(embryo_idx,:));
            hold all
        end
        legend(legend_array{:})
        ylabel('Number of MT contacts in visible cortex');
        xlabel(sprintf('Time from furrow ingression onset (s)'));
        title(given_set)
        namePlot3 = ['Superposition_cortexCountersAverage-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set '.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        % finalize
        amp(:,nb_embryo+1) = histo(:,1);
        amp(amp(:,nb_embryo+1)==0,:) = [];
        amp(:,nb_embryo+2) = nanmean(amp(:,1:nb_embryo),2);
        amp(:,nb_embryo+3) = 0;
        for j = 1 : size(amp,1)
            for i = 1 : nb_embryo
                if ~isnan(amp(j,i))
                    amp(j,nb_embryo+3) = amp(j,nb_embryo+3) + (amp(j,i) - amp(j,nb_embryo+2)).^2;
                end
            end
        end
        for j = 1 : size(amp,1)
            amp(j,nb_embryo+3) = sqrt(amp(j,nb_embryo+3)/amp(j,nb_embryo+1));
        end
        
        % time / average / sem
        histo_f=finalise_histo_helper_ln(histo,half_width,steps); % to get correct time reference
        final_values(:,1) = histo_f(:,1);
        final_values(:,2) = amp(:,nb_embryo+2);
        final_values(:,3) = amp(:,nb_embryo+3);
        ext3 = ['threeRegions'];
        averageRegions.(ext3).(name1).brut.(given_set) = final_values;
        clear amp histo_f maxamp histo final_values
        
    end
    
    %save data
    name_ = ['averageCountersWholeEmbryo' param.extra 'nb' num2str(nb_embryo) '_' name1 '.mat'];
    save(fullfile(main_path,name_), '-struct', 'averageRegions');
    
    
    %-------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------
    % averaging of instantaneous counters over different embryos (smooth one)
    
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        %-----------------------------------------
        % all tracks
        figure
        histo = zeros(2*ceil(half_width/steps)+1,3);
        for embryo_idx=1:nb_embryo
            ext2 = ['number' num2str(embryo_idx)];
            time = transpose(embryo.(ext2).time);
            data = embryo.(ext2).(name1).countersSmoothed.(given_set);
            time_reference = embryo.(ext2).onsetCyto; % in time unit
            [histo,maxamp]=add_histo_helper_ln(cat(2,time,data),histo,steps,half_width,time_reference); %to average
            amp(:,embryo_idx) = maxamp.data(:,1);
            plot( (time -time_reference) , smooth(data,10),'LineWidth',0.5,'color',col(embryo_idx,:));
            hold all
        end
        
        % finalize
        amp(:,nb_embryo+1) = histo(:,1);
        amp(amp(:,nb_embryo+1)==0,:) = [];
        amp(:,nb_embryo+2) = nanmean(amp(:,1:nb_embryo),2);
        amp(:,nb_embryo+3) = 0;
        for j = 1 : size(amp,1)
            for i = 1 : nb_embryo
                if ~isnan(amp(j,i))
                    amp(j,nb_embryo+3) = amp(j,nb_embryo+3) + (amp(j,i) - amp(j,nb_embryo+2)).^2;
                end
            end
        end
        for j = 1 : size(amp,1)
            amp(j,nb_embryo+3) = sqrt(amp(j,nb_embryo+3)/amp(j,nb_embryo+1));
        end
        
        % time / average / sem
        histo_f=finalise_histo_helper_ln(histo,half_width,steps); % to get correct time reference
        final_values(:,1) = histo_f(:,1);
        final_values(:,2) = amp(:,nb_embryo+2);
        final_values(:,3) = amp(:,nb_embryo+3);
        errorbar(final_values(1:end,1),final_values(1:end,2),final_values(1:end,3),'k');
        legend(legend_array{:} , 'mean')
        ylabel('Number of MT contacts in visible cortex');
        xlabel(sprintf('Time from furrow ingression onset (s)'));
        title(given_set)
        namePlot3 = ['Superposition_cortexCountersAverage-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set  '_smooth.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        clear amp histo_f maxamp histo final_values
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plotting of averaged results
    
    % plot area normalized average counters
    
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        ext3 = ['threeRegions'];
        time = averageRegions.(ext3).(name1).area_normalized.(given_set)(:,1);
        delta_furrow_ana = param.delta_furrowDetection_anaphaseOnset; % in sec
        for i = 1 : length(time)
            time2(i) = time(i) + delta_furrow_ana;
        end
        
        % plot
        figure
        plot(time2,averageRegions.(ext3).(name1).area_normalized.(given_set)(1:end,2));
        ylabel('Number of MTs per \mum^2 ');
        xlabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N='  num2str(nb_embryo) ' embryos)' '-' given_set ]);
        namePlot3 = ['Mean_cortexCountersAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set '.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        % plot with error bars
        figure
        errorbar(time2,averageRegions.(ext3).(name1).area_normalized.(given_set)(:,2),averageRegions.(ext3).(name1).area_normalized.(given_set)(:,3));
        ylabel('Number of MTs per \mum^2 ');
        xlabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N='  num2str(nb_embryo) ' embryos)' '-' given_set ]);
        namePlot3 = ['Mean_cortexCountersAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set '_withError.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        
        %----------------------------------------------------------------------------------
        % plot counters
        
        time = averageRegions.(ext3).(name1).brut.(given_set)(:,1);
        delta_furrow_ana = param.delta_furrowDetection_anaphaseOnset; % in sec
        for i = 1 : length(time)
            time2(i) = time(i) + delta_furrow_ana;
        end
        
        % plot
        figure
        plot(time2,averageRegions.(ext3).(name1).brut.(given_set)(1:end,2));
        ylabel('Number of MT contacts in visible cortex');
        xlabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N='  num2str(nb_embryo) ' embryos)''-' given_set ]);
        namePlot3 = ['Mean_cortexCountersAverage-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set '.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
        
        % plot with error bars
        figure
        errorbar(time2,averageRegions.(ext3).(name1).brut.(given_set)(1:end,2),averageRegions.(ext3).(name1).brut.(given_set)(:,3));
        ylabel('Numbers of MTs contacts in whole cortex ');
        xlabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N='  num2str(nb_embryo) ' embryos)' '-' given_set ]);
        namePlot3 = ['Mean_cortexCountersAverage-refAna' param.extra 'nb' num2str(nb_embryo) '_' name1 '-' given_set '_withError.fig'];
        saveas(gcf,[save_stem namePlot3]);
        close(gcf)
       
        
    end
    
end


%% to get density map of each assigned population

steps = 10; % in sec
half_width = 1000; % in time unit

if compute_3regions == 1 && compute_wholeEmbryo == 1
    density_map_label = {'counters_nbR10' 'counters_nbR10_3regions'};
    density_map_label_mean = {'averaged_counters_nbR10' 'averaged_counters_nbR10_3regions'};
    density_map_label_mean_f = {'averaged_counters_nbR10_final' 'averaged_counters_nbR10_3regions_final'};
elseif compute_3regions == 0 && compute_wholeEmbryo == 1
    density_map_label = {'counters_nbR10'};
    density_map_label_mean = {'averaged_counters_nbR10'};
    density_map_label_mean_f = {'averaged_counters_nbR10_final'};
elseif compute_3regions == 1 && compute_wholeEmbryo == 0
    density_map_label = {'counters_nbR10_3regions'};
    density_map_label_mean = {'averaged_counters_nbR10_3regions'};
    density_map_label_mean_f = {'averaged_counters_nbR10_3regions_final'};
end

% averaging of counters over different embryos

for iMap = 1 : numel(density_map_label)
    
    map_name = density_map_label{iMap};
    averaged_map_name = density_map_label_mean{iMap};
    averaged_map_name_f = density_map_label_mean_f{iMap};
    
    local_conditions = {'notAssigned' 'short_duration' 'long_duration'};
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        % average over embryos for each of the ten regions
        nbRegions = 10;
        for region_idx = 1 : nbRegions
            
            histo = zeros(2*ceil(half_width/steps)+1,3);
            for embryo_idx=1:nb_embryo
                
                ext2 = ['number' num2str(embryo_idx)];
                time = transpose(embryo.(ext2).time);
                data = embryo.(ext2).(map_name).(given_set);
                time_reference = embryo.(ext2).onsetCyto; % in time unit
                data2 = data(:,region_idx); %nb of  MTs in given region per image
                
                [histo,averaged_data] = add_histo_helper_ln(cat(2,time,data2),histo,steps,half_width,time_reference); %to average along 10s
                
                if region_idx == 1
                    embryo.(ext2).(averaged_map_name).(given_set)(:,1) = averaged_data.time; % average over 10 s duration
                end
                embryo.(ext2).(averaged_map_name).(given_set)(:,region_idx+1) = averaged_data.data;
                %~isnan(maxamp.data)
                clear averaged_data
                
            end
            
            histo_f=finalise_histo_helper_ln(histo,half_width,steps); % to get correct time reference and average over embryos
            
            if region_idx == 1
                mean3Dcounters_nbR10.(map_name).(given_set)(:,1) = histo_f(:,1); % to get time in first column of the matrix
                std3Dcounters_nbR10.(map_name).(given_set)(:,1) = histo_f(:,1); % to get time in first column of the matrix
            end
            
            % use add_Histo_helper_ln + finalise_histo_helper_ln
            % averages of counters over the different embryos for each region
            % 1/((param.resol/1000)^2) = 51.757 to put in um2 unit since 1pixel = 139nm
            
            mean3Dcounters_nbR10.(map_name).(given_set)(:,region_idx +1)= 1/((param.resol/1000)^2) * histo_f(:,2);
            % standard deviation of the averaged counters over the different embryos for each region
            std3Dcounters_nbR10.(map_name).(given_set)(:,region_idx +1)= 1/((param.resol/1000)^2) * histo_f(:,3);
            
            ext3 = ['region' num2str(region_idx)];
            averageRegions.(map_name).(given_set).(ext3).original.mean3Dcounters = mean3Dcounters_nbR10.(map_name).(given_set); % mean of count in given region index, from 1 to 10
            averageRegions.(map_name).(given_set).(ext3).original.std3Dcounters = std3Dcounters_nbR10.(map_name).(given_set);
            averageRegions.(map_name).(given_set).(ext3).original.raw = histo_f;
            
        end
        
        % to get individual map "cleaned" fron NaN values
        for embryo_idx=1:nb_embryo
            ext2 = ['number' num2str(embryo_idx)];
            embryo.(ext2).(averaged_map_name).(given_set)(:,2:end) = 1/((param.resol/1000)^2) * embryo.(ext2).(averaged_map_name).(given_set)(:,2:end); % to have results in proper unit
            index_selection = ~isnan(embryo.(ext2).(averaged_map_name).(given_set)(:,2));
            embryo.(ext2).(averaged_map_name_f).(given_set) = embryo.(ext2).(averaged_map_name).(given_set)(index_selection,:);
            clear index_selection
        end
        
    end
    
end

%save data of averaged counters over all embryos
name2 = ['averagecounters_nbR10Regions-area' param.extra 'nb' num2str(nb_embryo) '.mat'];
save(fullfile(main_path,name2), '-struct', 'averageRegions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot density with ref time is anaphase onset - original_colormap

for iMap = 1 : numel(density_map_label)
    map_name = density_map_label{iMap};
    
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        space = [1/2 : 1 : (nbRegions-1/2)];
        time = mean3Dcounters_nbR10.(map_name).(given_set)(:,1);
        delta_furrow_ana = param.delta_furrowDetection_anaphaseOnset; % in sec
        for i = 1 : length(time)
            time2(i) = time(i) + delta_furrow_ana;
        end
        
        % savind density map displayed below
        
        save([save_stem given_set '_space_vector'], 'space');
        save([save_stem given_set '_time_vector'], 'time');
        save([save_stem given_set '_time_vector2'], 'time2');
        save([save_stem given_set '_meanCounters_all'], 'mean3Dcounters_nbR10');
        
        %-----------------------------------------------------------
        % plot 3D
        
        figure
        mesh(space,time2,mean3Dcounters_nbR10.(map_name).(given_set)(1:end,2:end));
        set(gca,'xTick',[0 :1: nbRegions]);
        colormap(gca,'jet')
        zlabel('Numbers of MTs per \mum^2 ');
        xlabel(sprintf('Spatial region\nnumbered from anterior to posterior'));
        ylabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N='  num2str(nb_embryo) ' embryos)']);
        colorbar;
        namePlot3 = [map_name '_3DAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '_' given_set '.fig'];
        saveas(gcf,[save_stem namePlot3]);
        
    end
end

%---------------------------------------------------
% ask user the ylim range for this analysis
min_ylim = input_perso(['what is the minimum of ylim (time in sec from anaphase onset)?' ], -200);
max_ylim = input_perso(['what is the maximum of ylim (time in sec from anaphase onset)?' ], 200);

nb_classes = input_perso(['Number of classes used to have density map? '],10);
max_density = input_perso(['Give max density in the scale bar for averaged map (all tracks)? '],0.12);
max_density_short = input_perso(['Give max density in the scale bar for averaged map (short_duration)? '],0.04);
max_density_long = input_perso(['Give max density in the scale bar for averaged map (long_durations)? '],0.08);


% find index for which time equal to minYlim and maxYlim
tmp1 = abs(time2 - min_ylim);
tmp2 = abs(time2 - max_ylim);
[V_min Idx_min] = min(tmp1);
if isempty(Idx_min)
    Idx_min = 1;
end
[V_max Idx_max] = min(tmp2);
if isempty(Idx_max)
    Idx_max = length(time);
end

for iMap = 1 : numel(density_map_label)
    map_name = density_map_label{iMap};
    for iSet = 1 : numel(local_conditions)
        given_set = local_conditions{iSet};
        
        %-------------------------------------------------
        % plot 3D counters in 2D
        
        figure
        contourf(space,time2,mean3Dcounters_nbR10.(map_name).(given_set)(1:end,2:end),nb_classes);
        xlim([0 10])
        colormap(gca,'jet')
        xlabel(sprintf('Spatial region numbered from anterior to posterior'));
        ylabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N=' num2str(nb_embryo) ' embryos)']);
        colorbar;
        caxis([0 max_density]);
        namePlot4 = [map_name '_3Din2DAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '-setMax' num2str(max_density) '-class' num2str(nb_classes) '_' given_set '.tif'];
        saveas(gcf,[save_stem namePlot4]);
        namePlot5 = [map_name '_3Din2DAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '-setMax' num2str(max_density) '-class' num2str(nb_classes) '_' given_set '.fig'];
        saveas(gcf,[save_stem namePlot5]);
        
        
        %----------------------------------------------------------
        % plot 3D counters in 2D - reduced time
        
        figure
        contourf(space,time2(Idx_min:Idx_max),mean3Dcounters_nbR10.(map_name).(given_set)(Idx_min:Idx_max,2:end),nb_classes);
        xlim([0 10])
        colormap(gca,'jet')
        xlabel(sprintf('Spatial region numbered from anterior to posterior'));
        ylabel(sprintf('Time from anaphase onset (s)'));
        title(['Evolution of the MT contacts (N=' num2str(nb_embryo) ' embryos)']);
        colorbar;
        caxis([0 max_density]);
        namePlot4 = [map_name '_3Din2DAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '-setMax' num2str(max_density) '-ylim-class' num2str(nb_classes) '_' given_set '.tif'];
        saveas(gcf,[save_stem namePlot4]);
        namePlot5 = [map_name '_3Din2DAverage-area-refAna' param.extra 'nb' num2str(nb_embryo) '-setMax' num2str(max_density) '-ylim-class' num2str(nb_classes) '_' given_set '.fig'];
        saveas(gcf,[save_stem namePlot5]);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to show individual map

for iMap = 1 : numel(density_map_label)
    
    map_name = density_map_label{iMap};
    averaged_map_name = density_map_label_mean{iMap};
    averaged_map_name_f = density_map_label_mean_f{iMap};
    
    if visualize_individual_embryo == 1
        
        disp('INDIVIDUAL EMBRYO MAP WITH AUTOSCALE');
        
        for embryo_idx = 1 : nb_embryo
            
            ext2 = ['number' num2str(embryo_idx)];
            
            time_vector = embryo.(ext2).(averaged_map_name_f).(given_set)(:,1) + param.delta_furrowDetection_anaphaseOnset;
            
            for iSet = 1 : numel(local_conditions)
                given_set = local_conditions{iSet};
                
                data = embryo.(ext2).(averaged_map_name_f).(given_set)(:,2:end);
                
                % find index for which time equal to minYlim and maxYlim
                tmp1 = abs(time_vector - min_ylim);
                tmp2 = abs(time_vector - max_ylim);
                [V_min Idx_min] = min(tmp1);
                if isempty(Idx_min)
                    Idx_min = 1;
                end
                [V_max Idx_max] = min(tmp2);
                if isempty(Idx_max)
                    Idx_max = length(time_vector);
                end
                
                figure
                contourf(space,time_vector(Idx_min:Idx_max),data(Idx_min:Idx_max,:),nb_classes);
                xlim([0 10])
                colormap(gca,'jet')
                xlabel(sprintf('Spatial region numbered from anterior to posterior'));
                ylabel(sprintf('Time from anaphase onset (s)'));
                title(['Evolution of the MT contacts : embryo index ',num2str(embryo_idx) ]);
                colorbar;
                %caxis([0 max_density]);
                namePlot4 = [map_name '_3Din2DAverage-area-refAna' param.extra 'embryo' num2str(embryo_idx) '-autoscale-ylim-class' num2str(nb_classes) '_' given_set '.tif'];
                saveas(gcf,[save_stem namePlot4]);
                namePlot5 = [map_name '_3Din2DAverage-area-refAna' param.extra 'embryo' num2str(embryo_idx) '-autoscale-ylim-class' num2str(nb_classes) '_' given_set '.fig'];
                saveas(gcf,[save_stem namePlot5]);
                clear data
                
            end
        end
        
        disp('INDIVIDUAL EMBRYO MAP WITH SAME SCALE');
        max_density_ = input_perso(['Give max density in the scale bar common to all INDIVIDUAL maps (all tracks)? '],0.12);
        
        for embryo_idx = 1 : nb_embryo
            
            ext2 = ['number' num2str(embryo_idx)];
            
            for iSet = 1 : numel(local_conditions)
                given_set = local_conditions{iSet};
                
                time_vector = embryo.(ext2).(averaged_map_name_f).(given_set)(:,1) + param.delta_furrowDetection_anaphaseOnset;
                data = embryo.(ext2).(averaged_map_name_f).(given_set)(:,2:end);
                
                % find index for which time equal to minYlim and maxYlim
                tmp1 = abs(time_vector - min_ylim);
                tmp2 = abs(time_vector - max_ylim);
                [V_min Idx_min] = min(tmp1);
                if isempty(Idx_min)
                    Idx_min = 1;
                end
                [V_max Idx_max] = min(tmp2);
                if isempty(Idx_max)
                    Idx_max = length(time_vector);
                end
                
                figure
                contourf(space,time_vector(Idx_min:Idx_max),data(Idx_min:Idx_max,:),nb_classes);
                xlim([0 10])
                colormap(gca,'jet')
                xlabel(sprintf('Spatial region numbered from anterior to posterior'));
                ylabel(sprintf('Time from anaphase onset (s)'));
                title(['Evolution of the MT contacts : embryo index ',num2str(embryo_idx) ]);
                colorbar;
                caxis([0 max_density_]);
                namePlot4 = [map_name '_3Din2DAverage-area-refAna' param.extra 'embryo' num2str(embryo_idx) '-setMax' num2str(max_density) '-ylim-class' num2str(nb_classes) '_' given_set '.tif'];
                saveas(gcf,[save_stem namePlot4]);
                namePlot5 = [map_name '_3Din2DAverage-area-refAna' param.extra 'embryo' num2str(embryo_idx) '-setMax' num2str(max_density) '-ylim-class' num2str(nb_classes) '_' given_set '.fig'];
                saveas(gcf,[save_stem namePlot5]);
                clear data
                
            end
            
        end
        
    end
end

close all

end

