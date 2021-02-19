function [ dataTracks_new ] = to_get_dataTracks_for_global_residencyTime_analysis...
    ( dataTracks,furrow_position,regionXlength,regionXlimit,conditions1,conditions2,folder_tag,main_path,polarity_metaphase_50,vector_datum,time_reference_choice )

% function to get proper informations necessary for residency time analysis
% using the final measures saved in param

if nargin < 10
    vector_datum = [];
end
if nargin < 11
    time_reference_choice = 0;
end

global param
global general_param

if  isfield(dataTracks,'dataTracks_rotated') == 1
  %  dataTracks_new.entireRecording.entireEmbryo = dataTracks.dataTracks_rotated.entireEmbryo;
    dataTracks_new.numTimePoints = dataTracks.dataTracks_rotated.numTimePoints;
    dataTracks_.entireEmbryo = dataTracks.dataTracks_rotated.entireEmbryo;
    dataTracks_.numTimePoints = dataTracks.dataTracks_rotated.numTimePoints;
    clear dataTracks
    dataTracks = dataTracks_;
    clear dataTracks_
elseif isfield(dataTracks.entireEmbryo,'entireRecording') == 1 % case of INRIA pipeline whitout classification
   % dataTracks_new.entireRecording.entireEmbryo = dataTracks.entireEmbryo.entireRecording;
    dataTracks_new.numTimePoints = dataTracks.entireEmbryo.entireRecording.numTimePoints;  
    dataTracks_.entireEmbryo = dataTracks.entireEmbryo.entireRecording;
    dataTracks_.numTimePoints = dataTracks.entireEmbryo.entireRecording.numTimePoints;
    clear dataTracks
    dataTracks = dataTracks_;
    clear dataTracks_
else
   % dataTracks_new.entireRecording.entireEmbryo = dataTracks.entireEmbryo;
    dataTracks_new.numTimePoints = dataTracks.numTimePoints;
end

if ismember('late',conditions1)
    late_separation = 1;
else
    late_separation = 0;
end

if ismember('prophase',conditions1) % before NEBD
    early_separation = 1;
else
    early_separation = 0;
end

% if dataTracks_new.entireRecording.entireEmbryo.numTracks == 0
%     dataTracks_new.entireRecording.entireEmbryo.numTracks = length(dataTracks_new.entireRecording.entireEmbryo.lengthTracks);
% end

%% case where geenal_param.corttex_analysis.minLength has been changed

index_above = [];
for i = 1 : dataTracks.entireEmbryo.numTracks
    if dataTracks.entireEmbryo.lengthTracks(i) >= general_param.cortex_analysis.minLength * param.sp6 / param.decimate
        index_above = [  [ index_above ] i ];
    end
end
dataTracks_new.entireRecording.entireEmbryo.numTracks = length(index_above);
%dataTracks_new.numTimePoints = dataTracks.numTimePoints;
dataTracks_new.entireRecording.entireEmbryo.lengthTracks = dataTracks.entireEmbryo.lengthTracks(index_above);
dataTracks_new.entireRecording.entireEmbryo.indexXEnd = dataTracks.entireEmbryo.indexXEnd(index_above);
dataTracks_new.entireRecording.entireEmbryo.indexXStart = dataTracks.entireEmbryo.indexXStart(index_above);
dataTracks_new.entireRecording.entireEmbryo.tracksX = dataTracks.entireEmbryo.tracksX(:,index_above);
dataTracks_new.entireRecording.entireEmbryo.tracksY  = dataTracks.entireEmbryo.tracksY(:,index_above);


%% temporal separation of the tracks along cell division

if time_reference_choice == 0
    
    numTracks_after = 0; % anaphase
    numTracks_before = 0; %prometaphase
    numTracks_late = 0; % late anaphase
    numTracks_after_ = 0; % metaphase
    numTracks_before_ = 0; % prophase
    
    for i = 1 : dataTracks.entireEmbryo.numTracks
        if dataTracks.entireEmbryo.lengthTracks(i) >= general_param.cortex_analysis.minLength * param.sp6 / param.decimate
            
        if dataTracks.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 / param.decimate ) )
            if early_separation == 1  % to look at prophase / metaphase
                if dataTracks.entireEmbryo.indexXStart(1,i) < ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset *param.sp6 )/param.decimate ) ...
                        && dataTracks.entireEmbryo.indexXStart(1,i) >= ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        + param.delta_early_metaphase *param.sp6 ) / param.decimate ) % metaphase
                    numTracks_after_ = numTracks_after_ +1;
                    lengthTracks_after_(numTracks_after_) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after_(numTracks_after_) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after_(numTracks_after_) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after_(:,numTracks_after_) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after_(:,numTracks_after_) = dataTracks.entireEmbryo.tracksY(:,i);
                elseif dataTracks.entireEmbryo.indexXStart(1,i) < ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        + param.delta_early_metaphase *param.sp6 )/ param.decimate ) % prophase
                    numTracks_before_ = numTracks_before_ +1;
                    lengthTracks_before_(numTracks_before_) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_before_(numTracks_before_) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_before_(numTracks_before_) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_before_(:,numTracks_before_) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_before_(:,numTracks_before_) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            else %prometaphase
                numTracks_before = numTracks_before +1;
                lengthTracks_before(numTracks_before) = dataTracks.entireEmbryo.lengthTracks(i);
                indexXEnd_before(numTracks_before) = dataTracks.entireEmbryo.indexXEnd(i);
                indexXStart_before(numTracks_before) = dataTracks.entireEmbryo.indexXStart(i);
                tracksX_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksX(:,i);
                tracksY_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksY(:,i);
            end
        else
            if late_separation == 1 % to look at anaphase / late mitosis
                if dataTracks.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 )/ param.decimate ) ...
                        && dataTracks.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        - param.delta_late_anaphase *param.sp6 )/ param.decimate )
                    numTracks_after = numTracks_after +1;
                    lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
                elseif dataTracks.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        - param.delta_late_anaphase *param.sp6 )/ param.decimate ) % late_anaphase
                    numTracks_late = numTracks_late +1;
                    lengthTracks_late(numTracks_late) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_late(numTracks_late) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_late(numTracks_late) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            else
                if dataTracks.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6 )/ param.decimate )
                    numTracks_after = numTracks_after +1;
                    lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            end
        end
    end
    
    end
elseif time_reference_choice == 1
    
    numTracks_after = 0; % anaphase
    numTracks_before = 0; %prometaphase
    numTracks_late = 0; % late anaphase
    numTracks_after_ = 0; % metaphase
    numTracks_before_ = 0; % prophase
    
    for i = 1 : dataTracks.entireEmbryo.numTracks
        if dataTracks.entireEmbryo.lengthTracks(i) >= general_param.cortex_analysis.minLength * param.sp6 / param.decimate
            
        if dataTracks.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase )*param.sp6 )/ param.decimate
            if early_separation == 1  % to look at prophase / metaphase
                if dataTracks.entireEmbryo.indexXStart(1,i) < ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6/ param.decimate ) ...
                        && dataTracks.entireEmbryo.indexXStart(1,i) >= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage )*param.sp6/ param.decimate ) % metaphase
                    numTracks_after_ = numTracks_after_ +1;
                    lengthTracks_after_(numTracks_after_) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after_(numTracks_after_) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after_(numTracks_after_) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after_(:,numTracks_after_) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after_(:,numTracks_after_) = dataTracks.entireEmbryo.tracksY(:,i);
                elseif dataTracks.entireEmbryo.indexXStart(1,i) < ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage )*param.sp6/ param.decimate ) % prophase
                    numTracks_before_ = numTracks_before_ +1;
                    lengthTracks_before_(numTracks_before_) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_before_(numTracks_before_) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_before_(numTracks_before_) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_before_(:,numTracks_before_) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_before_(:,numTracks_before_) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            else %prometaphase
                numTracks_before = numTracks_before +1;
                lengthTracks_before(numTracks_before) = dataTracks.entireEmbryo.lengthTracks(i);
                indexXEnd_before(numTracks_before) = dataTracks.entireEmbryo.indexXEnd(i);
                indexXStart_before(numTracks_before) = dataTracks.entireEmbryo.indexXStart(i);
                tracksX_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksX(:,i);
                tracksY_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksY(:,i);
            end
        else
            if late_separation == 1 % to look at anaphase / late mitosis
                if dataTracks.entireEmbryo.indexXStart(1,i) > ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6/ param.decimate ) ...
                        && dataTracks.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase + ...
                        param.delta_late_anaphase ) *param.sp6/ param.decimate ) % anaphase
                    numTracks_after = numTracks_after +1;
                    lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
                elseif dataTracks.entireEmbryo.indexXStart(1,i) > ( param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase + ...
                        param.delta_late_anaphase ) *param.sp6/ param.decimate  % late_anaphase
                    numTracks_late = numTracks_late +1;
                    lengthTracks_late(numTracks_late) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_late(numTracks_late) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_late(numTracks_late) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            else
                if dataTracks.entireEmbryo.indexXStart(1,i) > ( param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6/ param.decimate
                    
                    numTracks_after = numTracks_after +1;
                    lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
                end
            end
        end
        end
    
    end
end

%save data from tracking in structure file
if ~isnan(time_reference_choice)
    
    if numTracks_after > 0 && ~isempty(find(contains('anaphase',conditions1)))
        dataTracks_anaphase.entireEmbryo.tracksX = tracksX_after;
        dataTracks_anaphase.entireEmbryo.tracksY = tracksY_after;
        dataTracks_anaphase.entireEmbryo.lengthTracks = lengthTracks_after;
        dataTracks_anaphase.entireEmbryo.indexXStart = indexXStart_after;
        dataTracks_anaphase.entireEmbryo.indexXEnd = indexXEnd_after;
        dataTracks_anaphase.entireEmbryo.numTracks = numTracks_after;
    else
        dataTracks_anaphase.entireEmbryo.numTracks = 0;
    end
    
    if numTracks_before > 0 && ~isempty(find(contains('prometaphase',conditions1)))
        dataTracks_prometaphase.entireEmbryo.tracksX = tracksX_before;
        dataTracks_prometaphase.entireEmbryo.tracksY = tracksY_before;
        dataTracks_prometaphase.entireEmbryo.lengthTracks = lengthTracks_before;
        dataTracks_prometaphase.entireEmbryo.indexXStart = indexXStart_before;
        dataTracks_prometaphase.entireEmbryo.indexXEnd = indexXEnd_before;
        dataTracks_prometaphase.entireEmbryo.numTracks = numTracks_before;
    else
        dataTracks_prometaphase.entireEmbryo.numTracks = 0;
    end
    
    if numTracks_after_ > 0 && ~isempty(find(contains('metaphase',conditions1)))
        dataTracks_metaphase.entireEmbryo.tracksX = tracksX_after_;
        dataTracks_metaphase.entireEmbryo.tracksY = tracksY_after_;
        dataTracks_metaphase.entireEmbryo.lengthTracks = lengthTracks_after_;
        dataTracks_metaphase.entireEmbryo.indexXStart = indexXStart_after_;
        dataTracks_metaphase.entireEmbryo.indexXEnd = indexXEnd_after_;
        dataTracks_metaphase.entireEmbryo.numTracks = numTracks_after_;
    else
        dataTracks_metaphase.entireEmbryo.numTracks = 0;
    end
    
    if numTracks_before_ > 0 && ~isempty(find(contains('prophase',conditions1)))
        dataTracks_prophase.entireEmbryo.tracksX = tracksX_before_;
        dataTracks_prophase.entireEmbryo.tracksY = tracksY_before_;
        dataTracks_prophase.entireEmbryo.lengthTracks = lengthTracks_before_;
        dataTracks_prophase.entireEmbryo.indexXStart = indexXStart_before_;
        dataTracks_prophase.entireEmbryo.indexXEnd = indexXEnd_before_;
        dataTracks_prophase.entireEmbryo.numTracks = numTracks_before_;
    else
        dataTracks_prophase.entireEmbryo.numTracks = 0;
    end
    
    if numTracks_late > 0 && ~isempty(find(contains('late',conditions1)))
        dataTracks_late.entireEmbryo.tracksX = tracksX_late;
        dataTracks_late.entireEmbryo.tracksY = tracksY_late;
        dataTracks_late.entireEmbryo.lengthTracks = lengthTracks_late;
        dataTracks_late.entireEmbryo.indexXStart = indexXStart_late;
        dataTracks_late.entireEmbryo.indexXEnd = indexXEnd_late;
        dataTracks_late.entireEmbryo.numTracks = numTracks_late;
    else
        dataTracks_late.entireEmbryo.numTracks = 0;
    end
    
    
    clear tracksX_before
    clear lengthTracks_before
    clear tracksY_before
    clear lengthTracks_before
    clear indexXStart_before
    clear indexXEnd_before
    
    clear tracksX_after
    clear tracksY_after
    clear lengthTracks_after
    clear indexXStart_after
    clear indexXEnd_after
    
    clear tracksX_before_
    clear lengthTracks_before_
    clear tracksY_before_
    clear lengthTracks_before_
    clear indexXStart_before_
    clear indexXEnd_before_
    
    clear tracksX_after_
    clear tracksY_after_
    clear lengthTracks_after_
    clear indexXStart_after_
    clear indexXEnd_after_
    
    clear tracksX_late
    clear tracksY_late
    clear lengthTracks_late
    clear indexXStart_late
    clear indexXEnd_late
    
else
    
    numTracks_after = 0; % anaphase
    numTracks_before = 0; %prometaphase
    numTracks_late = 0; % late anaphase
    numTracks_after_ = 0; % metaphase
    numTracks_before_ = 0; % prophase
    
end

%% save in the proper format for later use

if ~isempty(find(contains('prophase',conditions1)))
    if dataTracks_prophase.entireEmbryo.numTracks > 0
        dataTracks_new.prophase.entireEmbryo = dataTracks_prophase.entireEmbryo;
    else
        dataTracks_new.prophase.entireEmbryo.numTracks = dataTracks_prophase.entireEmbryo.numTracks;
    end
end

if ~isempty(find(contains('metaphase',conditions1))) 
    if dataTracks_metaphase.entireEmbryo.numTracks > 0
        dataTracks_new.metaphase.entireEmbryo = dataTracks_metaphase.entireEmbryo;
    else
        dataTracks_new.metaphase.entireEmbryo.numTracks = dataTracks_metaphase.entireEmbryo.numTracks;
    end
end

if ~isempty(find(contains('prometaphase',conditions1))) 
    if dataTracks_prometaphase.entireEmbryo.numTracks > 0
        dataTracks_new.prometaphase.entireEmbryo = dataTracks_prometaphase.entireEmbryo;
    else
        dataTracks_new.prometaphase.entireEmbryo.numTracks = dataTracks_prometaphase.entireEmbryo.numTracks;
    end
end

if ~isempty(find(contains('anaphase',conditions1)))
    if dataTracks_anaphase.entireEmbryo.numTracks > 0
        dataTracks_new.anaphase.entireEmbryo = dataTracks_anaphase.entireEmbryo;
    else
        dataTracks_new.anaphase.entireEmbryo.numTracks = dataTracks_anaphase.entireEmbryo.numTracks;
    end
end

if ~isempty(find(contains('late',conditions1)))
    if dataTracks_late.entireEmbryo.numTracks > 0
        dataTracks_new.late.entireEmbryo = dataTracks_late.entireEmbryo;
    else
        dataTracks_new.late.entireEmbryo.numTracks = dataTracks_late.entireEmbryo.numTracks;
    end
end


%% display the tracks in the different time of mitosis

if general_param.cortex_analysis.checking_tracks == 1
    
        % in prophase
    if ~isempty(find(contains('prophase',conditions1)))
        if dataTracks_new.prophase.entireEmbryo.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prophase.entireEmbryo.numTracks
                plot(dataTracks_new.prophase.entireEmbryo.tracksX(:,i),dataTracks_new.prophase.entireEmbryo.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prophase-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prophase-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    % in metaphase
    if ~isempty(find(contains('metaphase',conditions1))) 
        if dataTracks_new.metaphase.entireEmbryo.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.metaphase.entireEmbryo.numTracks
                plot(dataTracks_new.metaphase.entireEmbryo.tracksX(:,i),dataTracks_new.metaphase.entireEmbryo.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_metaphase-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_metaphase-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
        % in prometaphase
    if ~isempty(find(contains('prometaphase',conditions1))) 
        if dataTracks_new.prometaphase.entireEmbryo.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prometaphase.entireEmbryo.numTracks
                plot(dataTracks_new.prometaphase.entireEmbryo.tracksX(:,i),dataTracks_new.prometaphase.entireEmbryo.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prometaphase-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prometaphase-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    % in late anaphase
    if ~isempty(find(contains('late',conditions1)))
        if dataTracks_new.late.entireEmbryo.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.late.entireEmbryo.numTracks
                plot(dataTracks_new.late.entireEmbryo.tracksX(:,i),dataTracks_new.late.entireEmbryo.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_late-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_late-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    % first 2 min anaphase
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.entireEmbryo.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.entireRecording.entireEmbryo.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6/ param.decimate  ) > dataTracks_new.entireRecording.entireEmbryo.indexXStart(i) && ...
                        dataTracks_new.entireRecording.entireEmbryo.indexXStart(1,i) >= ( furrow_position.image_start_detection - ( param.delta_furrowDetection_anaphaseOnset *param.sp6/ param.decimate ) )
                    plot(dataTracks_new.entireRecording.entireEmbryo.tracksX(:,i),dataTracks_new.entireRecording.entireEmbryo.tracksY(:,i))
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_first-2min-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_first-2min-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    % next 2 min in anaphse
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.entireEmbryo.numTracks > 0
            count = 0;
            figure
            for i = 1 : dataTracks_new.entireRecording.entireEmbryo.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 240) *param.sp6/ param.decimate  ) > dataTracks_new.entireRecording.entireEmbryo.indexXStart(i) && ...
                        dataTracks_new.entireRecording.entireEmbryo.indexXStart(1,i) >= ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6/ param.decimate )
                    plot(dataTracks_new.entireRecording.entireEmbryo.tracksX(:,i),dataTracks_new.entireRecording.entireEmbryo.tracksY(:,i))
                    count = count +1;
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_last-2min-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_last-2min-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
end


%% spatial separation of the tracks

 if ~isempty(find(contains('anterior',conditions2))) || ~isempty(find(contains('posterior',conditions2)))
    
    
    % focus on prometaphase / anaphase
    if numTracks_before > 0 && numTracks_after > 0 && ~isempty(find(contains('prometaphase',conditions1))) && ~isempty(find(contains('anaphase',conditions1)))
        
        [dataTracks,dataTracks_anaphase,dataTracks_prometaphase] = separate_anterior_posterior_tracks...
            (dataTracks,dataTracks_anaphase,dataTracks_prometaphase,regionXlength,regionXlimit,furrow_position,2,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.entireRecording.anterior = dataTracks.anterior;
        dataTracks_new.entireRecording.posterior = dataTracks.posterior;
        dataTracks_new.anaphase.anterior = dataTracks_anaphase.anterior;
        dataTracks_new.anaphase.posterior = dataTracks_anaphase.posterior;
        dataTracks_new.prometaphase.anterior = dataTracks_prometaphase.anterior;
        dataTracks_new.prometaphase.posterior = dataTracks_prometaphase.posterior;
        
    elseif numTracks_before > 0 && numTracks_after == 0 && ~isempty(find(contains('prometaphase',conditions1)))
        
        [dataTracks,~,dataTracks_prometaphase] = separate_anterior_posterior_tracks...
            (dataTracks,[],dataTracks_prometaphase,regionXlength,regionXlimit,furrow_position,2,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.entireRecording.anterior = dataTracks.anterior;
        dataTracks_new.entireRecording.posterior = dataTracks.posterior;
        dataTracks_new.prometaphase.anterior = dataTracks_prometaphase.anterior;
        dataTracks_new.prometaphase.posterior = dataTracks_prometaphase.posterior;
        
    elseif numTracks_before == 0 && numTracks_after > 0 && ~isempty(find(contains('anaphase',conditions1)))
        
        [dataTracks,dataTracks_anaphase] = separate_anterior_posterior_tracks...
            (dataTracks,dataTracks_anaphase,[],regionXlength,regionXlimit,furrow_position,2,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.entireRecording.anterior = dataTracks.anterior;
        dataTracks_new.entireRecording.posterior = dataTracks.posterior;
        dataTracks_new.anaphase.anterior = dataTracks_anaphase.anterior;
        dataTracks_new.anaphase.posterior = dataTracks_anaphase.posterior;
        
    else
        
        [dataTracks] = separate_anterior_posterior_tracks...
            (dataTracks,[],[],regionXlength,regionXlimit,furrow_position,0,vector_datum);
        
        dataTracks_new.entireRecording.anterior = dataTracks.anterior;
        dataTracks_new.entireRecording.posterior = dataTracks.posterior;
        
    end
    
    % focus on late
    if numTracks_late > 0 && ~isempty(find(contains('late',conditions1)))
        
        [dataTracks,dataTracks_late] = separate_anterior_posterior_tracks...
            (dataTracks,dataTracks_late,[],regionXlength,regionXlimit,furrow_position,3,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.late.anterior = dataTracks_late.anterior;
        dataTracks_new.late.posterior = dataTracks_late.posterior;        
        
    end

    % focus on prophase/metaphase
    if numTracks_before_ > 0 && numTracks_after_ > 0 && ~isempty(find(contains('metaphase',conditions1))) && ~isempty(find(contains('prophase',conditions1)))
        
        [dataTracks,dataTracks_metaphase,dataTracks_prophase] = separate_anterior_posterior_tracks...
            (dataTracks,dataTracks_metaphase,dataTracks_prophase,regionXlength,regionXlimit,furrow_position,4,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.entireRecording.anterior = dataTracks.anterior;
        dataTracks_new.entireRecording.posterior = dataTracks.posterior;
        dataTracks_new.metaphase.anterior = dataTracks_metaphase.anterior;
        dataTracks_new.metaphase.posterior = dataTracks_metaphase.posterior;
        dataTracks_new.prophase.anterior = dataTracks_prophase.anterior;
        dataTracks_new.prophase.posterior = dataTracks_prophase.posterior;      
        
    elseif numTracks_before_ > 0 && numTracks_after_ == 0 && ~isempty(find(contains('prophase',conditions1)))
        
        [dataTracks,~,dataTracks_prophase] = separate_anterior_posterior_tracks...
            (dataTracks,[],dataTracks_prophase,regionXlength,regionXlimit,furrow_position,4,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.prophase.anterior = dataTracks_prophase.anterior;
        dataTracks_new.prophase.posterior = dataTracks_prophase.posterior;  
        
    elseif numTracks_before_ == 0 && numTracks_after_ > 0 && ~isempty(find(contains('metaphase',conditions1)))
        
        [dataTracks,dataTracks_metaphase] = separate_anterior_posterior_tracks...
            (dataTracks,dataTracks_metaphase,[],regionXlength,regionXlimit,furrow_position,4,polarity_metaphase_50,1,vector_datum);
        
        dataTracks_new.metaphase.anterior = dataTracks_metaphase.anterior;
        dataTracks_new.metaphase.posterior = dataTracks_metaphase.posterior; 
           
    end
    
end

