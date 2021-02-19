function [ tracks_duration_histo ] = to_get_tracks_duration_histo_finalAnalysis( dataTracks,nbEmbryo,...
    minLength_tracks, tracks_duration_histo,late_fitting,metaphase_fitting,furrow_position,area_blastomere, ...
    conditions1,conditions2,main_path,prophase_fitting,prometaphase_fitting,last_bin_kept,regionArea, maximal_timePoint)

% to get trcaks duration histo when metaphase and anaphase can be studied
% + add late anaphase possibility

if nargin < 14
    last_bin_kept = 1;
end

if nargin < 15
    regionArea = [];
end

if nargin < 16
    maximal_timePoint = NaN;
end

global param
global general_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get duration of each temporal section

if ~isnan(maximal_timePoint)
    timeDuration_entireRecording = (maximal_timePoint - param.sp2 + 1) / param.sp6; % sec (correct?)
else
    timeDuration_entireRecording = dataTracks.numTimePoints * param.decimate / param.sp6; % sec
end
%timeDuration_entireRecording_bis = (param.sp3 - param.sp2 +1 - general_param.cortex_analysis.kalman_shift) /param.sp6;

if ~isempty(furrow_position) % if no time reference (case cold chok), furrow_position is empty
    
if  furrow_position.image_start_detection *param.decimate/param.sp6 > param.delta_furrowDetection_anaphaseOnset % meaning we have the complete anaphase
    %timeDuration_metaphase_bis = (furrow_position.image_start_detection - param.sp2 +1 - general_param.cortex_analysis.kalman_shift) /param.sp6 - ...
    %param.delta_furrowDetection_anaphaseOnset ;
    if late_fitting == 0 
        timeDuration_anaphase = (param.sp3 - furrow_position.image_start_detection * param.decimate   +1 ) /param.sp6 + param.delta_furrowDetection_anaphaseOnset ;
        if prophase_fitting == 0 
            timeDuration_prometaphase = timeDuration_entireRecording - timeDuration_anaphase;
        else
            timeDuration_metaphase = param.delta_early_metaphase;
            timeDuration_prophase = timeDuration_entireRecording - timeDuration_anaphase - timeDuration_metaphase;
        end
    else
        timeDuration_late = (param.sp3 - furrow_position.image_start_detection* param.decimate )/param.sp6 + param.delta_furrowDetection_anaphaseOnset - param.delta_late_anaphase;
        timeDuration_anaphase = param.delta_late_anaphase;
        if prophase_fitting == 0
            timeDuration_prometaphase = timeDuration_entireRecording - timeDuration_late - timeDuration_anaphase;
        else
            % modified 4/11/20
            if ( timeDuration_entireRecording - timeDuration_late - timeDuration_anaphase ) > param.delta_early_metaphase
                timeDuration_metaphase = param.delta_early_metaphase;
                timeDuration_prophase = timeDuration_entireRecording - timeDuration_anaphase - timeDuration_metaphase - timeDuration_late;
            else
                timeDuration_metaphase = timeDuration_entireRecording - timeDuration_late - timeDuration_anaphase ;
                timeDuration_prophase = 0;
            end
        end
    end
else % meaning that we have partial anaphase
    if late_fitting == 0
        timeDuration_prometaphase = 0;
        timeDuration_prophase = 0;
        timeDuration_metaphase = 0;
        timeDuration_anaphase = (param.sp3 - param.sp2 +1 ) /param.sp6 ;
    else
        timeDuration_anaphase = param.delta_late_anaphase - (param.delta_furrowDetection_anaphaseOnset*param.sp6 - furrow_position.image_start_detection* param.decimate )/param.sp6;
        timeDuration_metaphase = 0;
        timeDuration_prometaphase = 0;
        timeDuration_prophase = 0;
        timeDuration_late = timeDuration_entireRecording - timeDuration_anaphase;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the area of each region studied, and along time
% furrow_position.image_start_detection is in image from start recording


if ~isempty(furrow_position) % if no time reference (case cold chok), furrow_position is empty
    
nb_images = size(area_blastomere.area_blastomere,1);

mean_entireRecording = mean(area_blastomere.area_blastomere);
area_anterior.entireRecording = mean_entireRecording(1,1);
area_posterior.entireRecording = mean_entireRecording(1,2);
area_anterior_reduced.entireRecording = mean_entireRecording(1,3);
area_entireEmbryo.entireRecording = area_anterior.entireRecording + area_posterior.entireRecording;

if ~isnan(furrow_position.image_start_detection)
    if  (furrow_position.image_start_detection* param.decimate - param.sp2)/param.sp6 > param.delta_furrowDetection_anaphaseOnset % cyto onset image above duration between this time and m2a
        
        if prophase_fitting == 0
            last_image_prometaphase = round( (timeDuration_prometaphase*param.sp6 + param.sp2 -1)/ param.decimate ); % when 1 is the first recorded frame
            mean_prometaphase = mean(area_blastomere.area_blastomere(1:(last_image_prometaphase- round(param.sp2/ param.decimate) ),:)); % shift from sp2
            area_anterior.prometaphase = mean_prometaphase(1,1) ;
            area_posterior.prometaphase = mean_prometaphase(1,2) ;
            area_anterior_reduced.prometaphase = mean_prometaphase(1,3) ;
            area_entireEmbryo.prometaphase =  area_anterior.prometaphase + area_posterior.prometaphase;
        else
            last_image_metaphase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6/param.decimate  -1 );
            last_image_prophase = round(furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset *param.sp6 + param.delta_early_metaphase *param.sp6)/ param.decimate  -1 );
            if last_image_prophase < param.sp2/param.decimate
                last_image_prophase = round(param.sp2/param.decimate);
                area_anterior.prophase = [];
                area_posterior.prophase = [];
                area_anterior_reduced.prophase = [];
                area_entireEmbryo.prophase = [];
            else
                mean_prophase = mean(area_blastomere.area_blastomere(1:(last_image_prophase-round(param.sp2/param.decimate)),:)); % shift from sp2
                area_anterior.prophase = mean_prophase(1,1) ;
                area_posterior.prophase = mean_prophase(1,2) ;
                area_anterior_reduced.prophase = mean_prophase(1,3) ;
                area_entireEmbryo.prophase =  area_anterior.prophase + area_posterior.prophase;
            end
            mean_metaphase = mean(area_blastomere.area_blastomere((last_image_prophase-round(param.sp2/param.decimate)+1):(last_image_metaphase-round(param.sp2/param.decimate)),:)); % shift from sp2
            area_anterior.metaphase = mean_metaphase(1,1) ;
            area_posterior.metaphase = mean_metaphase(1,2) ;
            area_anterior_reduced.metaphase = mean_metaphase(1,3) ;
            area_entireEmbryo.metaphase =  area_anterior.metaphase + area_posterior.metaphase;            
            
        end
        if late_fitting == 0
            if prophase_fitting == 0
                mean_anaphase = mean(area_blastomere.area_blastomere(last_image_prometaphase+1-round(param.sp2/param.decimate):end,:));
            else
                mean_anaphase = mean(area_blastomere.area_blastomere(last_image_metaphase+1-round(param.sp2/param.decimate):end,:));
            end
            area_anterior.anaphase = mean_anaphase(1,1) ;
            area_posterior.anaphase = mean_anaphase(1,2) ;
            area_anterior_reduced.anaphase = mean_anaphase(1,3) ;
            area_entireEmbryo.anaphase =  area_anterior.anaphase + area_posterior.anaphase;
        else
            last_image_anaphase = round(furrow_position.image_start_detection  - (param.delta_furrowDetection_anaphaseOnset *param.sp6 - param.delta_late_anaphase *param.sp6)/param.decimate -1 );
            if last_image_anaphase > param.sp3/param.decimate
                last_image_anaphase = round(param.sp3/param.decimate);
            end
            if last_image_anaphase > size(area_blastomere.area_blastomere,1) + round(param.sp2/param.decimate)- 1
                last_image_anaphase = size(area_blastomere.area_blastomere,1) + round(param.sp2/param.decimate) - 1;
            end
            if prophase_fitting == 0
                mean_anaphase = mean(area_blastomere.area_blastomere(last_image_prometaphase+1-round(param.sp2/param.decimate):last_image_anaphase-round(param.sp2/param.decimate),:));
            else
                mean_anaphase = mean(area_blastomere.area_blastomere(last_image_metaphase+1-round(param.sp2/param.decimate):last_image_anaphase-round(param.sp2/param.decimate),:));
            end
            area_anterior.anaphase = mean_anaphase(1,1) ;
            area_posterior.anaphase = mean_anaphase(1,2) ;
            area_anterior_reduced.anaphase = mean_anaphase(1,3) ;
            area_entireEmbryo.anaphase =  area_anterior.anaphase + area_posterior.anaphase;
            if size( area_blastomere.area_blastomere(last_image_anaphase+1-round(param.sp2/param.decimate):end,:) , 1 ) > 1
                mean_late = mean(area_blastomere.area_blastomere(last_image_anaphase+1-round(param.sp2/param.decimate):end,:));
            else % only 1 frame in late period
                mean_late = area_blastomere.area_blastomere(last_image_anaphase+1-round(param.sp2/param.decimate):end,:);
            end
            area_anterior.late = mean_late(1,1) ;
            area_posterior.late = mean_late(1,2) ;
            area_anterior_reduced.late = mean_late(1,3) ;
            area_entireEmbryo.late =  area_anterior.late + area_posterior.late;
        end
        
    else % no metaphase
        first_image_anaphase = round( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset*param.sp6/param.decimate ); %when 1 is the first image studied (param.sp2)
        if first_image_anaphase < 1
            first_image_anaphase = 1;
        end
        if late_fitting == 0
            mean_anaphase = mean(area_blastomere.area_blastomere(first_image_anaphase:end,:));
            area_anterior.anaphase = mean_anaphase(1,1) ;
            area_posterior.anaphase = mean_anaphase(1,2) ;
            area_anterior_reduced.anaphase = mean_anaphase(1,3) ;
            area_entireEmbryo.anaphase =  area_anterior.anaphase + area_posterior.anaphase;
        else
            last_image_anaphase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6/param.decimate + ...
                param.delta_late_anaphase *param.sp6/param.decimate -1 ); % when 1 is the first image recorded
            if last_image_anaphase > param.sp3/param.decimate
                last_image_anaphase = round(param.sp3/param.decimate);
            end
            if (last_image_anaphase - param.sp2/param.decimate) > size(area_blastomere.area_blastomere,1)
                last_image_anaphase = size(area_blastomere.area_blastomere,1);
            end
            mean_anaphase = mean(area_blastomere.area_blastomere(first_image_anaphase:last_image_anaphase-round(param.sp2/param.decimate),:));
            area_anterior.anaphase = mean_anaphase(1,1) ;
            area_posterior.anaphase = mean_anaphase(1,2) ;
            area_anterior_reduced.anaphase = mean_anaphase(1,3) ;
            area_entireEmbryo.anaphase =  area_anterior.anaphase + area_posterior.anaphase;
            mean_late = mean(area_blastomere.area_blastomere(last_image_anaphase+1-round(param.sp2/param.decimate):end,:));
            area_anterior.late = mean_late(1,1) ;
            area_posterior.late = mean_late(1,2) ;
            area_anterior_reduced.late = mean_late(1,3) ;
            area_entireEmbryo.late =  area_anterior.late + area_posterior.late;
        end
        
    end
end

clear area_blastomere

else % ex cold chok, no time reference, use regionArea.nbR1 to have area_entireEmbryo.entireRecording 
    area_entireEmbryo.entireRecording = mean( regionArea.entireEmbryo.nbR1 );
end


%% T0 DISPLAY TRACK OF GIVEN LENGTH (for illustration)

if general_param.cortex_analysis.checking_tracks == 1
    figure,
    for i = 1 : dataTracks.entireRecording.entireEmbryo.numTracks
        if dataTracks.entireRecording.entireEmbryo.lengthTracks(i) == 3
            plot(dataTracks.entireRecording.entireEmbryo.tracksX(:,i),dataTracks.entireRecording.entireEmbryo.tracksY(:,i),'b','LineWidth',3);
            hold all
        end
    end
    xlabel ('x coordinate');
    ylabel ('y coordinate');
    figureName = strcat('Display-tracks_3frames-', short_name, '.fig');
    saveas(gcf,fullfile(main_path,figureName));
    figureName = strcat('Display-tracks_3frames-', short_name, '.tif');
    saveas(gcf,fullfile(main_path,figureName));
    close(gcf)
    
    figure,
    for i = 1 : dataTracks.entireRecording.entireEmbryo.numTracks
        if dataTracks.entireRecording.entireEmbryo.lengthTracks(i) == 10
            plot(dataTracks.entireRecording.entireEmbryo.tracksX(:,i),dataTracks.entireRecording.entireEmbryo.tracksY(:,i),'r','LineWidth',3);
            hold all
        end
    end
    xlabel ('x coordinate');
    ylabel ('y coordinate');
    figureName = strcat('Display-tracks_10frames-', short_name, '.fig');
    saveas(gcf,fullfile(main_path,figureName));
    figureName = strcat('Display-tracks_10frames-', short_name, '.tif');
    saveas(gcf,fullfile(main_path,figureName));
    close(gcf)
    
    figure,
    for i = 1 : dataTracks.entireRecording.entireEmbryo.numTracks
        if dataTracks.entireRecording.entireEmbryo.lengthTracks(i) == 30
            plot(dataTracks.entireRecording.entireEmbryo.tracksX(:,i),dataTracks.entireRecording.entireEmbryo.tracksY(:,i),'g','LineWidth',3);
            hold all
        end
    end
    xlabel ('x coordinate');
    ylabel ('y coordinate');
    figureName = strcat('Display-tracks_30frames-', short_name, '.fig');
    saveas(gcf,fullfile(main_path,figureName));
    figureName = strcat('Display-tracks_30frames-', short_name, '.tif');
    saveas(gcf,fullfile(main_path,figureName));
    close(gcf)
end

%% TO GET HISTOGRAM OF TRACKS DURATION

if late_fitting == 0 
    conditions1 = conditions1(~strcmp(conditions1,'late'));
end
if metaphase_fitting == 0
    conditions1 = conditions1(~strcmp(conditions1,'metaphase'));
end
if prophase_fitting == 0
    conditions1 = conditions1(~strcmp(conditions1,'prophase'));
end
if prometaphase_fitting == 0
    conditions1 = conditions1(~strcmp(conditions1,'prometaphase'));
end

n_condi1 = numel(conditions1);
n_condi2 = numel(conditions2);

for iCondition=1:n_condi1
    
    name1 = conditions1{iCondition};
    
    for jCondition=1:n_condi2
        
        name2 = conditions2{jCondition};
        
        % allocate local variables to 0
        % number of tracks above min length from user choice
        numberTracksAboveMinLength = 0;
        % length of each track whose length above min
        lengthTracksAboveMinLength = [];
        
        if dataTracks.(name1).(name2).numTracks > 0
            % get number, index and respective lengths of tracks above min
            for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
                if ~isnan(maximal_timePoint)
                    if dataTracks.(name1).(name2).lengthTracks(yourNumberTrack) >= minLength_tracks * param.sp6/param.decimate && ...
                            dataTracks.(name1).(name2).indexXEnd(1,yourNumberTrack) <= maximal_timePoint
                        numberTracksAboveMinLength = numberTracksAboveMinLength + 1;
                        lengthTracksAboveMinLength= [[ lengthTracksAboveMinLength] dataTracks.(name1).(name2).lengthTracks(yourNumberTrack)];
                    end
                else
                    if dataTracks.(name1).(name2).lengthTracks(yourNumberTrack) >= minLength_tracks * param.sp6/param.decimate
                        numberTracksAboveMinLength = numberTracksAboveMinLength + 1;
                        lengthTracksAboveMinLength= [[ lengthTracksAboveMinLength] dataTracks.(name1).(name2).lengthTracks(yourNumberTrack)];
                    end
                end
            end
            
            if strcmp('entireRecording',name1)
                timeDuration = timeDuration_entireRecording;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.entireRecording; % area in squared um
                elseif strcmp('anterior',name2)
                    area = area_anterior.entireRecording;
                elseif strcmp('posterior',name2)
                    area = area_posterior.entireRecording;
                elseif strcmp('anterior_reduced',name2)
                    area = area_anterior_reduced.entireRecording;
                end
            elseif strcmp('prometaphase',name1)
                timeDuration = timeDuration_prometaphase;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.prometaphase;
                elseif strcmp('anterior',name2)
                    area = area_anterior.prometaphase;
                elseif strcmp('posterior',name2)
                    area = area_posterior.prometaphase;
                elseif strcmp('anterior_reduced',name2)
                    area = area_anterior_reduced.prometaphase;
                end
             elseif strcmp('prophase',name1)
                timeDuration = timeDuration_prophase;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.prophase;
                elseif strcmp('anterior',name2)
                    area = area_anterior.prophase;
                elseif strcmp('posterior',name2)
                    area = area_posterior.prophase;
                elseif strcmp('anterior_reduced',name2)
                    area = area_anterior_reduced.prophase;
                end
            elseif strcmp('metaphase',name1)
                if prophase_fitting == 0
                    timeDuration = timeDuration_prometaphase;
                    if strcmp('entireEmbryo',name2)
                        area = area_entireEmbryo.prometaphase;
                    elseif strcmp('anterior',name2)
                        area = area_anterior.prometaphase;
                    elseif strcmp('posterior',name2)
                        area = area_posterior.prometaphase;
                    elseif strcmp('anterior_reduced',name2)
                        area = area_anterior_reduced.prometaphase;
                    end
                elseif prophase_fitting == 1
                    timeDuration = timeDuration_metaphase;
                    if strcmp('entireEmbryo',name2)
                        area = area_entireEmbryo.metaphase;
                    elseif strcmp('anterior',name2)
                        area = area_anterior.metaphase;
                    elseif strcmp('posterior',name2)
                        area = area_posterior.metaphase;
                    elseif strcmp('anterior_reduced',name2)
                        area = area_anterior_reduced.metaphase;
                    end
                end
            elseif strcmp('anaphase',name1)
                timeDuration = timeDuration_anaphase;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.anaphase;
                elseif strcmp('anterior',name2)
                    area = area_anterior.anaphase;
                elseif strcmp('posterior',name2)
                    area = area_posterior.anaphase;
                elseif strcmp('anterior_reduced',name2)
                    area = area_anterior_reduced.anaphase;
                end
            elseif strcmp('late',name1)
                timeDuration = timeDuration_late;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.late;
                elseif strcmp('anterior',name2)
                    area = area_anterior.late;
                elseif strcmp('posterior',name2)
                    area = area_posterior.late;
                elseif strcmp('anterior_reduced',name2)
                    area = area_anterior_reduced.late;
                end
            end
            
            
            maxLength= max(lengthTracksAboveMinLength);
            minLength= min(lengthTracksAboveMinLength);
            binranges = [minLength : 1 : maxLength];
            
            if maxLength > 100 % 10 sec
                index_to_remove = [(100-minLength+1):(maxLength-minLength+1)];
                binranges(index_to_remove) = [];
            end
            
            [bincounts] = hist(lengthTracksAboveMinLength,binranges);
            
            nb_zero = 0;
            index_zero = [];
            for i = 1 : length(bincounts)
                if bincounts(i) ==0
                    nb_zero = nb_zero+1;
                    index_zero = [[index_zero] i ];
                end
                if nb_zero == 10
                    max_index = i;
                    break
                end
            end
            
            if nb_zero < 10
                max_index = max(index_zero);
            end
            
            max_max_index = length(bincounts);
            index_to_remove2 = [(max_index):(max_max_index)];
            binranges(index_to_remove2) = [];
            bincounts(index_to_remove2) = [];
            
            nb_one = 0;
            index_one = [];
            for i = 1 : length(bincounts)
                if bincounts(i) ==1
                    nb_one = nb_one+1;
                    index_one = [[index_one] i ];
                end
                if nb_one == 3
                    max_index_one = i;
                    break
                end
            end
            
            % addition to avoid troubles if values of 1 is in the first
            % bins
            
            late_index = index_one(index_one> 50);
            if sum(late_index) > 100
                if nb_one < 3
                    max_index_one = max(index_one);
                end
                
                max_max_index_one = length(bincounts);
                
                bincounts_above_max_one = 0;
                for index = max_index_one+1 : max_max_index_one
                    bincounts_above_max_one = bincounts_above_max_one + bincounts(index);
                end
                if bincounts_above_max_one >= 1
                    bincounts_above_max_one = bincounts_above_max_one / (max_max_index_one - max_index_one - 1);
                end
                
                if ~isnan(bincounts_above_max_one)
                    index_to_remove3 = [(max_index_one):(max_max_index_one)];
                    binranges(index_to_remove3) = [];
                    bincounts(index_to_remove3) = [];
                    %   binranges(max_index_one) = max_index_one + ( max_max_index_one - max_index_one)/2;
                    %   bincounts(max_index_one) = bincounts_above_max_one;
                end
            end
            
            % to remove bin/count for which count below treshold
            if general_param.cortex_analysis.use_threshold_count_fit == 1
                nb_belowTreshold = 0;
                index_belowThreshold = [];
                for i = 1 : length(bincounts)
                    if bincounts(i) <= general_param.cortex_analysis.binCount_threshold
                        nb_belowTreshold = nb_belowTreshold+1;
                        index_belowThreshold = [[index_belowThreshold] i ];
                    end
                end
                if nb_belowTreshold > 0
                    binranges(index_belowThreshold) = [];
                    bincounts(index_belowThreshold) = [];
                end
            end
            
            if general_param.cortex_analysis.use_smoothing_fit
                bincounts_smoothed = smooth(bincounts,general_param.cortex_analysis.binCount_smoothW_size);
                figure,
                plot(binranges,bincounts,'xk');
                set(gca,'yscale','log');
                hold all
                plot(binranges,bincounts_smoothed,'or');
                xlabel ('Duration (image nb.)');
                ylabel ('Count (a.u.)');
                figureName = strcat('SmoothedDistribution-', short_name, '.fig');
                saveas(gcf,fullfile(main_path,figureName));
                figureName = strcat('SmoothedDistribution-', short_name, '.tif');
                saveas(gcf,fullfile(main_path,figureName));
                close(gcf)
                clear bincounts
                bincounts(1,:) = bincounts_smoothed;
                clear bincounts_smoothed
            end
            
            % check for any inf in data
            [ ri] = find(~isfinite(bincounts));
            bincounts(ri)= [];
            binranges(ri)= [];
            
            if last_bin_kept == 1
                tracks_duration_histo{nbEmbryo}.(name1).(name2).bincounts = double(bincounts);
                tracks_duration_histo{nbEmbryo}.(name1).(name2).binranges = double( binranges*param.decimate./param.sp6 ); % in sec
                tracks_duration_histo{nbEmbryo}.(name1).(name2).number_tracks = sum(bincounts);
            elseif last_bin_kept == 0
                tracks_duration_histo{nbEmbryo}.(name1).(name2).bincounts = double(bincounts(1:end-1));
                tracks_duration_histo{nbEmbryo}.(name1).(name2).binranges = double( binranges(1:end-1)*param.decimate./param.sp6 ); % in sec
                tracks_duration_histo{nbEmbryo}.(name1).(name2).number_tracks = sum(bincounts(1:end-1));
            end
            tracks_duration_histo{nbEmbryo}.(name1).(name2).timeDuration_phase = timeDuration;
            tracks_duration_histo{nbEmbryo}.(name1).(name2).lengthTracks = lengthTracksAboveMinLength;
            tracks_duration_histo{nbEmbryo}.(name1).(name2).area = area;
            
            clear maxLength
            clear minimLength
            clear lengthTracksAboveMinLength
            clear numberTracksAboveMinLength
            clear bincounts
            clear binranges
            clear timeDuration
            clear area
            
        end
        
    end
    
end

tracks_duration_histo{nbEmbryo}.name_embryo = short_name;

clear area_entireEmbryo area_anterior area_posterior area_anterior_reduced

end

