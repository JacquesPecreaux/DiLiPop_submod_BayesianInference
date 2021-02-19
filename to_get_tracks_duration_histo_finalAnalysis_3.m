function [ tracks_duration_histo ] = to_get_tracks_duration_histo_finalAnalysis_3( dataTracks,nbEmbryo,...
    minLength_tracks, tracks_duration_histo,late_fitting,anaphase_fitting,metaphase_fitting,furrow_position,area_3regions, ...
    conditions1,conditions2,main_path,prophase_fitting,prometaphase_fitting)

% to get trcaks duration histo when metaphase and anaphase can be studied
% + add late anaphase possibility


global param
global general_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get duration of each temporal section

timeDuration_entireRecording = dataTracks.numTimePoints / param.sp6;
%timeDuration_entireRecording_bis = (param.sp3 - param.sp2 +1 - general_param.cortex_analysis.kalman_shift) /param.sp6;

if  furrow_position.image_start_detection/param.sp6 > param.delta_furrowDetection_anaphaseOnset % meaning we have the complete anaphase
    %timeDuration_metaphase_bis = (furrow_position.image_start_detection - param.sp2 +1 - general_param.cortex_analysis.kalman_shift) /param.sp6 - ...
    %param.delta_furrowDetection_anaphaseOnset ;
    if late_fitting == 0 
        timeDuration_anaphase = (param.sp3 - furrow_position.image_start_detection +1 ) /param.sp6 + param.delta_furrowDetection_anaphaseOnset ;
        if prophase_fitting == 0 
            timeDuration_prometaphase = timeDuration_entireRecording - timeDuration_anaphase;
        else
            timeDuration_metaphase = param.delta_early_metaphase;
            timeDuration_prophase = timeDuration_entireRecording - timeDuration_anaphase - timeDuration_metaphase;
        end
    else
        timeDuration_late = (param.sp3 - furrow_position.image_start_detection)/param.sp6 + param.delta_furrowDetection_anaphaseOnset - param.delta_late_anaphase;
        timeDuration_anaphase = param.delta_late_anaphase;
        if prophase_fitting == 0
            timeDuration_prometaphase = timeDuration_entireRecording - timeDuration_late - timeDuration_anaphase;
        else
            timeDuration_metaphase = param.delta_early_metaphase;
            timeDuration_prophase = timeDuration_entireRecording - timeDuration_anaphase - timeDuration_metaphase;
        end
    end
else % meaning that we have partial anaphase
    if late_fitting == 0
        timeDuration_prometaphase = 0;
        timeDuration_prophase = 0;
        timeDuration_metaphase = 0;
        timeDuration_anaphase = (param.sp3 - param.sp2 +1 ) /param.sp6 ;
    else
        timeDuration_anaphase = param.delta_late_anaphase - (param.delta_furrowDetection_anaphaseOnset*param.sp6 - furrow_position.image_start_detection)/param.sp6;
        timeDuration_metaphase = 0;
        timeDuration_prometaphase = 0;
        timeDuration_prophase = 0;
        timeDuration_late = timeDuration_entireRecording - timeDuration_anaphase;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_images = length(area_3regions.region1);  % in pixels**2

area_region1.entireRecording = mean(area_3regions.region1).* ( (param.resol/1000)^2 );  % area in squared um
area_region2.entireRecording = mean(area_3regions.region2).* ( (param.resol/1000)^2 );  % area in squared um
area_region3.entireRecording = mean(area_3regions.region3).* ( (param.resol/1000)^2 );  % area in squared um
area_entireEmbryo.entireRecording = area_region1.entireRecording + area_region2.entireRecording + area_region3.entireRecording;

if  (furrow_position.image_start_detection-param.sp2)/param.sp6 > param.delta_furrowDetection_anaphaseOnset % cyto onset image above duration between this time and m2a
       
    if prophase_fitting == 0
        last_image_prometaphase = round(timeDuration_prometaphase*param.sp6 + param.sp2 -1); % when 1 is the first recorded frame
        area_region1.prometaphase = mean(area_3regions.region1(1:last_image_prometaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.prometaphase = mean(area_3regions.region2(1:last_image_prometaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.prometaphase = mean(area_3regions.region3(1:last_image_prometaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.prometaphase = area_region1.prometaphase + area_region2.prometaphase + area_region3.prometaphase;               
    else
        last_image_metaphase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 -1 );
        last_image_prophase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 - param.delta_early_metaphase *param.sp6 -1 );
        if last_image_prophase < param.sp2
            last_image_prophase = param.sp2;
        end       
        area_region1.prophase = mean(area_3regions.region1(1:last_image_prophase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.prophase = mean(area_3regions.region2(1:last_image_prophase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.prophase = mean(area_3regions.region3(1:last_image_prophase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.prophase = area_region1.prophase + area_region2.prophase + area_region3.prophase;
        
        area_region1.metaphase = mean(area_3regions.region1((last_image_prophase-param.sp2+1):(last_image_metaphase-param.sp2))).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.metaphase = mean(area_3regions.region2((last_image_prophase-param.sp2+1):(last_image_metaphase-param.sp2))).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.metaphase = mean(area_3regions.region3((last_image_prophase-param.sp2+1):(last_image_metaphase-param.sp2))).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.metaphase = area_region1.metaphase + area_region2.metaphase + area_region3.metaphase;
        
    end
    
    if late_fitting == 0
        if prophase_fitting == 0
            area_region1.anaphase = mean(area_3regions.region1(last_image_prometaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region2.anaphase = mean(area_3regions.region2(last_image_prometaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region3.anaphase = mean(area_3regions.region3(last_image_prometaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;
        else
            area_region1.anaphase = mean(area_3regions.region1(last_image_metaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region2.anaphase = mean(area_3regions.region2(last_image_metaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region3.anaphase = mean(area_3regions.region3(last_image_metaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
            area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;
        end        
    else       
        last_image_anaphase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 + param.delta_late_anaphase *param.sp6 -1 );
        if last_image_anaphase > param.sp3
            last_image_anaphase = param.sp3;
        end
        if last_image_anaphase > length(area_3regions.region1) + param.sp2 - 1
            last_image_anaphase = length(area_3regions.region1) + param.sp2 - 1;
        end        
        if prophase_fitting == 0           
            area_region1.anaphase = mean(area_3regions.region1(last_image_prometaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region2.anaphase = mean(area_3regions.region2(last_image_prometaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um;
            area_region3.anaphase = mean(area_3regions.region3(last_image_prometaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
            area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;
            
        else            
            area_region1.anaphase = mean(area_3regions.region1(last_image_metaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
            area_region2.anaphase = mean(area_3regions.region2(last_image_metaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um;
            area_region3.anaphase = mean(area_3regions.region3(last_image_metaphase+1-param.sp2:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
            area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;            
        end        
        area_region1.late = mean(area_3regions.region1(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.late = mean(area_3regions.region2(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.late = mean(area_3regions.region3(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.late =  area_region1.late + area_region2.late + area_region3.late;        
    end
       
else % no metaphase
    first_image_anaphase = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset*param.sp6; %when 1 is the first image studied (param.sp2)
    if first_image_anaphase < 1
        first_image_anaphase = 1;
    end
    if late_fitting == 0
        area_region1.anaphase = mean(area_3regions.region1(first_image_anaphase:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.anaphase = mean(area_3regions.region2(first_image_anaphase:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.anaphase = mean(area_3regions.region3(first_image_anaphase:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;
    else
        last_image_anaphase = round(furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 + ...
            param.delta_late_anaphase *param.sp6 -1 ); % when 1 is the first image recorded
        if last_image_anaphase > param.sp3
            last_image_anaphase = param.sp3;
        end
        if (last_image_anaphase - param.sp2) > length(area_3regions.region1)
            last_image_anaphase = length(area_3regions.region1);
        end
        area_region1.anaphase = mean(area_3regions.region1(first_image_anaphase:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.anaphase = mean(area_3regions.region2(first_image_anaphase:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.anaphase = mean(area_3regions.region3(first_image_anaphase:last_image_anaphase-param.sp2)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.anaphase =  area_region1.anaphase + area_region2.anaphase + area_region3.anaphase;
        area_region1.late = mean(area_3regions.region1(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region2.late = mean(area_3regions.region2(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_region3.late = mean(area_3regions.region3(last_image_anaphase+1-param.sp2:end)).* ( (param.resol/1000)^2 );  % area in squared um
        area_entireEmbryo.late =  area_region1.late + area_region2.late + area_region3.late;
    end
    
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
if anaphase_fitting == 0
    conditions1 = conditions1(~strcmp(conditions1,'anaphase'));
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
        
        % get number, index and respective lengths of tracks above min
        for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
            if dataTracks.(name1).(name2).lengthTracks(yourNumberTrack) >= minLength_tracks * param.sp6
                numberTracksAboveMinLength = numberTracksAboveMinLength + 1;
                lengthTracksAboveMinLength= [[ lengthTracksAboveMinLength] dataTracks.(name1).(name2).lengthTracks(yourNumberTrack)];
            end
        end
        
        if strcmp('entireRecording',name1)
            timeDuration = timeDuration_entireRecording;
            if strcmp('entireEmbryo',name2)
                area = area_entireEmbryo.entireRecording; % area in squared um
            elseif strcmp('region1',name2)
                area = area_region1.entireRecording;
            elseif strcmp('region2',name2)
                area = area_region2.entireRecording;
            elseif strcmp('region3',name2)
                area = area_region3.entireRecording;
            end
        elseif strcmp('prometaphase',name1)
            timeDuration = timeDuration_prometaphase;
            if strcmp('entireEmbryo',name2)
                area = area_entireEmbryo.prometaphase;
            elseif strcmp('region1',name2)
                area = area_region1.prometaphase;
            elseif strcmp('region2',name2)
                area = area_region2.prometaphase;
            elseif strcmp('region3',name2)
                area = area_region3.prometaphase;
            end
        elseif strcmp('prophase',name1)
            timeDuration = timeDuration_prophase;
            if strcmp('entireEmbryo',name2)
                area = area_entireEmbryo.prophase;
            elseif strcmp('region1',name2)
                area = area_region1.prophase;
            elseif strcmp('region2',name2)
                area = area_region2.prophase;
            elseif strcmp('region3',name2)
                area = area_region3.prophase;
            end
        elseif strcmp('metaphase',name1)
            if prophase_fitting == 0
                timeDuration = timeDuration_prometaphase;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.prometaphase;
                elseif strcmp('region1',name2)
                    area = area_region1.prometaphase;
                elseif strcmp('region2',name2)
                    area = area_region2.prometaphase;
                elseif strcmp('region3',name2)
                    area = area_region3.prometaphase;
                end
            elseif prophase_fitting == 1
                timeDuration = timeDuration_metaphase;
                if strcmp('entireEmbryo',name2)
                    area = area_entireEmbryo.metaphase;
                elseif strcmp('region1',name2)
                    area = area_region1.metaphase;
                elseif strcmp('region2',name2)
                    area = area_region2.metaphase;
                elseif strcmp('region3',name2)
                    area = area_region3.metaphase;
                end
            end
        elseif strcmp('anaphase',name1)
            timeDuration = timeDuration_anaphase;
            if strcmp('entireEmbryo',name2)
                area = area_entireEmbryo.anaphase;
            elseif strcmp('region1',name2)
                area = area_region1.anaphase;
            elseif strcmp('region2',name2)
                area = area_region2.anaphase;
            elseif strcmp('region3',name2)
                area = area_region3.anaphase;
            end
        elseif strcmp('late',name1)
            timeDuration = timeDuration_late;
            if strcmp('entireEmbryo',name2)
                area = area_entireEmbryo.late;
            elseif strcmp('region1',name2)
                area = area_region1.late;
            elseif strcmp('region2',name2)
                area = area_region2.late;
            elseif strcmp('region3',name2)
                area = area_region3.late;
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
%            binranges(max_index_one) = max_index_one + ( max_max_index_one - max_index_one)/2;
%            bincounts(max_index_one) = bincounts_above_max_one;
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
        
        tracks_duration_histo{nbEmbryo}.(name1).(name2).bincounts = double(bincounts);
        tracks_duration_histo{nbEmbryo}.(name1).(name2).binranges = double(binranges./param.sp6); % in sec
        tracks_duration_histo{nbEmbryo}.(name1).(name2).number_tracks = sum(bincounts);
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

tracks_duration_histo{nbEmbryo}.name_embryo = short_name;

clear area_entireEmbryo area_anterior area_posterior area_anterior_reduced

end

