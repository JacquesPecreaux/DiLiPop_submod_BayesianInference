function [ dataTracks_new,area_3regions,xStart_3regions,xEnd_3regions ] = to_get_dataTracks_for_global_residencyTime_analysis_3...
    ( dataTracks,furrow_position,regionXlength,regionXlimit,regionArea,conditions1,main_path,vector_datum,time_reference_choice )

% function to get proper informations necessary for residency time analysis
% using the final measures saved in param

if nargin < 8
    vector_datum = [];
end
if nargin < 9
    time_reference_choice = 0;
end

global param
global general_param

if  isfield(dataTracks,'dataTracks_rotated') == 1
    dataTracks_new.entireRecording.entireEmbryo = dataTracks.dataTracks_rotated.entireEmbryo;
    dataTracks_new.numTimePoints = dataTracks.dataTracks_rotated.numTimePoints;
    dataTracks_.entireEmbryo = dataTracks.dataTracks_rotated.entireEmbryo;
    dataTracks_.numTimePoints = dataTracks.dataTracks_rotated.numTimePoints;
    clear dataTracks
    dataTracks = dataTracks_;
    clear dataTracks_
else
    dataTracks_new.entireRecording.entireEmbryo = dataTracks.entireEmbryo;
    dataTracks_new.numTimePoints = dataTracks.numTimePoints;
end
% to correct when saving of dataTracks done as dataTracks.entireEmbryo.entireRecording    
if  isfield(dataTracks.entireEmbryo,'entireRecording') == 1 % account for data not save properly.
    dataTracks_.entireEmbryo = dataTracks.entireEmbryo.entireRecording;
    dataTracks_.numTimePoints = dataTracks.numTimePoints;
    clear dataTracks
    dataTracks = dataTracks_;
    clear dataTracks_
end

limit1 = general_param.cortex_analysis.analysis_3SpatialRegions_limit1;
limit2 = general_param.cortex_analysis.analysis_3SpatialRegions_limit2;
[ area_3regions,xStart_3regions,xEnd_3regions ] = get_infos_3regions( regionXlimit,regionXlength,regionArea,limit1,limit2 ); % area in squared pixels

if ismember('late',conditions1)
    late_separation = 1;
else
    late_separation = 0;
end

if ismember('prophase',conditions1)
    early_separation = 1;
else
    early_separation = 0;
end

% will consider datatracks.entireEmbryo.tracksX (and not dataTracks.entireEmbryo.entireRecordin.tracksX)
[dataTracks_rotated,dataTracks_rotated_afterFurrow,dataTracks_rotated_beforeFurrow,dataTracks_rotated_anaphase,dataTracks_rotated_metaphase,...
    dataTracks_rotated_late, dataTracks_rotated_prometaphase, dataTracks_rotated_prophase] ...
    = separate_tracks_3regions(dataTracks,xStart_3regions,xEnd_3regions,furrow_position,late_separation,early_separation,vector_datum,time_reference_choice);

dataTracks_new.entireRecording = dataTracks_rotated;
dataTracks_new.anaphase = dataTracks_rotated_anaphase;
dataTracks_new.metaphase = dataTracks_rotated_metaphase;
dataTracks_new.late = dataTracks_rotated_late;
dataTracks_new.prometaphase = dataTracks_rotated_prometaphase;
dataTracks_new.prophase = dataTracks_rotated_prophase;
   

%% display the tracks in the different time of mitosis

if general_param.cortex_analysis.checking_tracks == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in metaphase
    if ~isempty(find(contains('metaphase',conditions1)))
        if dataTracks_new.metaphase.region1.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.metaphase.region1.numTracks
                plot(dataTracks_new.metaphase.region1.tracksX(:,i),dataTracks_new.metaphase.region1.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_metaphase_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_metaphase_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('metaphase',conditions1)))
        if dataTracks_new.metaphase.region2.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.metaphase.region2.numTracks
                plot(dataTracks_new.metaphase.region2.tracksX(:,i),dataTracks_new.metaphase.region2.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_metaphase_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_metaphase_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('metaphase',conditions1)))
        if dataTracks_new.metaphase.region3.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.metaphase.region3.numTracks
                plot(dataTracks_new.metaphase.region3.tracksX(:,i),dataTracks_new.metaphase.region3.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_metaphase_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_metaphase_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in prometaphase
    if ~isempty(find(contains('prometaphase',conditions1)))
        if dataTracks_new.prometaphase.region1.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prometaphase.region1.numTracks
                plot(dataTracks_new.prometaphase.region1.tracksX(:,i),dataTracks_new.prometaphase.region1.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prometaphase_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prometaphase_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('prometaphase',conditions1)))
        if dataTracks_new.prometaphase.region2.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prometaphase.region2.numTracks
                plot(dataTracks_new.prometaphase.region2.tracksX(:,i),dataTracks_new.prometaphase.region2.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prometaphase_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prometaphase_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('prometaphase',conditions1)))
        if dataTracks_new.prometaphase.region3.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prometaphase.region3.numTracks
                plot(dataTracks_new.prometaphase.region3.tracksX(:,i),dataTracks_new.prometaphase.region3.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prometaphase_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prometaphase_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in prophase
    if ~isempty(find(contains('prophase',conditions1)))
        if dataTracks_new.prophase.region1.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prophase.region1.numTracks
                plot(dataTracks_new.prophase.region1.tracksX(:,i),dataTracks_new.prophase.region1.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prophase_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prophase_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('prophase',conditions1)))
        if dataTracks_new.prophase.region2.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prophase.region2.numTracks
                plot(dataTracks_new.prophase.region2.tracksX(:,i),dataTracks_new.prophase.region2.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prophase_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prophase_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('prophase',conditions1)))
        if dataTracks_new.prophase.region3.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.prophase.region3.numTracks
                plot(dataTracks_new.prophase.region3.tracksX(:,i),dataTracks_new.prophase.region3.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_prophase_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_prophase_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in anaphase
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region1.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.anaphase.region1.numTracks
                plot(dataTracks_new.anaphase.region1.tracksX(:,i),dataTracks_new.anaphase.region1.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region2.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.anaphase.region2.numTracks
                plot(dataTracks_new.anaphase.region2.tracksX(:,i),dataTracks_new.anaphase.region2.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region3.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.anaphase.region3.numTracks
                plot(dataTracks_new.anaphase.region3.tracksX(:,i),dataTracks_new.anaphase.region3.tracksY(:,i))
                hold all
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first 2 min anaphase
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region1.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.entireRecording.region1.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6  ) > dataTracks_new.entireRecording.region1.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region1.indexXStart(1,i) >= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 )
                    plot(dataTracks_new.entireRecording.region1.tracksX(:,i),dataTracks_new.entireRecording.region1.tracksY(:,i))
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_first-2min_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_first-2min_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region2.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.entireRecording.region2.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6  ) > dataTracks_new.entireRecording.region2.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region2.indexXStart(1,i) >= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 )
                    plot(dataTracks_new.entireRecording.region2.tracksX(:,i),dataTracks_new.entireRecording.region2.tracksY(:,i))
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_first-2min_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_first-2min_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region3.numTracks > 0
            
            figure
            for i = 1 : dataTracks_new.entireRecording.region3.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6  ) > dataTracks_new.entireRecording.region3.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region3.indexXStart(1,i) >= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 )
                    plot(dataTracks_new.entireRecording.region3.tracksX(:,i),dataTracks_new.entireRecording.region3.tracksY(:,i))
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_first-2min_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_first-2min_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % next 2 min in anaphse
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region1.numTracks > 0
            count = 0;
            figure
            for i = 1 : dataTracks_new.entireRecording.region1.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 240) *param.sp6  ) > dataTracks_new.entireRecording.region1.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region1.indexXStart(1,i) >= ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6 )
                    plot(dataTracks_new.entireRecording.region1.tracksX(:,i),dataTracks_new.entireRecording.region1.tracksY(:,i))
                    count = count +1;
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_last-2min_region1-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_last-2min_region1-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region2.numTracks > 0
            count = 0;
            figure
            for i = 1 : dataTracks_new.entireRecording.region2.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 240) *param.sp6  ) > dataTracks_new.entireRecording.region2.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region2.indexXStart(1,i) >= ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6 )
                    plot(dataTracks_new.entireRecording.region2.tracksX(:,i),dataTracks_new.entireRecording.region2.tracksY(:,i))
                    count = count +1;
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_last-2min_region2-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_last-2min_region2-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
    if ~isempty(find(contains('anaphase',conditions1)))
        if dataTracks_new.anaphase.region1.numTracks > 0
            count = 0;
            figure
            for i = 1 : dataTracks_new.entireRecording.region3.numTracks
                if  ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 240) *param.sp6  ) > dataTracks_new.entireRecording.region3.indexXStart(i) && ...
                        dataTracks_new.entireRecording.region3.indexXStart(1,i) >= ( furrow_position.image_start_detection - (param.delta_furrowDetection_anaphaseOnset - 120) *param.sp6 )
                    plot(dataTracks_new.entireRecording.region3.tracksX(:,i),dataTracks_new.entireRecording.region3.tracksY(:,i))
                    count = count +1;
                    hold all
                end
            end
            xlabel ('x coordinate');
            ylabel ('y coordinate');
            figureName = strcat('tracksChecking_anaphase_last-2min_region3-', short_name, '.fig');
            saveas(gcf,fullfile(main_path,figureName));
            figureName = strcat('tracksChecking_anaphase_last-2min_region3-', short_name, '.tif');
            saveas(gcf,fullfile(main_path,figureName));
            close(gcf)
            
        end
    end
    
end

    
end

