function [ familyAssignement,dataTracks,tracks_duration_histo ] = to_assign_MTs_to_subPopulations_shortOrLong...
    ( dataTracks,best_model,fitting_results,name1,name2,familyAssignement,proba_threshold,tracks_duration_histo,iEmbryo,main_path,frequency,given_number )

global param
global general_param

if nargin < 11
    frequency = param.sp6;
end


%% perform assignment

if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
    
    [ familyAssignement,probability_long_duration,probability_short_duration,probability_veryshort_duration,fitting_results ] = ...
        to_assign_MTs_toGivenPopulation_usingProba_mle_tripleExpo_...
        ( dataTracks,fitting_results,name1,name2,familyAssignement,proba_threshold,given_number,frequency );
    
elseif strcmp(best_model,'DoubleExpo') || strcmp(best_model,'DoubleExpo_fixedT0')  || strcmp(best_model,'DoubleExpo_fixedT0P0')
    
    [ familyAssignement,probability_long_duration,probability_short_duration,fitting_results ] = ...
        to_assign_MTs_toGivenPopulation_usingProba_mle_doubleExpo__...
        ( dataTracks,fitting_results,name1,name2,familyAssignement,proba_threshold,given_number,frequency  );
    
else
    
    familyAssignement.short_duration.number = 0;
    familyAssignement.long_duration.number = 0;
    
end


%% saving in dataTracks

% very-short duration
if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
    dataTracks.(name1).(name2).veryshort_duration.lengthTracks = ...
        dataTracks.(name1).(name2).lengthTracks(familyAssignement.veryshort_duration.index);
    dataTracks.(name1).(name2).veryshort_duration.indexXStart = ...
        dataTracks.(name1).(name2).indexXStart(familyAssignement.veryshort_duration.index);
    dataTracks.(name1).(name2).veryshort_duration.indexXEnd = ...
        dataTracks.(name1).(name2).indexXEnd(familyAssignement.veryshort_duration.index);
    dataTracks.(name1).(name2).veryshort_duration.numTracks = familyAssignement.veryshort_duration.number;
    dataTracks.(name1).(name2).veryshort_duration.tracksX = ...
        dataTracks.(name1).(name2).tracksX(:, familyAssignement.veryshort_duration.index);
end

% short duration
if familyAssignement.short_duration.number > 0
    dataTracks.(name1).(name2).short_duration.lengthTracks = ...
        dataTracks.(name1).(name2).lengthTracks(familyAssignement.short_duration.index);
    dataTracks.(name1).(name2).short_duration.indexXStart = ...
        dataTracks.(name1).(name2).indexXStart(familyAssignement.short_duration.index);
    dataTracks.(name1).(name2).short_duration.indexXEnd = ...
        dataTracks.(name1).(name2).indexXEnd(familyAssignement.short_duration.index);
    dataTracks.(name1).(name2).short_duration.numTracks = familyAssignement.short_duration.number;
    dataTracks.(name1).(name2).short_duration.tracksX = ...
        dataTracks.(name1).(name2).tracksX(:, familyAssignement.short_duration.index);
else
    dataTracks.(name1).(name2).short_duration.numTracks = familyAssignement.short_duration.number;
    dataTracks.(name1).(name2).short_duration.lengthTracks = [];
    dataTracks.(name1).(name2).short_duration.indexXStart = [];
    dataTracks.(name1).(name2).short_duration.indexXEnd = [];
    dataTracks.(name1).(name2).short_duration.tracksX = [];
end

% long duration
if familyAssignement.long_duration.number > 0
    dataTracks.(name1).(name2).long_duration.lengthTracks = ...
        dataTracks.(name1).(name2).lengthTracks(familyAssignement.long_duration.index);
    dataTracks.(name1).(name2).long_duration.indexXStart = ...
        dataTracks.(name1).(name2).indexXStart(familyAssignement.long_duration.index);
    dataTracks.(name1).(name2).long_duration.indexXEnd = ...
        dataTracks.(name1).(name2).indexXEnd(familyAssignement.long_duration.index);
    dataTracks.(name1).(name2).long_duration.numTracks = familyAssignement.long_duration.number;
    dataTracks.(name1).(name2).long_duration.tracksX = ...
        dataTracks.(name1).(name2).tracksX(:, familyAssignement.long_duration.index);
else
    dataTracks.(name1).(name2).long_duration.numTracks = familyAssignement.long_duration.number;
    dataTracks.(name1).(name2).long_duration.lengthTracks = [];
    dataTracks.(name1).(name2).long_duration.indexXStart = [];
    dataTracks.(name1).(name2).long_duration.indexXEnd = [];
    dataTracks.(name1).(name2).long_duration.tracksX = [];
end

% tag as all tracks the ones not assigned
dataTracks.(name1).(name2).notAssigned.numTracks = dataTracks.(name1).(name2).numTracks;
dataTracks.(name1).(name2).notAssigned.lengthTracks = dataTracks.(name1).(name2).lengthTracks;
dataTracks.(name1).(name2).notAssigned.indexXStart = dataTracks.(name1).(name2).indexXStart;
dataTracks.(name1).(name2).notAssigned.indexXEnd = dataTracks.(name1).(name2).indexXEnd;
dataTracks.(name1).(name2).notAssigned.tracksX = dataTracks.(name1).(name2).tracksX;


if ~isempty(tracks_duration_histo)
    
%% get histograms
    
    if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
        local_conditions = {'veryshort_duration' 'short_duration' 'long_duration'};
    elseif strcmp(best_model,'DoubleExpo') || strcmp(best_model,'DoubleExpo_fixedT0')  || strcmp(best_model,'DoubleExpo_fixedT0P0')
        local_conditions = {'short_duration' 'long_duration'};
    end
    
    for iSet = 1 : numel(local_conditions)
        
        given_set = local_conditions{iSet};
        
        if familyAssignement.(given_set).number > 0
            
            maxLength= max(dataTracks.(name1).(name2).(given_set).lengthTracks);
            minLength= min(dataTracks.(name1).(name2).(given_set).lengthTracks);
            binranges = [minLength : 1 : maxLength];
            
            if maxLength > 100 && minLength < 100
                index_to_remove = [(100-minLength+1):(maxLength-minLength+1)];
                binranges(index_to_remove) = [];
            end
            
            [bincounts] = hist(dataTracks.(name1).(name2).(given_set).lengthTracks,binranges);
            
            if length(binranges) == length(bincounts)
                binranges(bincounts==0) = [];
                bincounts(bincounts==0) = [];
            else
                bincounts(bincounts==0) = [];
            end
            
            if  strcmp(given_set,'long_duration')
                % removed too much count equal to 1
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
                    binranges(max_index_one) = max_index_one + ( max_max_index_one - max_index_one)/2;
                    bincounts(max_index_one) = bincounts_above_max_one;
                end
            end
            
            if strcmp(given_set,'short_duration') || strcmp(given_set,'long_duration')
                
                % to perform smoothing of the distribution
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
            end
            
            % check for any inf in data
            [ ri] = find(~isfinite(bincounts));
            bincounts(ri)= [];
            binranges(ri)= [];
            
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).bincounts = bincounts;
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).binranges = binranges./param.sp6; % in sec
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).number_tracks = sum(bincounts);
            
            clear maxLength
            clear minimLength
            clear bincounts
            clear binranges
            
        else
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).number_tracks = 0;
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).binranges = []; % in sec
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).bincounts = [];
            tracks_duration_histo{iEmbryo}.(name1).(name2).(given_set).number_tracks = [];
            
        end
        
    end
    
    
    %% to plot results of the assignement
    
    cell_legend = {};
    local_count = 0;
    
    figure
    
    %all count
    plot(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges,tracks_duration_histo{iEmbryo}.(name1).(name2).bincounts, 'ok');
    hold all
    local_count = local_count +1;
    cell_legend{local_count} = 'Raw distribution';
    % count of each sub-population
    if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
        if tracks_duration_histo{iEmbryo}.(name1).(name2).veryshort_duration.number_tracks > 0
            plot(tracks_duration_histo{iEmbryo}.(name1).(name2).veryshort_duration.binranges,...
                tracks_duration_histo{iEmbryo}.(name1).(name2).veryshort_duration.bincounts, 'or');
            local_count = local_count +1;
            cell_legend{local_count} = 'veryshort duration pop. count';
        end
    end
    if tracks_duration_histo{iEmbryo}.(name1).(name2).short_duration.number_tracks > 0
        plot(tracks_duration_histo{iEmbryo}.(name1).(name2).short_duration.binranges,...
            tracks_duration_histo{iEmbryo}.(name1).(name2).short_duration.bincounts, 'og');
        local_count = local_count +1;
        cell_legend{local_count} = 'short duration pop. count';
    end
    if tracks_duration_histo{iEmbryo}.(name1).(name2).long_duration.number_tracks > 0
        plot(tracks_duration_histo{iEmbryo}.(name1).(name2).long_duration.binranges,...
            tracks_duration_histo{iEmbryo}.(name1).(name2).long_duration.bincounts, 'ob');
        local_count = local_count +1;
        cell_legend{local_count} = 'long duration pop. count';
    end
    % function of best model
    if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
        if iEmbryo == 0
            f_TripleExpo = @(x) fitting_results.size_population_raw  * ...
                ( fitting_results.TripleExpo.PP1 / fitting_results.TripleExpo.R_normalization.triple(1,1) *  exp(-x/fitting_results.TripleExpo.TT1)+ ...
                fitting_results.TripleExpo.PP2 / fitting_results.TripleExpo.R_normalization.triple(2,1) * exp(-x/fitting_results.TripleExpo.TT2)+ ...
                fitting_results.TripleExpo.PP3 / fitting_results.TripleExpo.R_normalization.triple(3,1) * exp(-x/fitting_results.TripleExpo.TT3)  );
        else
            name0 = ['embryo', num2str(iEmbryo)];
            f_TripleExpo = @(x) fitting_results.size_population.(name0).raw * ...
                ( fitting_results.TripleExpo.PP1 / fitting_results.TripleExpo.R_normalization.(name0).triple(1,1) *  exp(-x/fitting_results.TripleExpo.TT1)+ ...
                fitting_results.TripleExpo.PP2 / fitting_results.TripleExpo.R_normalization.(name0).triple(2,1) * exp(-x/fitting_results.TripleExpo.TT2)+ ...
                fitting_results.TripleExpo.PP3 / fitting_results.TripleExpo.R_normalization.(name0).triple(3,1) * exp(-x/fitting_results.TripleExpo.TT3)  );
        end
        fplot(f_TripleExpo,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-k');
        local_count = local_count +1;
        cell_legend{local_count} = 'Triple expo. fit';
    elseif strcmp(best_model,'DoubleExpo') || strcmp(best_model,'DoubleExpo_fixedT0')  || strcmp(best_model,'DoubleExpo_fixedT0P0')
        if iEmbryo == 0
            f_DoubleExpo = @(x) fitting_results.size_population_raw * ...
                ( fitting_results.DoubleExpo.P1 / fitting_results.DoubleExpo.R_normalization.double(1,1) *  exp(-x/fitting_results.DoubleExpo.T1)+ ...
                fitting_results.DoubleExpo.P2 / fitting_results.DoubleExpo.R_normalization.double(2,1) * exp(-x/fitting_results.DoubleExpo.T2) );
            fplot(f_DoubleExpo,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-k');
        else
            name0 = ['embryo', num2str(iEmbryo)];
            f_DoubleExpo = @(x) fitting_results.size_population.(name0).raw * ...
                ( fitting_results.DoubleExpo.P1 / fitting_results.DoubleExpo.R_normalization.(name0).double(1,1) *  exp(-x/fitting_results.DoubleExpo.T1)+ ...
                fitting_results.DoubleExpo.P2 / fitting_results.DoubleExpo.R_normalization.(name0).double(2,1) * exp(-x/fitting_results.DoubleExpo.T2) );
        end
        fplot(f_DoubleExpo,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-k');
        local_count = local_count +1;
        cell_legend{local_count} = 'Double expo. fit';
    end
    % fonction of sub-populations
    if strcmp(best_model,'TripleExpo') || strcmp(best_model,'TripleExpo_fixedT0')  || strcmp(best_model,'TripleExpo_fixedT0P0')
        if familyAssignement.veryshort_duration.number > 0
            fplot(probability_veryshort_duration,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-r');
            local_count = local_count +1;
            cell_legend{local_count} = 'very-short-lived pop. fit';
        end
    end
    if familyAssignement.short_duration.number > 0
        fplot(probability_short_duration,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-g');
        local_count = local_count +1;
        cell_legend{local_count} = 'short-lived pop. fit';
    end
    
    if familyAssignement.long_duration.number > 0
        fplot(probability_long_duration,[0 max(tracks_duration_histo{iEmbryo}.(name1).(name2).binranges)], '-b');
        local_count = local_count +1;
        cell_legend{local_count} = 'long-lived pop. fit';
    end
    % layout
    ylim([1 max(tracks_duration_histo{iEmbryo}.(name1).(name2).bincounts)+500])
    set(gca,'yscale','log');
    xlabel('Track duration (s)')
    xlabel('Count (a.u.)')
    legend(cell_legend);
    name0 =  [ 'embryo', num2str(iEmbryo) ];
    namePlot = strcat('Count-distribution_track_duration_distribution_subpopulations_proba_',...
        num2str(proba_threshold), '_', name0, '_', name1, '_', name2, '.fig');
    saveas(gcf,fullfile(main_path, namePlot));
    close(gcf)
    
end

end

