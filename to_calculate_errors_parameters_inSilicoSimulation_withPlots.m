function [ fitting_results ] = to_calculate_errors_parameters_inSilicoSimulation_withPlots...
    ( fitting_results,best_model,timing_phase,blastomere,nbEmbryo_givenCondition,save_stem_error,...
    tmin,fixed_short_lifetime,fixed_short_percent,plot_visualization,algo_fitting,nb_simulations )

% perform in silico simulations to get errors on parameters

if ~exist('plot_visualization','var')
    plot_visualization = 0;
end

global general_param

if nargin < 12
    nb_simulations = general_param.cortex_analysis.bootstap_nb_for_errors_determination;
end

models = {best_model};

%% allocate empty variables

if strcmp(best_model,'MonoExpo')   
    T_allSimulations = [];  
elseif strcmp(best_model,'DoubleExpo')   
    T1_allSimulations = [];
    T2_allSimulations = [];
    P1_allSimulations = [];
    P2_allSimulations = [];   
elseif strcmp(best_model,'DoubleExpo_fixedT0')  
    T2_2_allSimulations = [];
    P2_1_allSimulations = [];
    P2_2_allSimulations = [];   
elseif strcmp(best_model,'DoubleExpo_fixedT0P0')  
    T2_2_allSimulations_ = []; 
elseif strcmp(best_model,'TripleExpo')
    TT1_allSimulations = [];
    TT2_allSimulations = [];
    TT3_allSimulations = [];
    PP1_allSimulations = [];
    PP2_allSimulations = [];
    PP3_allSimulations = [];
elseif strcmp(best_model,'TripleExpo_fixedT0')  
    T1_allSimulations_ = [];
    T2_allSimulations_ = [];
    P0_allSimulations_ = [];
    P1_allSimulations_ = [];
    P2_allSimulations_ = [];   
elseif strcmp(best_model,'TripleExpo_fixedT0P0')  
    T1_allSimulations__ = [];
    T2_allSimulations__ = [];
    P1_allSimulations__ = [];
    P2_allSimulations__ = [];      
elseif strcmp(best_model,'MonoExpo_stretched')
    Ts_allSimulations = [];
    power_allSimulations = [];
elseif strcmp(best_model,'QuadroExpo')
    TTT1_allSimulations = [];
    TTT2_allSimulations = [];
    TTT3_allSimulations = [];
    TTT4_allSimulations = [];
    PPP1_allSimulations = [];
    PPP2_allSimulations = [];
    PPP3_allSimulations = [];
    PPP4_allSimulations = [];
end

%% get distributions + fitting part

for iSimulation = 1 : nb_simulations
    
    disp(['Parameter errors in progress: Bootsrap #' num2str(iSimulation)]);
    
    if nbEmbryo_givenCondition == 1
        
        data_simulation = to_get_inSilicoSets_of_tracksDurationHisto2_...
            ( fitting_results,models,blastomere,timing_phase,nbEmbryo_givenCondition,tmin,fixed_short_lifetime,fixed_short_percent );
        name_embryo = ['embryo' num2str(1)];
        data_singleSimulation = data_simulation.(name_embryo).data;
        size_population_ref = nansum(data_singleSimulation(:,2));
        while size_population_ref == 0
            data_simulation = to_get_inSilicoSets_of_tracksDurationHisto2_...
                ( fitting_results,models,blastomere,timing_phase,nbEmbryo_givenCondition,tmin,fixed_short_lifetime,fixed_short_percent );
            name_embryo = ['embryo' num2str(1)];
            data_singleSimulation = data_simulation.(name_embryo).data;
            size_population_ref = nansum(data_singleSimulation(:,2));
        end
        
        if algo_fitting == 1
            [ fitting_simulation ] =  to_calculate_parameters_using_fmincon_mle_...
                (data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent );
        elseif algo_fitting == 2
            [ fitting_simulation ] =  to_calculate_parameters_using_fminsearch_mle_...
                (data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent );
        elseif algo_fitting == 3
            [ fitting_simulation ] =  to_calculate_parameters_using_fminunc_mle_...
                (data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent );            
        elseif algo_fitting == 4 
            [ fitting_simulation] = to_calculate_parameters_using_fminsearch_fminunc_mle...
                (data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent );
        elseif algo_fitting == 5 
            %[ fitting_simulation] = to_calculate_parameters_using_fmincon_fminunc_mle...
                %(data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent ); 
                disp('not implemented yet');
        elseif algo_fitting == 6 
           % [ fitting_simulation] = to_calculate_parameters_using_fmincon_fminsearch_mle...
            %    (data_singleSimulation,models,timing_phase,blastomere,[],fixed_short_lifetime,fixed_short_percent );   
            disp('not implemented yet');
        end
        
        if strcmp(best_model,'MonoExpo')
            if fitting_simulation.(blastomere).(timing_phase).MonoExpo.flag > 0
                T_allSimulations = [ [T_allSimulations] fitting_simulation.(blastomere).(timing_phase).MonoExpo.T ];
            end
        elseif strcmp(best_model,'DoubleExpo')
            if fitting_simulation.(blastomere).(timing_phase).DoubleExpo.flag > 0
                T1_allSimulations = [ [T1_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo.T1 ];
                T2_allSimulations = [ [T2_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo.T2 ];
                P1_allSimulations = [ [P1_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo.P1 ];
                P2_allSimulations = [ [P2_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo.P2 ];
            end
        elseif strcmp(best_model,'DoubleExpo_fixedT0')
            if fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0.flag > 0
                T2_2_allSimulations = [ [T2_2_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0.T2 ];
                P2_1_allSimulations = [ [P2_1_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0.P1 ];
                P2_2_allSimulations = [ [P2_2_allSimulations] fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0.P2 ];
            end     
        elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
            if fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0P0.flag > 0
                T2_2_allSimulations_ = [ [T2_2_allSimulations_] fitting_simulation.(blastomere).(timing_phase).DoubleExpo_fixedT0P0.T2 ];
            end              
        elseif strcmp(best_model,'TripleExpo')
            if fitting_simulation.(blastomere).(timing_phase).TripleExpo.flag > 0
                TT1_allSimulations = [ [TT1_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.TT1 ];
                TT2_allSimulations = [ [TT2_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.TT2 ];
                TT3_allSimulations = [ [TT3_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.TT3 ];
                PP1_allSimulations = [ [PP1_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.PP1 ];
                PP2_allSimulations = [ [PP2_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.PP2 ];
                PP3_allSimulations = [ [PP3_allSimulations] fitting_simulation.(blastomere).(timing_phase).TripleExpo.PP3 ];
            end
        elseif strcmp(best_model,'TripleExpo_fixedT0')
            if fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.flag > 0
                T1_allSimulations_ = [ [T1_allSimulations_] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.T1 ];
                T2_allSimulations_ = [ [T2_allSimulations_] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.T2 ];
                P0_allSimulations_ = [ [P0_allSimulations_] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.P0 ];
                P1_allSimulations_ = [ [P1_allSimulations_] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.P1 ];
                P2_allSimulations_ = [ [P2_allSimulations_] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0.P2 ];
            end
        elseif strcmp(best_model,'TripleExpo_fixedT0P0')
            if fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0P0.flag > 0
                T1_allSimulations__ = [ [T1_allSimulations__] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0P0.T1 ];
                T2_allSimulations__ = [ [T2_allSimulations__] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0P0.T2 ];
                P1_allSimulations__ = [ [P1_allSimulations__] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0P0.P1 ];
                P2_allSimulations__ = [ [P2_allSimulations__] fitting_simulation.(blastomere).(timing_phase).TripleExpo_fixedT0P0.P2 ];
            end              
        elseif strcmp(best_model,'MonoExpo_stretched')
            if fitting_simulation.(blastomere).(timing_phase).MonoExpo_stretched.flag > 0
                Ts_allSimulations = [ [Ts_allSimulations] fitting_simulation.(blastomere).(timing_phase).MonoExpo_stretched.Ts ];
                power_allSimulations = [ [power_allSimulations] fitting_simulation.(blastomere).(timing_phase).MonoExpo_stretched.power ];
            end
         elseif strcmp(best_model,'QuadroExpo')
            if fitting_simulation.(blastomere).(timing_phase).QuadroExpo.flag > 0
                TTT1_allSimulations = [ [TTT1_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.TTT1 ];
                TTT2_allSimulations = [ [TTT2_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.TTT2 ];
                TTT3_allSimulations = [ [TTT3_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.TTT3 ];
                TTT4_allSimulations = [ [TTT4_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.TTT4 ];
                PPP1_allSimulations = [ [PPP1_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.PPP1 ];
                PPP2_allSimulations = [ [PPP2_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.PPP2 ];
                PPP3_allSimulations = [ [PPP3_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.PPP3 ];
                PPP4_allSimulations = [ [PPP4_allSimulations] fitting_simulation.(blastomere).(timing_phase).QuadroExpo.PPP4 ];
            end           
        end
        
        clear fitting_simulation
        clear data_simulation
        
    else
        
        data_simulation = to_get_inSilicoSets_of_tracksDurationHisto...
            ( fitting_results,models,nbEmbryo_givenCondition,tmin,fixed_short_lifetime,fixed_short_percent  );
        
        if algo_fitting == 1
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fmincon_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );
        elseif algo_fitting == 2
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fminsearch_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );
       elseif algo_fitting == 3
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fminunc_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );            
        elseif algo_fitting == 4 
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fminsearch_fminunc_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  ); 
        elseif algo_fitting == 5 
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fmincon_fminunc_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );           
       elseif algo_fitting == 6 
            [ data_simulation,fitting_simulation] = to_calculate_parameters_using_fmincon_fminsearch_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );   
        elseif algo_fitting == 7
            [ data_simulation,fitting_simulation]  = to_calculate_parameters_using_patternsearch_sum_mle...
                ( data_simulation,models,nbEmbryo_givenCondition,1,fixed_short_lifetime,fixed_short_percent  );
        end
        
        if strcmp(best_model,'MonoExpo')
            if fitting_simulation.MonoExpo.flag > 0
                T_allSimulations = [ [T_allSimulations] fitting_simulation.MonoExpo.T ];
            end
        elseif strcmp(best_model,'DoubleExpo')
            if fitting_simulation.DoubleExpo.flag > 0
                T1_allSimulations = [ [T1_allSimulations] fitting_simulation.DoubleExpo.T1 ];
                T2_allSimulations = [ [T2_allSimulations] fitting_simulation.DoubleExpo.T2 ];
                P1_allSimulations = [ [P1_allSimulations] fitting_simulation.DoubleExpo.P1 ];
                P2_allSimulations = [ [P2_allSimulations] fitting_simulation.DoubleExpo.P2 ];
            end
        elseif strcmp(best_model,'DoubleExpo_fixedT0')
            if fitting_simulation.DoubleExpo_fixedT0.flag > 0
                T2_2_allSimulations = [ [T2_2_allSimulations] fitting_simulation.DoubleExpo_fixedT0.T2 ];
                P2_1_allSimulations = [ [P2_1_allSimulations] fitting_simulation.DoubleExpo_fixedT0.P1 ];
                P2_2_allSimulations = [ [P2_2_allSimulations] fitting_simulation.DoubleExpo_fixedT0.P2 ];
            end     
        elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
            if fitting_simulation.DoubleExpo_fixedT0P0.flag > 0
                T2_2_allSimulations_ = [ [T2_2_allSimulations_] fitting_simulation.DoubleExpo_fixedT0P0.T2 ];
            end               
        elseif strcmp(best_model,'TripleExpo')
            if fitting_simulation.TripleExpo.flag > 0
                TT1_allSimulations = [ [TT1_allSimulations] fitting_simulation.TripleExpo.TT1 ];
                TT2_allSimulations = [ [TT2_allSimulations] fitting_simulation.TripleExpo.TT2 ];
                TT3_allSimulations = [ [TT3_allSimulations] fitting_simulation.TripleExpo.TT3 ];
                PP1_allSimulations = [ [PP1_allSimulations] fitting_simulation.TripleExpo.PP1 ];
                PP2_allSimulations = [ [PP2_allSimulations] fitting_simulation.TripleExpo.PP2 ];
                PP3_allSimulations = [ [PP3_allSimulations] fitting_simulation.TripleExpo.PP3 ];
            end
        elseif strcmp(best_model,'TripleExpo_fixedT0')
            if fitting_simulation.TripleExpo_fixedT0.flag > 0
                T1_allSimulations_ = [ [T1_allSimulations_] fitting_simulation.TripleExpo_fixedT0.T1 ];
                T2_allSimulations_ = [ [T2_allSimulations_] fitting_simulation.TripleExpo_fixedT0.T2 ];
                P0_allSimulations_ = [ [P0_allSimulations_] fitting_simulation.TripleExpo_fixedT0.P0 ];
                P1_allSimulations_ = [ [P1_allSimulations_] fitting_simulation.TripleExpo_fixedT0.P1 ];
                P2_allSimulations_ = [ [P2_allSimulations_] fitting_simulation.TripleExpo_fixedT0.P2 ];
            end
        elseif strcmp(best_model,'TripleExpo_fixedT0P0')
            if fitting_simulation.TripleExpo_fixedT0P0.flag > 0
                T1_allSimulations__ = [ [T1_allSimulations__] fitting_simulation.TripleExpo_fixedT0P0.T1 ];
                T2_allSimulations__ = [ [T2_allSimulations__] fitting_simulation.TripleExpo_fixedT0P0.T2 ];
                P1_allSimulations__ = [ [P1_allSimulations__] fitting_simulation.TripleExpo_fixedT0P0.P1 ];
                P2_allSimulations__ = [ [P2_allSimulations__] fitting_simulation.TripleExpo_fixedT0P0.P2 ];    
            end            
        elseif strcmp(best_model,'MonoExpo_stretched')
            if fitting_simulation.MonoExpo_stretched.flag > 0
                Ts_allSimulations = [ [Ts_allSimulations] fitting_simulation.MonoExpo_stretched.Ts ];
                power_allSimulations = [ [power_allSimulations] fitting_simulation.MonoExpo_stretched.power ];
            end
        elseif strcmp(best_model,'QuadroExpo')
            if fitting_simulation.QuadroExpo.flag > 0
                TTT1_allSimulations = [ [TTT1_allSimulations] fitting_simulation.QuadroExpo.TTT1 ];
                TTT2_allSimulations = [ [TTT2_allSimulations] fitting_simulation.QuadroExpo.TTT2 ];
                TTT3_allSimulations = [ [TTT3_allSimulations] fitting_simulation.QuadroExpo.TTT3 ];
                TTT4_allSimulations = [ [TTT4_allSimulations] fitting_simulation.QuadroExpo.TTT4 ];
                PPP1_allSimulations = [ [PPP1_allSimulations] fitting_simulation.QuadroExpo.PPP1 ];
                PPP2_allSimulations = [ [PPP2_allSimulations] fitting_simulation.QuadroExpo.PPP2 ];
                PPP3_allSimulations = [ [PPP3_allSimulations] fitting_simulation.QuadroExpo.PPP3 ];
                PPP4_allSimulations = [ [PPP4_allSimulations] fitting_simulation.QuadroExpo.PPP4 ];
            end            
        end
        clear fitting_simulation
        clear data_simulation
    end
    
end


%% display plots to visualize/undertsnad why errors huge sometimes

if plot_visualization == 1
    if strcmp(best_model,'MonoExpo')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,T_allSimulations,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],nbEmbryo_givenCondition,save_stem_error,...
            best_model,50,1);
    elseif strcmp(best_model,'DoubleExpo')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,[],T1_allSimulations,T2_allSimulations,P1_allSimulations,P2_allSimulations,[],[],[],[],[],[],[],[],[],[],...
            [],[],nbEmbryo_givenCondition,save_stem_error,best_model,50,1);
    elseif strcmp(best_model,'DoubleExpo_fixedT0')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,[],[],[],[],[],T2_2_allSimulations,P2_1_allSimulations,P2_2_allSimulations,[],[],[],[],[],[],[],...
            [],[],nbEmbryo_givenCondition,save_stem_error,best_model,50,1);   
    elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,[],[],[],[],[],[],[],[],T2_2_allSimulations_,[],[],[],[],[],[],...
            [],[],nbEmbryo_givenCondition,save_stem_error,best_model,50,1);          
    elseif strcmp(best_model,'TripleExpo')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,[],[],[],[],[],[],[],[],[],TT1_allSimulations,TT2_allSimulations,TT3_allSimulations,PP1_allSimulations,...
            PP2_allSimulations,PP3_allSimulations,[],[],nbEmbryo_givenCondition,save_stem_error,best_model,50,1);
    elseif strcmp(best_model,'MonoExpo_stretched')
        to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],Ts_allSimulations,power_allSimulations,...
            nbEmbryo_givenCondition,save_stem_error,best_model,50,1);
    end
end

%% get the errors as the std or sem of the obtained fitting values set
%% + study an eventual correlation between T1/T2 or P1/P2 or T1/P1 or P2/T2

if nbEmbryo_givenCondition == 1
    
    if strcmp(best_model,'MonoExpo')
        
        fitting_results.(blastomere).(timing_phase).MonoExpo.T_se = std(T_allSimulations);
        
    elseif strcmp(best_model,'DoubleExpo')
        
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T1_se = std(T1_allSimulations);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T2_se = std(T2_allSimulations);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P1_se = std(P1_allSimulations);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P2_se = std(P2_allSimulations);
        
        %-----------------
        % calculate correlation
        
        R2_T1T2 = corr(T1_allSimulations(:), T2_allSimulations(:));
        R2_T1P1 = corr(T1_allSimulations(:), P1_allSimulations(:));
        R2_T2P2 = corr(T2_allSimulations(:), P2_allSimulations(:));
        
        for j_ = 1 : length(T1_allSimulations)
            ratio_T(j_) = T1_allSimulations(j_) ./ T2_allSimulations(j_) ;
        end
        R2_ratios = corr(ratio_T(:), P1_allSimulations(:));
        
        if plot_visualization == 1
            
            %------------------
            % plot correlation between parameters
            
            figure
            plot(T1_allSimulations, T2_allSimulations, 'o');
            xlabel('tracks duration 1 (s)','FontSize',10);
            ylabel('tracks duration 2 (s)','FontSize',10);
            title(['Correlation T1/T2 :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_T1T2,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T1T2-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T1_allSimulations, P1_allSimulations, 'o');
            xlabel('tracks duration 1 (s)','FontSize',10);
            ylabel('ratio pop 1 (a.u.)','FontSize',10);
            title(['Correlation T1/P1 :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_T1P1,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T1P1-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T2_allSimulations, P2_allSimulations, 'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 2 (a.u.)','FontSize',10);
            title(['Correlation T2/P2 :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_T2P2,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P2-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(ratio_T, P1_allSimulations,'o');
            xlabel('Duration ratio (a.u.)','FontSize',10);
            ylabel('Pop ratio 2 (a.u.)','FontSize',10);
            title(['Correlation ratios :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_ratios,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_ratios-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end
        
    elseif strcmp(best_model,'DoubleExpo_fixedT0')
        
        fitting_results.(blastomere).(timing_phase).DoubleExpo_fixedT0.T2_se = 1.96 * std(T2_2_allSimulations);
        fitting_results.(blastomere).(timing_phase).DoubleExpo_fixedT0.P1_se = 1.96 * std(P2_1_allSimulations);
        fitting_results.(blastomere).(timing_phase).DoubleExpo_fixedT0.P2_se = 1.96 * std(P2_2_allSimulations);
        
        %-----------------
        % calculate correlation
        
        R2_T2P2_ = corr(T2_2_allSimulations(:), P2_2_allSimulations(:));
        R2_T2P1_ = corr(T2_2_allSimulations(:), P2_1_allSimulations(:));        
        
        if plot_visualization == 1
            
            %------------------
            % plot correlation between parameters
            
            figure
            plot(T2_2_allSimulations, P2_1_allSimulations, 'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 1 (a.u.)','FontSize',10);
            title(['Correlation T2/P1 :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_T2P1_,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P1_fixedT1-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T2_2_allSimulations, P2_2_allSimulations, 'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 2 (a.u.)','FontSize',10);
            title(['Correlation T2/P2 :' blastomere ' and ' timing_phase ]);
            string = ['R2 = ' num2str(round2(R2_T2P2_,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P2_fiexdT1-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end        
        
    elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
        fitting_results.(blastomere).(timing_phase).DoubleExpo_fixedT0P0.T2_se = 1.96 * std(T2_2_allSimulations_);
            
    elseif strcmp(best_model,'TripleExpo')
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1_se = std(TT1_allSimulations);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2_se = std(TT2_allSimulations);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3_se = std(TT3_allSimulations);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1_se = std(PP1_allSimulations);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2_se = std(PP2_allSimulations);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3_se = std(PP3_allSimulations);
        
    elseif strcmp(best_model,'TripleExpo_fixedT0')
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T1_se = 1.96 * std(T1_allSimulations_);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T2_se = 1.96 * std(T2_allSimulations_);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P0_se = 1.96 * std(P0_allSimulations_);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P1_se = 1.96 * std(P1_allSimulations_);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P2_se = 1.96 * std(P2_allSimulations_);
        
    elseif strcmp(best_model,'TripleExpo_fixedT0P0')
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0P0.T1_se = 1.96 * std(T1_allSimulations__);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0P0.T2_se = 1.96 * std(T2_allSimulations__);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0P0.P1_se = 1.96 * std(P1_allSimulations__);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0P0.P2_se = 1.96 * std(P2_allSimulations__);         
        
    elseif strcmp(best_model,'MonoExpo_stretched')
        fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.T_se = std(Ts_allSimulations);
        fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.power_se = std(power_allSimulations);
        
    elseif strcmp(best_model,'QuadroExpo')
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT1_se = std(TTT1_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT2_se = std(TTT2_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT3_se = std(TTT3_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT4_se = std(TTT4_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP1_se = std(PPP1_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP2_se = std(PPP2_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP3_se = std(PPP3_allSimulations);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP4_se = std(PPP4_allSimulations);
        
    end
    
else
    
    if strcmp(best_model,'MonoExpo')
        
        fitting_results.MonoExpo.T_se = std(T_allSimulations);
        
    elseif strcmp(best_model,'DoubleExpo')
        
        fitting_results.DoubleExpo.T1_se = std(T1_allSimulations);
        fitting_results.DoubleExpo.T2_se = std(T2_allSimulations);
        fitting_results.DoubleExpo.P1_se = std(P1_allSimulations);
        fitting_results.DoubleExpo.P2_se = std(P2_allSimulations);
        
        %-----------------
        % calculate correlation
        
        R2_T1T2 = corr(T1_allSimulations(:), T2_allSimulations(:));
        R2_T1P1 = corr(T1_allSimulations(:), P1_allSimulations(:));
        R2_T2P2 = corr(T2_allSimulations(:), P2_allSimulations(:));
        
        for j_ = 1 : length(T1_allSimulations)
            ratio_T(j_) = T1_allSimulations(j_) ./ T2_allSimulations(j_) ;
        end
        R2_ratios = corr(ratio_T(:), P1_allSimulations(:));
        
        
        if plot_visualization == 1
            %------------------
            % plot correlation between parameters
            
            figure
            plot(T1_allSimulations, T2_allSimulations,'o');
            xlabel('tracks duration 1 (s)','FontSize',10);
            ylabel('tracks duration 2 (s)','FontSize',10);
            title(['Correlation T1/T2 :' ]);
            string = ['R2 = ' num2str(round2(R2_T1T2,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T1T2.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T1_allSimulations, P1_allSimulations,'o');
            xlabel('tracks duration 1 (s)','FontSize',10);
            ylabel('ratio pop 1 (a.u.)','FontSize',10);
            title(['Correlation T1/P1 :' ]);
            string = ['R2 = ' num2str(round2(R2_T1P1,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T1P1.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T2_allSimulations, P2_allSimulations,'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 2 (a.u.)','FontSize',10);
            title(['Correlation T2/P2 :' ]);
            string = ['R2 = ' num2str(round2(R2_T2P2,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P2.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(ratio_T, P1_allSimulations,'o');
            xlabel('Duration ratio (a.u.)','FontSize',10);
            ylabel('Pop ratio 1 (a.u.)','FontSize',10);
            title(['Correlation ratios :']);
            string = ['R2 = ' num2str(round2(R2_ratios,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_ratios.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end
         
        
    elseif strcmp(best_model,'DoubleExpo_fixedT0')
        
        fitting_results.DoubleExpo_fixedT0.T2_se = 1.96 * std(T2_2_allSimulations);
        fitting_results.DoubleExpo_fixedT0.P1_se = 1.96 * std(P2_1_allSimulations);
        fitting_results.DoubleExpo_fixedT0.P2_se = 1.96 * std(P2_2_allSimulations);
        
        %-----------------
        % calculate correlation
        
        R2_T2P2_ = corr(T2_2_allSimulations(:), P2_2_allSimulations(:));
        R2_T2P1_ = corr(T2_2_allSimulations(:), P2_1_allSimulations(:));  
        
        
        if plot_visualization == 1
            
            %------------------
            % plot correlation between parameters
            
            figure
            plot(T2_2_allSimulations, P2_1_allSimulations, 'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 1 (a.u.)','FontSize',10);
            title(['Correlation T2/P1 ' ]);
            string = ['R2 = ' num2str(round2(R2_T2P1_,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P1_fixedT1-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            plot(T2_2_allSimulations, P2_2_allSimulations, 'o');
            xlabel('tracks duration 2 (s)','FontSize',10);
            ylabel('ratio pop 2 (a.u.)','FontSize',10);
            title(['Correlation T2/P2 ' ]);
            string = ['R2 = ' num2str(round2(R2_T2P2_,1e-2)) ];
            text(25,50,string,'Units','pixels')
            namePlot = strcat('Correlation_T2P2_fiexdT1-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end              
    elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
        fitting_results.DoubleExpo_fixedT0P0.T2_se = 1.96 * std(T2_2_allSimulations_);
        
    elseif strcmp(best_model,'TripleExpo')
        fitting_results.TripleExpo.TT1_se = std(TT1_allSimulations);
        fitting_results.TripleExpo.TT2_se = std(TT2_allSimulations);
        fitting_results.TripleExpo.TT3_se = std(TT3_allSimulations);
        fitting_results.TripleExpo.PP1_se = std(PP1_allSimulations);
        fitting_results.TripleExpo.PP2_se = std(PP2_allSimulations);
        fitting_results.TripleExpo.PP3_se = std(PP3_allSimulations);
        
    elseif strcmp(best_model,'TripleExpo_fixedT0')
        fitting_results.TripleExpo_fixedT0.T1_se = 1.96 * std(T1_allSimulations_);
        fitting_results.TripleExpo_fixedT0.T2_se = 1.96 * std(T2_allSimulations_);
        fitting_results.TripleExpo_fixedT0.P0_se = 1.96 * std(P0_allSimulations_);
        fitting_results.TripleExpo_fixedT0.P1_se = 1.96 * std(P1_allSimulations_);
        fitting_results.TripleExpo_fixedT0.P2_se = 1.96 * std(P2_allSimulations_);
        
    elseif strcmp(best_model,'TripleExpo_fixedT0P0')
        fitting_results.TripleExpo_fixedT0P0.T1_se = 1.96 * std(T1_allSimulations__);
        fitting_results.TripleExpo_fixedT0P0.T2_se = 1.96 * std(T2_allSimulations__);
        fitting_results.TripleExpo_fixedT0P0.P1_se = 1.96 * std(P1_allSimulations__);
        fitting_results.TripleExpo_fixedT0P0.P2_se = 1.96 * std(P2_allSimulations__);           
        
    elseif strcmp(best_model,'MonoExpo_stretched')
        fitting_results.MonoExpo_stretched.T_se = std(Ts_allSimulations);
        fitting_results.MonoExpo_stretched.power_se = std(power_allSimulations);
        
    elseif strcmp(best_model,'QuadroExpo')
        fitting_results.QuadroExpo.TTT1_se = std(TTT1_allSimulations);
        fitting_results.QuadroExpo.TTT2_se = std(TTT2_allSimulations);
        fitting_results.QuadroExpo.TTT3_se = std(TTT3_allSimulations);
        fitting_results.QuadroExpo.TTT4_se = std(TTT4_allSimulations);
        fitting_results.QuadroExpo.PPP1_se = std(PPP1_allSimulations);
        fitting_results.QuadroExpo.PPP2_se = std(PPP2_allSimulations);
        fitting_results.QuadroExpo.PPP3_se = std(PPP3_allSimulations);
        fitting_results.QuadroExpo.PPP4_se = std(PPP4_allSimulations);
         
    end
    
end


end

