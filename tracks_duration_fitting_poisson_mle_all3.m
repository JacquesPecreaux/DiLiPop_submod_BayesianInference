
function [fitting_results,model_choice] = tracks_duration_fitting_poisson_mle_all3...
    (models,name1,name2,nbEmbryo_givenCondition,tracks_duration_histo,index_embryo_metaphase,save_stem,...
    use_loglikelihood_to_calculate_parameters,index_embryo_late,classification_performed,algo_choice,error_computation,...
    choice_error_estimate,index_embryo_prophase,index_embryo_prometaphase,index_embryo_anaphase)

global general_param

if nargin < 10
    classification_performed = 0;
end

if nargin < 11
    algo_choice = 1; %fmincon = 1
end

if nargin < 12
    error_computation = 1;
end

if nargin < 13
    choice_error_estimate = 1;
end

for iiEmbryo = 1:nbEmbryo_givenCondition
    
    if strcmp(name1,'metaphase')
        iEmbryo = index_embryo_metaphase(iiEmbryo);
    elseif strcmp(name1,'late')
        iEmbryo = index_embryo_late(iiEmbryo);
    elseif strcmp(name1,'prometaphase')
        iEmbryo = index_embryo_prometaphase(iiEmbryo);
    elseif strcmp(name1,'prophase')
        iEmbryo = index_embryo_prophase(iiEmbryo);
    elseif strcmp(name1,'anaphase')
        iEmbryo = index_embryo_anaphase(iiEmbryo);
    else
        iEmbryo = iiEmbryo;
    end
    
    name_embryo = ['embryo' num2str(iiEmbryo)];
    
    input.(name_embryo).data = double( transpose( cat(1, tracks_duration_histo{iEmbryo}.(name1).(name2).binranges, ...
        tracks_duration_histo{iEmbryo}.(name1).(name2).bincounts ) ) );
    input.(name_embryo).duration_phase = tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase;
    input.(name_embryo).area = tracks_duration_histo{iEmbryo}.(name1).(name2).area; % area in squared um
    input.(name_embryo).raw_data = tracks_duration_histo{iEmbryo}.(name1).(name2).lengthTracks;
    
    if error_computation == 1 && choice_error_estimate == 2
            input.(name_embryo).raw_data = tracks_duration_histo{iEmbryo}.(name1).(name2).lengthTracks ;
    end
    
end

%% fitting using fmincon function to get models parameters
% possibility to use jacobian matrix as well as hesian matix to have a
% better fit
% maybe also add condition

if algo_choice == 2
    [ input,fitting_results] = to_calculate_parameters_using_patternsearch_sum_mle...
        ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
elseif algo_choice == 1
    [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
        ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );
end
%[ input,fitting_results] = to_calculate_parameters_using_fminsearch_sum_mle...
%    ( input,models,nbEmbryo_givenCondition,0,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent );

%% other possibility to calculate model parameters based on the min of loglikelihood

if use_loglikelihood_to_calculate_parameters == 1
    [ fitting_results ] = to_calculate_parameters_using_min_loglikelihood_sumEmbryo...
        ( input,fitting_results,models,nbEmbryo_givenCondition,name1,name2,save_stem);
end

%% find best model from the unique parameters set obtained from fitting all embryos individually

[ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models );


%% plot the residues for each embryo in comparison to the best theoretical model

%to_plot_residuals_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models,best_model,name1,name2,save_stem,classification_performed );


%% plot the fit for all the models with the raw data of each embryo

% no log scale
%to_plot_data_withAllFits_sum_mle( input,fitting_results,models,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,0 );

% x log scale
%to_plot_data_withAllFits_sum_mle( input,fitting_results,models,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,1 );

% y log scale
to_plot_data_withAllFits_sum_mle( input,fitting_results,models,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );


%% two possibilities to evaluate the errors on the fitting parameters

% 1st possibility by use of Camer Rao bound, problem give a minimal error

%  P(A,T/Ci) = PI mi ^ ci exp(-mi) /ci!;
% -2ln(L) = 2 sum_i (mi-ci*ln(mi)) +cste

% variance of any unbiaised estimator is bounded by the reciprocal of the fisher information I(par)
% var(estimated par) > = 1/I(par)
% Fisher information : I(par) = - E[ d² l(x,par) / d par²]
% l(x,par) is the natural log of the likelihood function and E denotes teh
% expected value over x

%[ fitting_results ] = to_calculate_errors_parameters_using_CamerRaoBound_all_new2( input,fitting_results,models,nbEmbryo_givenCondition );
%[ fitting_results ] = to_calculate_errors_parameters_using_CamerRaoBound_all( input,fitting_results,models,nbEmbryo_givenCondition );

%------------------------------------------------------------------------------------------
% 2nd possibility and third possibilities below by using respectively
% likelihood profile and boostrapping

if error_computation == 1
    if choice_error_estimate == 1
        if strcmp(best_model,'MonoExpo_stretched') == 0 && strcmp(best_model,'TripleExpo') == 0
            [ fitting_results ] = to_calculate_errors_parameters_using_ProfileLikelihood...
                ( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,fval_mono,fval_double,fval_triple );
        elseif strcmp(best_model,'TripleExpo') == 1
            fitting_results.TripleExpo.TT1_se = NaN;
            fitting_results.TripleExpo.TT2_se = NaN;
            fitting_results.TripleExpo.TT3_se = NaN;
            fitting_results.TripleExpo.PP1_se = NaN;
            fitting_results.TripleExpo.PP2_se = NaN;
            fitting_results.TripleExpo.PP3_se = NaN;
        end
    elseif choice_error_estimate == 2
        if strcmp(best_model,'MonoExpo_stretched') == 0
            plot_visualization = 1;
            tmin = general_param.cortex_analysis.minLength;
            [ fitting_results ] = to_calculate_errors_parameters_using_Boostrapping...
                ( input,models,fitting_results,best_model,name1,name2,nbEmbryo_givenCondition,...
                tmin,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,plot_visualization );
        end
    elseif choice_error_estimate == 3
        if strcmp(best_model,'MonoExpo_stretched') == 0
            [ fitting_results ] = to_calculate_errors_parameters_using_ProfileLikelihood...
                ( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,fval_mono,fval_double,fval_triple );
            plot_visualization = 1;
            tmin = general_param.cortex_analysis.minLength;
            [ fitting_results ] = to_calculate_errors_parameters_using_Boostrapping...
                ( input,models,fitting_results,best_model,name1,name2,nbEmbryo_givenCondition,...
                tmin,general_param.cortex_analysis.fixed_short_lifetime,general_param.cortex_analysis.fixed_short_percent,save_stem,plot_visualization );
        end
    end
end


%% plot the fit of the best model with the raw data of each embryo

% no log scale
%to_plot_data_withBestFit_sum_mle( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,0 );

% x log scale
%to_plot_data_withBestFit_sum_mle( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,1 );

% y log scale
if error_computation == 1
    if ( choice_error_estimate == 1 && strcmp(best_model,'MonoExpo_stretched')== 0 && strcmp(best_model,'TripleExpo') == 0 ) 
        to_plot_data_withBestFit_sum_mle( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        
    elseif ( choice_error_estimate == 2 && strcmp(best_model,'MonoExpo_stretched')== 0 ) 
        to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        
    elseif ( choice_error_estimate == 3 && strcmp(best_model,'MonoExpo_stretched')== 0 ) 
        if strcmp(best_model,'TripleExpo') == 0
            to_plot_data_withBestFit_sum_mle( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        end
        to_plot_data_withBestFit_sum_mle_bootstrap( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        
        %to_plot_data_withBestFit_sum_mle2( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        %to_plot_data_withBestFit_sum_mle3( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
        
    else
        to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
    end
else
    to_plot_data_withBestFit_sum_mle_noError( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,2,classification_performed );
end


%% clear

clear input


end

