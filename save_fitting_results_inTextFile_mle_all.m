function save_fitting_results_inTextFile_mle_all( final_fitting_all,models,save_stem,analysis_3regions,nbEmbryo_metaphase,nbEmbryo,...
    use_loglikelihood_to_calculate_parameters,minLength_tracks,nbEmbryo_late,conditions1,conditions2,folder_tag,choice_polarity,...
    polarity_lost_choice,polarity_metaphase_50,algo_choice,error_computation,choice_error_estimate,nbEmbryo_prophase,nbEmbryo_prometaphase,...
    nbEmbryo_anaphase,classification_analysis,threshold_value_max_speed,last_bin_kept,time_reference_choice)

global general_param

if nargin < 12
    folder_tag = 'no tag';
end

if nargin < 17
    error_computation = 0;
end

if nargin < 18
    choice_error_estimate = NaN;
end

if nargin < 19
    nbEmbryo_prophase = 0;
end

if nargin < 20
    nbEmbryo_prometaphase = 0;
end
if nargin < 21
    nbEmbryo_anaphase = nbEmbryo;
end

if nargin < 22
    classification_analysis = 0;
end
if nargin < 24
    last_bin_kept = 1;
end
if nargin < 25
    time_reference_choice = 0;
end

name_fid = [save_stem 'tracks_duration_histo_values-BayesianInference.txt'];
fid = fopen(name_fid ,'wt');

fprintf(fid,'total number of embryos: \n');
fprintf(fid,'%f\n',nbEmbryo );
if any(strcmp(conditions1,'prophase'))
    fprintf(fid,'prophase number of embryos: \n');
    fprintf(fid,'%f\n',nbEmbryo_prophase );
end
if any(strcmp(conditions1,'metaphase'))
    fprintf(fid,'metaphase number of embryos: \n');
    fprintf(fid,'%f\n',nbEmbryo_metaphase );
end
if any(strcmp(conditions1,'prometaphase'))
    fprintf(fid,'prometaphase number of embryos: \n');
    fprintf(fid,'%f\n',nbEmbryo_prometaphase );
end
if any(strcmp(conditions1,'anaphase'))
    fprintf(fid,'anaphase number of embryos: \n');
    fprintf(fid,'%f\n',nbEmbryo_anaphase );
end
if any(strcmp(conditions1,'late'))
    fprintf(fid,'late number of embryos: \n');
    fprintf(fid,'%f\n',nbEmbryo_late );
end
fprintf(fid,'total temporal reference: \n');
if time_reference_choice == 0
    fprintf(fid,'%s\n', 'cytokinesis furrowing' );
elseif time_reference_choice == 1
    fprintf(fid,'%s\n', 'pseudo cleavage' );
end
fprintf(fid,'minimum length kept for fitting: \n');
fprintf(fid,'%f\n',minLength_tracks );
fprintf(fid,'use of threshold count in duration distribution: \n');
fprintf(fid,'%f\n',general_param.cortex_analysis.use_threshold_count_fit);
fprintf(fid,'if used, threshold in count kept: \n');
fprintf(fid,'%f\n',general_param.cortex_analysis.binCount_threshold );
fprintf(fid,'use of smoothing duration distribution: \n');
fprintf(fid,'%f\n',general_param.cortex_analysis.use_smoothing_fit );
fprintf(fid,'if used, size of smoothing window: \n');
fprintf(fid,'%f\n',general_param.cortex_analysis.binCount_smoothW_size );
fprintf(fid,'tag of the folder: \n');
fprintf(fid,'%s\n',folder_tag );
fprintf(fid,'polarity selection: \n');
fprintf(fid,'%f\n',choice_polarity );
if choice_polarity == 1
    fprintf(fid,'polarity lost: \n');
    fprintf(fid,'%f\n',polarity_lost_choice );
end
fprintf(fid,'algo chosen (1 = fmincon, 2 = patternsearch): \n');
fprintf(fid,'%f\n',algo_choice );
if analysis_3regions == 0
    fprintf(fid,'polarity metaphase set to 50: \n');
    fprintf(fid,'%f\n',polarity_metaphase_50 );
elseif analysis_3regions == 1
    fprintf(fid,'boundary #1 of the 3 regions: \n');
    fprintf(fid,'%f\n',general_param.cortex_analysis.analysis_3SpatialRegions_limit1 );
    fprintf(fid,'boundary #2 of the 3 regions: \n');
    fprintf(fid,'%f\n',general_param.cortex_analysis.analysis_3SpatialRegions_limit2 );
end

if classification_analysis == 1
        fprintf(fid,'classification according to max speed performed, threshold set to: \n');
    fprintf(fid,'%f\n',threshold_value_max_speed );
end
fprintf(fid,'last bin considered in the fitting: \n');
fprintf(fid,'%f\n',last_bin_kept );

% if nargin < 10
%     if nbEmbryo_metaphase < 2 && nbEmbryo_late < 2
%         conditions1 = {'entireRecording','anaphase'};
%     elseif nbEmbryo_metaphase < 2 && nbEmbryo_late >= 2
%         conditions1 = {'entireRecording','anaphase','late'};
%     elseif nbEmbryo_metaphase >= 2 && nbEmbryo_late < 2
%         conditions1 = {'entireRecording','metaphase','anaphase'};
%     else
%         conditions1 = {'entireRecording','metaphase','anaphase','late'};
%     end
% end

if nargin < 11
    if analysis_3regions == 0
        %conditions2 = {'entireEmbryo','anterior','posterior','anterior_reduced'};
        conditions2 = {'entireEmbryo','anterior','posterior'};
    elseif analysis_3regions == 1
        conditions2 = {'region1','region2','region3'};
    end
end

n_condi1 = numel(conditions1);
n_condi2 = numel(conditions2);
n_condi = n_condi1 * n_condi2;
conditions = cell(3,n_condi);
count_condition = 0;

for iCondition=1:n_condi1
    for jCondition=1:n_condi2
        count_condition = count_condition +1;
        conditions{1,count_condition} = strcat(conditions1{iCondition}, '_', conditions2{jCondition});
        conditions{2,count_condition} = conditions1{iCondition};
        conditions{3,count_condition} = conditions2{jCondition};
    end
end

fprintf(fid,'conditions \t');
for i = 1 : n_condi
    fprintf(fid,'%s ', conditions{1,i} );
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'\t\t\t\t');
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');    
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'model \t');
for i = 1 : n_condi
    for j = 1 : numel(models)
        name_model = models{j};
        fprintf(fid,'%s ', models{1,j} );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'\t\t\t\t');
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');              
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'probability AIC \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        fprintf(fid,'%f',round2( final_fitting_all.(name1).(name2).model.(name_model).PrM_AIC,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'\t\t\t\t');
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');              
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'probability AICc \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        fprintf(fid,'%f',round2( final_fitting_all.(name1).(name2).model.(name_model).PrM_AICc,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'\t\t\t\t');
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');              
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'probability BIC \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        fprintf(fid,'%f',round2( final_fitting_all.(name1).(name2).model.(name_model).PrM_BIC,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'\t\t\t\t');
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');              
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'number of contacts \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total) );
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total3) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total0) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total0) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total1) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total2) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total3) );
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total4) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'%f\t',round( final_fitting_all.(name1).(name2).fitting.(name_model).size_population.total) );           
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'residency time \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T,1e-2) );
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2,1e-2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT3,1e-2) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).Ts,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2,1e-2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TTT1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TTT2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TTT3,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TTT4,1e-2) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T,1e-2) );             
        end
    end
end
fprintf(fid,'\n');

if use_loglikelihood_to_calculate_parameters == 1 
    fprintf(fid,'residency time (likelihood) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T_lik,1e-2) );
            elseif strcmp(name_model,'DoubleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1_lik,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2_lik,1e-2) );
            elseif strcmp(name_model,'DoubleExpo_fixedT0')
                fprintf(fid,'\t\t');
            elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
                fprintf(fid,'\t\t');
            elseif strcmp(name_model,'TripleExpo')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0P0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'QuadroExpo')
                fprintf(fid,'\t\t\t\t');
            elseif strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )
    fprintf(fid,'residency time (bootstrapping) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T_mean_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t');
                end
            elseif strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T1_mean_bootstrap,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).T2_mean_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            elseif strcmp(name_model,'DoubleExpo_fixedT0')
                fprintf(fid,'\t\t');
            elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
                fprintf(fid,'\t\t');
            elseif strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT1_mean_bootstrap,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT2_mean_bootstrap,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).TT3_mean_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0P0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'QuadroExpo')
                fprintf(fid,'\t\t\t\t');
            elseif strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 1 || choice_error_estimate == 3 )
    fprintf(fid,'residency time se \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T_se,1e-2) );
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                 end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                end
            end
            if strcmp(name_model,'TripleExpo')
%                 if strcmp(best_model,'TripleExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT2_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT3_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                    %fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).Ts_se,1e-2) );
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
%                 if strcmp(best_model,'TripleExpo_fixedT0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                 end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
%                 if strcmp(best_model,'TripleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                 end
            end
            if strcmp(name_model,'QuadroExpo')
%                 if strcmp(best_model,'QuadroExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT2_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT3_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT4_se,1e-2) );
%                 else
                     fprintf(fid,'\t\t\t\t');
%                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )
    fprintf(fid,'residency time se (bootstrap) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T_se_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se_bootstrap,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                 end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT1_se_bootstrap,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT2_se_bootstrap,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TT3_se_bootstrap,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                    %fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).Ts_se,1e-2) );
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
%                 if strcmp(best_model,'TripleExpo_fixedT0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                 end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
%                 if strcmp(best_model,'TripleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).T2_se,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                 end
            end
            if strcmp(name_model,'QuadroExpo')
%                 if strcmp(best_model,'QuadroExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT1_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT2_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT3_se,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).TTT4_se,1e-2) );
%                 else
                     fprintf(fid,'\t\t\t\t');
%                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end            
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'percentage \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2*100,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2*100,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2*100,1e-2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP2*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP3*100,1e-2) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%s%f\t','power: ',round2( final_fitting_all.(name1).(name2).fitting.(name_model).power,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P0*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2*100,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P0*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2*100,1e-2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PPP1*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PPP2*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PPP3*100,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PPP4*100,1e-2) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'\t');
        end
    end
end
fprintf(fid,'\n');

if use_loglikelihood_to_calculate_parameters == 1
    fprintf(fid,'percentage (likelihood) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'DoubleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1_lik*100,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2_lik*100,1e-2) );
            elseif strcmp(name_model,'DoubleExpo_fixedT0')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0P0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'QuadroExpo')
                fprintf(fid,'\t\t\t\t');
            elseif strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end


if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )
    fprintf(fid,'percentage (bootstrapping) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P1_mean_bootstrap*100,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).P2_mean_bootstrap*100,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            elseif strcmp(name_model,'DoubleExpo_fixedT0')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
                fprintf(fid,'\t');
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP1_mean_bootstrap*100,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP2_mean_bootstrap*100,1e-2) );
                fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).PP3_mean_bootstrap*100,1e-2) );    
                else
                    fprintf(fid,'\t\t\t');
                end
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'TripleExpo_fixedT0P0')
                fprintf(fid,'\t\t\t');
            elseif strcmp(name_model,'QuadroExpo')
                fprintf(fid,'\t\t\t\t');
            elseif strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end


if error_computation == 1 && ( choice_error_estimate == 1 || choice_error_estimate == 3 )
    fprintf(fid,'percentage se \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
                if strcmp(best_model,'DoubleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
%                 if strcmp(best_model,'TripleExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
%                 if strcmp(best_model,'TripleExpo_fixedT0')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
%                 if strcmp(best_model,'TripleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'QuadroExpo')
%                 if strcmp(best_model,'QuadroExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t\t');
%                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )
    fprintf(fid,'percentage se (bootstrap) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se_bootstrap*100,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se_bootstrap*100,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
%                 if strcmp(best_model,'DoubleExpo_fixedT0')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t');
%                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
                if strcmp(best_model,'DoubleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se_bootstrap*100,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se_bootstrap*100,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se_bootstrap*100,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
%                 if strcmp(best_model,'TripleExpo_fixedT0')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
%                 if strcmp(best_model,'TripleExpo_fixedT0P0')
%                     fprintf(fid,'\t');
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t');
%                end
            end
            if strcmp(name_model,'QuadroExpo')
%                 if strcmp(best_model,'QuadroExpo')
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se*100,1e-2) );
%                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se*100,1e-2) );
%                 else
                    fprintf(fid,'\t\t\t\t');
%                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'nb MTs/s \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean,1e-2) );
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean3,1e-2) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean3,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean4,1e-2) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean,1e-2) );            
        end
    end
end
fprintf(fid,'\n');

if error_computation == 1 && ( choice_error_estimate == 1 || choice_error_estimate == 3 )
    fprintf(fid,'nb MTs/s se \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
                if strcmp(best_model,'TripleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
                if strcmp(best_model,'TripleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'QuadroExpo')
                if strcmp(best_model,'QuadroExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t\t');
                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )  % notice, before split into 2 the erro, when both error were generated, only the bootstrap one was seave in text file
    fprintf(fid,'nb MTs/s se (bootstrap) \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );                        
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se_bootstrap  * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se_bootstrap  * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );                        
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
                if strcmp(best_model,'TripleExpo_fixedT0')
                    fprintf(fid,'\t\t\t');
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
                if strcmp(best_model,'TripleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'QuadroExpo')
                if strcmp(best_model,'QuadroExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t\t');
                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'nb MTs/min \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean*60,1e-2) );
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean3*60,1e-2) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean*60,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean0*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean0*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean1*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean2*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean3*60,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean4*60,1e-2) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_time.mean*60,1e-2) );            
        end
    end
end
fprintf(fid,'\n');

if error_computation == 1
    fprintf(fid,'nb MTs/min se \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    if ( choice_error_estimate == 2 || choice_error_estimate == 3 )
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );                        
                    else
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    end
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
                if strcmp(best_model,'DoubleExpo_fixedT0P0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    if ( choice_error_estimate == 2 || choice_error_estimate == 3 )
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );                        
                    else
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    end
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
                if strcmp(best_model,'TripleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
                if strcmp(best_model,'TripleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'QuadroExpo')
                if strcmp(best_model,'QuadroExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_time.mean*60,1e-2) );
                else
                    fprintf(fid,'\t\t\t\t');
                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'nb MTs/min/um2 \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean,1e-2) );
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
        elseif strcmp(name_model,'DoubleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean3,1e-2) );
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
        elseif strcmp(name_model,'TripleExpo_fixedT0P0')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean0,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
        elseif strcmp(name_model,'QuadroExpo')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean1,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean2,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean3,1e-2) );
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean4,1e-2) );
        elseif strcmp(name_model,'Drift_diffusion')
            fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(name_model).size_population_normalized_timeAndArea.mean,1e-2) );            
        end
    end
end
fprintf(fid,'\n');

if error_computation == 1 && ( choice_error_estimate == 1 || choice_error_estimate == 3 )
    fprintf(fid,'nb MTs/min/um2 se \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
                if strcmp(best_model,'DoubleExpo_fixedT0P0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
                if strcmp(best_model,'TripleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
                if strcmp(best_model,'TripleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'QuadroExpo')
                if strcmp(best_model,'QuadroExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t\t');
                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

if error_computation == 1 && ( choice_error_estimate == 2 || choice_error_estimate == 3 )
    fprintf(fid,'nb MTs/min/um2 se (boostrap) \t'); % notice, before split into 2 the erro, when both error were generated, only the bootstrap one was seave in text file
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        best_model = final_fitting_all.(name1).(name2).model.best_model_BIC;
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                if strcmp(best_model,'MonoExpo')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'DoubleExpo')
                if strcmp(best_model,'DoubleExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );                        
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0')
                if strcmp(best_model,'DoubleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'DoubleExpo_fixedT0P0')
                if strcmp(best_model,'DoubleExpo_fixedT0P0')
                    fprintf(fid,'\t\t');
                else
                    fprintf(fid,'\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo')
                if strcmp(best_model,'TripleExpo')
                     fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP1_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP2_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PP3_se_bootstrap * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );                       
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'MonoExpo_stretched')
                if strcmp(best_model,'MonoExpo_stretched')
                    fprintf(fid,'\t');
                else
                    fprintf(fid,'\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0')
                if strcmp(best_model,'TripleExpo_fixedT0')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P0_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'TripleExpo_fixedT0P0')
                if strcmp(best_model,'TripleExpo_fixedT0P0')
                    fprintf(fid,'\t');
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).P2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t');
                end
            end
            if strcmp(name_model,'QuadroExpo')
                if strcmp(best_model,'QuadroExpo')
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP1_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP2_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP3_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                    fprintf(fid,'%f\t',round2( final_fitting_all.(name1).(name2).fitting.(best_model).PPP4_se * ...
                        final_fitting_all.(name1).(name2).fitting.(best_model).size_population_normalized_timeAndArea.mean,1e-2) );
                else
                    fprintf(fid,'\t\t\t\t');
                end
            end
            if strcmp(name_model,'Drift_diffusion')
                fprintf(fid,'\t');
            end
        end
    end
    fprintf(fid,'\n');
end

fclose(fid);

end

