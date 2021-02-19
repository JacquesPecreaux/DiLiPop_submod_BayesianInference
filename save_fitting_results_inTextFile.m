function save_fitting_results_inTextFile( tracks_duration_histo,save_stem,name_embryo,models,iEmbryo,several_embryos,exist_metaphase )

name_fid{iEmbryo} = [save_stem 'tracks_duration_histo_values-BayesianInference-', name_embryo,'.txt'];
fid = fopen(name_fid{iEmbryo} ,'wt');

if exist_metaphase ==1
    conditions1 = {'metaphase','anaphase'};
elseif exist_metaphase == 0
    conditions1 = {'anaphase'};
end
conditions2 = {'entireEmbryo','anterior','posterior'};

n_condi1 = numel(conditions1);
n_condi2 = numel(conditions2);
n_condi = n_condi1 * n_condi2;
conditions = cell(3,n_condi);
count_condition = 0;

for iCondition=1:n_condi1 ;
    for jCondition=1:n_condi2 ;
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
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
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
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
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
        fprintf(fid,'%f',round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.(name_model).PrM_AIC,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
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
        fprintf(fid,'%f',round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.(name_model).PrM_AICc,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
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
        fprintf(fid,'%f',round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.(name_model).PrM_BIC,1e-2) );
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'\t');
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'prefactor \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.A,1e-2));
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.B1,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.B2,1e-2));
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C1,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C2,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C3,1e-2));
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.D,1e-2));
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'se prefactor \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.A_se,1e-2));
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.B1_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.B2_se,1e-2));
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C1_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C2_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.C3_se,1e-2));
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.D_se,1e-2));
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'tau \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.T,1e-2));
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.T1,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.T2,1e-2));
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT1,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT2,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT3,1e-2));
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.T,1e-2));
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'se tau \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.T_se,1e-2));
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.T1_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.T2_se,1e-2));
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT1_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT2_se,1e-2));
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.TT3_se,1e-2));
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.T_se,1e-2));
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'size population \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.size_population));
        elseif strcmp(name_model,'DoubleExpo')
            
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2) );
            
        elseif strcmp(name_model,'TripleExpo')
            
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population2) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population3) );
            
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.size_population));
        end
    end
end
fprintf(fid,'\n');


fprintf(fid,'size population se \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.size_population_se));
            
        elseif strcmp(name_model,'DoubleExpo')
            
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1_se) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2_se) );
            
        elseif strcmp(name_model,'TripleExpo')
            
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1_se) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population2_se) );
            fprintf(fid,'%f\t', round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population3_se) );
            
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round(tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.size_population_se));
        end
    end
end
fprintf(fid,'\n');


fprintf(fid,'percentage \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.size_population / ...
                tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.size_population,1e-2) ); %??
            
        elseif strcmp(name_model,'DoubleExpo')
            
            fprintf(fid,'%f\t', round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population,1e-2) );
            fprintf(fid,'%f\t', round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population,1e-2) );
            
        elseif strcmp(name_model,'TripleExpo')
            
            fprintf(fid,'%f\t', round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population,1e-2) );
            fprintf(fid,'%f\t', round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population2 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population,1e-2) );
            fprintf(fid,'%f\t', round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population3 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population,1e-2) );
            
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( 100*tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.size_population / ...
                tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.size_population,1e-2) ); %??
        end
    end
end
fprintf(fid,'\n');


fprintf(fid,'se percentage \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo.percent_population_se,1e-2) );
            
        elseif strcmp(name_model,'DoubleExpo')
            
            fprintf(fid,'%f\t', round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.percent_population1_se,1e-2) );
            fprintf(fid,'%f\t', round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.percent_population2_se,1e-2) );
            
        elseif strcmp(name_model,'TripleExpo')
            
            fprintf(fid,'%f\t', round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.percent_population1_se,1e-2) );
            fprintf(fid,'%f\t', round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.percent_population2_se,1e-2) );
            fprintf(fid,'%f\t', round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.percent_population3_se,1e-2) );
            
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round2( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.MonoExpo_stretched.percent_population_se,1e-2) );
        end
    end
end
fprintf(fid,'\n');

fprintf(fid,'number of MTs \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) ); %??
            
        elseif strcmp(name_model,'DoubleExpo')
            
            fprintf(fid,'%f\t', round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population * tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
            fprintf(fid,'%f\t', round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population * tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
            
        elseif strcmp(name_model,'TripleExpo')
            
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population2 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population3 ...
                /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
            
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks) );
        end
    end
end
fprintf(fid,'\n');



if several_embryos == 1
    
    fprintf(fid,'number of MTs per embryo \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        for j = 1 : numel(models)
            name_model = models{j};
            if strcmp(name_model,'MonoExpo')
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks ...
                    / tracks_duration_histo{iEmbryo}.number_embryo) );
            elseif strcmp(name_model,'DoubleExpo')
                
                fprintf(fid,'%f\t', round ( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo) );
                fprintf(fid,'%f\t', round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo) );
                
            elseif strcmp(name_model,'TripleExpo')
                
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo) );
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population2 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo) );
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population3 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo) );
                
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks ...
                    / tracks_duration_histo{iEmbryo}.number_embryo) );
            end
        end
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'number of MTs per embryo per sec \t');
    for i = 1 : n_condi
        name1 = conditions{2,i};
        name2 = conditions{3,i};
        for j = 1 : numel(models)
            name_model = models{j};
            
            if strcmp(name_model,'MonoExpo')
                
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks ...
                    / tracks_duration_histo{iEmbryo}.number_embryo * tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                
            elseif strcmp(name_model,'DoubleExpo')
                
                fprintf(fid,'%f\t', round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                fprintf(fid,'%f\t', round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population2 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.DoubleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                
            elseif strcmp(name_model,'TripleExpo')
                
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population1 ...
                    /tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.TripleExpo.size_population* tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks / ...
                    tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
                
            elseif strcmp(name_model,'MonoExpo_stretched')
                fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).number_tracks ...
                    / tracks_duration_histo{iEmbryo}.number_embryo* tracks_duration_histo{iEmbryo}.(name1).(name2).timeDuration_phase_mean) );
            end
            
        end
    end
    fprintf(fid,'\n');
    
end




fprintf(fid,'power \t');
for i = 1 : n_condi
    name1 = conditions{2,i};
    name2 = conditions{3,i};
    for j = 1 : numel(models)
        name_model = models{j};
        if strcmp(name_model,'MonoExpo')
            fprintf(fid,'\t');
        elseif strcmp(name_model,'DoubleExpo')
            fprintf(fid,'\t\t');
        elseif strcmp(name_model,'TripleExpo')
            fprintf(fid,'\t\t\t');
        elseif strcmp(name_model,'MonoExpo_stretched')
            fprintf(fid,'%f\t',round( tracks_duration_histo{iEmbryo}.(name1).(name2).fitting.(name_model).power) );
        end
    end
end
fprintf(fid,'\n');

fclose(fid);


end

