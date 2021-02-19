function to_plot_data_withAllFits_sum_mle...
    ( input,fitting_results,models,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,log_scale,classification_performed, given_set )

% enables to plot for each embryo the distribution (possibility to choose
% scale) with all the fitting models

if nargin < 10
    classification_performed = 0;
end
if nargin < 11
    given_set = '';
end

%% choice of scale for the representation

if log_scale == 0
    name0 = 'nolog';
elseif log_scale == 1
    name0 = 'xlog';
elseif log_scale ==2
    name0 = 'ylog';
end


%% plot for each embryo

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    name_embryo = ['embryo' num2str(iEmbryo)];
    size_population.(name_embryo).raw = nansum( input.(name_embryo).data(:,2) );
    binranges = input.(name_embryo).data(:,1);
    bincounts = input.(name_embryo).data(:,2);
    
    figure
    plot(binranges,bincounts,'ko');
    xlim([0 max(binranges)])
    hold all

    %-------------
    % model = mono-expo
    if ~isempty(find(ismember(models,'MonoExpo'),1))
        
        f_MonoExpo = @(x) size_population.(name_embryo).raw / fitting_results.MonoExpo.R_normalization.(name_embryo).mono ...
            * exp(-x/fitting_results.MonoExpo.T);
        fplot(f_MonoExpo,[0 max(binranges)],'-r');
        hold all
        
    end

    %-------------
    % model = double-expo
    if ~isempty(find(ismember(models,'DoubleExpo'),1))
        
        f_DoubleExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.DoubleExpo.P1 / fitting_results.DoubleExpo.R_normalization.(name_embryo).double(1,1) ...
            * exp(-x/fitting_results.DoubleExpo.T1)+ fitting_results.DoubleExpo.P2 / fitting_results.DoubleExpo.R_normalization.(name_embryo).double(2,1) * ...
            exp(-x/fitting_results.DoubleExpo.T2) );
        fplot(f_DoubleExpo,[0 max(binranges)], '-b');
        hold all
        
    end

        %-------------
    % model = double-expo with fixed T0
    if ~isempty(find(ismember(models,'DoubleExpo_fixedT0'),1))
        
        f_DoubleExpo_fixedT0 = @(x) size_population.(name_embryo).raw * ( fitting_results.DoubleExpo_fixedT0.P1 / fitting_results.DoubleExpo_fixedT0.R_normalization.(name_embryo).double_(1,1) ...
            * exp(-x/fitting_results.DoubleExpo_fixedT0.T1)+ fitting_results.DoubleExpo_fixedT0.P2 / fitting_results.DoubleExpo_fixedT0.R_normalization.(name_embryo).double_(2,1) * ...
            exp(-x/fitting_results.DoubleExpo_fixedT0.T2) );
        fplot(f_DoubleExpo_fixedT0,[0 max(binranges)], '--b');
        hold all
        
    end

            %-------------
    % model = double-expo with fixed T0P0
    if ~isempty(find(ismember(models,'DoubleExpo_fixedT0P0'),1))
        
        f_DoubleExpo_fixedT0P0 = @(x) size_population.(name_embryo).raw * ( fitting_results.DoubleExpo_fixedT0P0.P1 / fitting_results.DoubleExpo_fixedT0P0.R_normalization.(name_embryo).double__(1,1) ...
            * exp(-x/fitting_results.DoubleExpo_fixedT0P0.T1)+ fitting_results.DoubleExpo_fixedT0P0.P2 / fitting_results.DoubleExpo_fixedT0P0.R_normalization.(name_embryo).double__(2,1) * ...
            exp(-x/fitting_results.DoubleExpo_fixedT0P0.T2) );
        fplot(f_DoubleExpo_fixedT0P0,[0 max(binranges)], '-.');
        hold all
        
    end
    
    %-------------
    % model = mono-expo stretched    
    if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
        
        f_MonoExpo_stretched = @(x) size_population.(name_embryo).raw /  ...
            fitting_results.MonoExpo_stretched.R_normalization.(name_embryo).monoS * ...
            exp(- (x/ fitting_results.MonoExpo_stretched.Ts ) .^(1/fitting_results.MonoExpo_stretched.power) );
        fplot(f_MonoExpo_stretched,[0 max(binranges)],'-k');
        hold all
        
    end

    %-------------
    % model = triple-expo
    if ~isempty(find(ismember(models,'TripleExpo'),1))
        
        f_TripleExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.TripleExpo.PP1 / fitting_results.TripleExpo.R_normalization.(name_embryo).triple(1,1) ...
            * exp(-x/fitting_results.TripleExpo.TT1)+ fitting_results.TripleExpo.PP2 / fitting_results.TripleExpo.R_normalization.(name_embryo).triple(2,1) * ...
            exp(-x/fitting_results.TripleExpo.TT2)+ fitting_results.TripleExpo.PP3 / fitting_results.TripleExpo.R_normalization.(name_embryo).triple(3,1) * ...
            exp(-x/fitting_results.TripleExpo.TT3) );
        fplot(f_TripleExpo,[0 max(binranges)], '-c');
        hold all
        
    end
    
     %-------------
    % model = triple-expo with fixed T0
    if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1))
        
        f_TripleExpo_fixedT0 = @(x) size_population.(name_embryo).raw * ( fitting_results.TripleExpo_fixedT0.P0 / fitting_results.TripleExpo_fixedT0.R_normalization.(name_embryo).triple_(1,1) ...
            * exp(-x/fitting_results.TripleExpo_fixedT0.T0)+ fitting_results.TripleExpo_fixedT0.P1 / fitting_results.TripleExpo_fixedT0.R_normalization.(name_embryo).triple_(2,1) * ...
            exp(-x/fitting_results.TripleExpo_fixedT0.T1)+ fitting_results.TripleExpo_fixedT0.P2 / fitting_results.TripleExpo_fixedT0.R_normalization.(name_embryo).triple_(3,1) * ...
            exp(-x/fitting_results.TripleExpo_fixedT0.T2) );
        fplot(f_TripleExpo_fixedT0,[0 max(binranges)], '--c');
        hold all
        
    end

         %-------------
    % model = triple-expo with fixed T0 and P0
    if ~isempty(find(ismember(models,'TripleExpo_fixedT0P0'),1))
        
        f_TripleExpo_fixedT0P0 = @(x) size_population.(name_embryo).raw * ( fitting_results.TripleExpo_fixedT0P0.P0 / fitting_results.TripleExpo_fixedT0P0.R_normalization.(name_embryo).triple__(1,1) ...
            * exp(-x/fitting_results.TripleExpo_fixedT0P0.T0)+ fitting_results.TripleExpo_fixedT0P0.P1 / fitting_results.TripleExpo_fixedT0P0.R_normalization.(name_embryo).triple__(2,1) * ...
            exp(-x/fitting_results.TripleExpo_fixedT0P0.T1)+ fitting_results.TripleExpo_fixedT0P0.P2 / fitting_results.TripleExpo_fixedT0P0.R_normalization.(name_embryo).triple__(3,1) * ...
            exp(-x/fitting_results.TripleExpo_fixedT0P0.T2) );
        fplot(f_TripleExpo_fixedT0P0,[0 max(binranges)], '-.c');
        hold all
        
    end
    
    %-------------
    % model = quadro-expo
    if ~isempty(find(ismember(models,'QuadroExpo'),1))
        
        f_QuadroExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.QuadroExpo.PPP1 / fitting_results.QuadroExpo.R_normalization.(name_embryo).quadro(1,1) ...
            * exp(-x/fitting_results.QuadroExpo.TTT1) + fitting_results.QuadroExpo.PPP2 / fitting_results.QuadroExpo.R_normalization.(name_embryo).quadro(2,1) * ...
            exp(-x/fitting_results.QuadroExpo.TTT2) + fitting_results.QuadroExpo.PPP3 / fitting_results.QuadroExpo.R_normalization.(name_embryo).quadro(3,1) * ...
            exp(-x/fitting_results.QuadroExpo.TTT3) + fitting_results.QuadroExpo.PPP4 / fitting_results.QuadroExpo.R_normalization.(name_embryo).quadro(4,1) * ...
            exp(-x/fitting_results.QuadroExpo.TTT4) );
        fplot(f_QuadroExpo,[0 max(binranges)], '-g');
        hold all
        
    end
    
    %-------------
    % model = diffusion - drift
    if ~isempty(find(ismember(models,'Drift_diffusion'),1))
        
        f_DD = @(x) size_population.(name_embryo).raw / fitting_results.Drift_diffusion.R_normalization.(name_embryo).dd * x.^(-3/2) ...
            .* exp(-x/fitting_results.Drift_diffusion.T);
        fplot(f_DD,[0 max(binranges)],'-k');
        hold all
        
    end
  
    
    %-------------
    % cosmetic stuff    
    if log_scale == 1
        set(gca,'xscale','log');
    elseif log_scale ==2
        set(gca,'yscale','log');
    end
    
    legend({'raw data' models{:}})
    xlabel('tracks duration in sec')
    
    if classification_performed == 0
            title(['Duration distributions : ', name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+allfits-', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 1
            title(['Duration distributions : ', name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution_removalT0+allfits-', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 2
            title(['Duration distributions : ', given_set, ', ', name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+allfits_', given_set, '_', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    end
    saveas(gcf,[save_stem namePlot]);
    
    close(gcf)
    
end


end

