function to_plot_residuals_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models,best_model,name1,name2,save_stem,...
    classification_performed,given_set )

% get for each embryo a plot that shows the residues of the different
% models that are studied

% the best model is expected to have low residues distributed around zero

if nargin < 9
    classification_performed = 0;
end

if nargin < 10
    given_set = '';
end

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    name_embryo = ['embryo' num2str(iEmbryo)];
    binranges = input.(name_embryo).data(:,1);
    
    figure
    cst_raw = zeros(length(binranges),1);
    plot(cst_raw,':k');
    hold all
    for j = 1 : numel(models)
        name_model = models{j};
        plot(transpose(fitting_results.(name_model).(name_embryo).residuals(:,1)),transpose(fitting_results.(name_model).(name_embryo).residuals(:,2)));
        hold all
    end
    legend('cst=0',{'data' models{:}})
    xlim([0 max(binranges)])
    ylim([-2 2])
    ylabel('Residues')
    xlabel('tracks duration in sec')
    string = ['flag (' best_model ')= ' num2str(fitting_results.(best_model).flag)];
    text(50,50,string,'Units','pixels')
    
    if classification_performed == 0
            title(['Residuals of the tracks duration distributions fitting: ', name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_residuals-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 1
            title(['Residuals of the tracks duration distributions fitting: ', name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_residuals_T0removal-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 2
            title(['Residuals of the duration distribution fit for: ', given_set, ', ',name1, ' and ', 10, name2, '( best model: ', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_residuals_', given_set, '_', name_embryo, '-', name1 , '-', name2, '.fig');
    end
    saveas(gcf,[save_stem namePlot]);
    close(gcf)
    
    clear bincounts
    
end


end

