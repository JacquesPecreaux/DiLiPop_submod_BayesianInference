function to_plot_data_withBestFit_sum_mle...
    ( input,fitting_results,best_model,nbEmbryo_givenCondition,name1,name2,save_stem,log_scale,classification_performed,given_set )

% enables to plot for each embryo the distribution (possibility to choose
% scale) with the best fitting model: values of the parameters are written
% on the plots

if nargin < 9
    classification_performed = 0;
end
if nargin < 10
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
    % best model = mono-expo
    
    if strcmp(best_model,'MonoExpo')
        
        f_MonoExpo = @(x) size_population.(name_embryo).raw / fitting_results.(best_model).R_normalization.(name_embryo).mono * exp(-x/fitting_results.(best_model).T);
        fplot(f_MonoExpo,[0 max(binranges)],'-r')
        
        string3 = ['residency time = ' num2str(round2(fitting_results.(best_model).T,1e-2)) '+/-' ...
            num2str(round2(fitting_results.(best_model).T_se,1e-2)) ' sec.'];
        text(25,50,string3,'Units','pixels')
        
        string1 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string1,'Units','pixels')
        
    %-------------
    % best model = double-expo        
    elseif strcmp(best_model,'DoubleExpo')
        
        f_DoubleExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).P1 / fitting_results.(best_model).R_normalization.(name_embryo).double(1,1) ...
            * exp(-x/fitting_results.(best_model).T1) +  fitting_results.(best_model).P2 / fitting_results.(best_model).R_normalization.(name_embryo).double(2,1) ...
            * exp(-x/fitting_results.(best_model).T2) );
        
        fplot(f_DoubleExpo,[0 max(binranges)], '-r');
        
        string4 = ['residency time = ' num2str(round2(fitting_results.(best_model).T1,1e-2)) '+/-' ...
            num2str(round2(fitting_results.(best_model).T1_se,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P1*100,1e-1)) ' % +/-' ...
            num2str(round2(fitting_results.(best_model).P1_se*100,1e-1)) ];
        text(25,50,string4,'Units','pixels')
        
        string6 = ['residency time = ' num2str(round2(fitting_results.(best_model).T2,1e-2))  '+/-' ...
            num2str(round2(fitting_results.(best_model).T2_se,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P2*100,1e-1)) ' % +/-' ...
            num2str(round2(fitting_results.(best_model).P2_se*100,1e-1)) ];
        text(25,25,string6,'Units','pixels')
        
        string2 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string2,'Units','pixels')

    %-------------
    % best model = double-expo with fixed T0        
    elseif strcmp(best_model,'DoubleExpo_fixedT0')
        
        f_DoubleExpo_fixedT0 = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).P1 / fitting_results.(best_model).R_normalization.(name_embryo).double_(1,1) ...
            * exp(-x/fitting_results.(best_model).T1) +  fitting_results.(best_model).P2 / fitting_results.(best_model).R_normalization.(name_embryo).double_(2,1) ...
            * exp(-x/fitting_results.(best_model).T2) );
        
        fplot(f_DoubleExpo_fixedT0,[0 max(binranges)], '-r');
        
        string4 = ['residency time = ' num2str(round2(fitting_results.(best_model).T1,1e-2)) ' s. and ' ...
            num2str(round2(fitting_results.(best_model).P1*100,1e-1)) ' +/- ' ...
            num2str(round2(fitting_results.(best_model).P1_se*100,1e-1)) '%'];
        text(25,50,string4,'Units','pixels')
        
        string6 = ['residency time = ' num2str(round2(fitting_results.(best_model).T2,1e-2))  '+/-' ...
            num2str(round2(fitting_results.(best_model).T2_se,1e-2)) ' s and ' ...
            num2str(round2(fitting_results.(best_model).P2*100,1e-1)) ' +/-' ...
            num2str(round2(fitting_results.(best_model).P2_se*100,1e-1)) '%' ];
        text(25,25,string6,'Units','pixels')
        
        string2 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string2,'Units','pixels')

            %-------------
    % best model = double-expo with fixed T0P0        
    elseif strcmp(best_model,'DoubleExpo_fixedT0P0')
        
        f_DoubleExpo_fixedT0P0 = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).P1 / fitting_results.(best_model).R_normalization.(name_embryo).double__(1,1) ...
            * exp(-x/fitting_results.(best_model).T1) +  fitting_results.(best_model).P2 / fitting_results.(best_model).R_normalization.(name_embryo).double__(2,1) ...
            * exp(-x/fitting_results.(best_model).T2) );
        
        fplot(f_DoubleExpo_fixedT0P0,[0 max(binranges)], '-r');
        
        string4 = ['residency time = ' num2str(round2(fitting_results.(best_model).T1,1e-2)) ' s. and ' ...
            num2str(round2(fitting_results.(best_model).P1*100,1e-1)) '%'];
        text(25,50,string4,'Units','pixels')
        
        string6 = ['residency time = ' num2str(round2(fitting_results.(best_model).T2,1e-2))  '+/-' ...
            num2str(round2(fitting_results.(best_model).T2_se,1e-2)) ' s and ' ...
            num2str(round2(fitting_results.(best_model).P2*100,1e-1)) '%' ];
        text(25,25,string6,'Units','pixels')
        
        string2 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string2,'Units','pixels')
        
    %-------------
    % best model = simple-expo stretched        
    elseif strcmp(best_model,'MonoExpo_stretched')
        
        f_MonoExpo_stretched = @(x) size_population.(name_embryo).raw / fitting_results.(best_model).R_normalization.(name_embryo).monoS...
            * exp(- (x/ fitting_results.(best_model).Ts ) .^(1/fitting_results.(best_model).power) );
        fplot(f_MonoExpo_stretched,[0 max(binranges)],'-r')
        
        string3 = ['residency time = ' num2str(round2(fitting_results.(best_model).Ts,1e-2)) ' sec. and power = '...
            num2str(round2(fitting_results.(best_model).power*100,1e-2))];
        text(25,50,string3,'Units','pixels')
        
        string1 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string1,'Units','pixels')

    %-------------
    % best model = triple-expo        
    elseif strcmp(best_model,'TripleExpo')
        
        f_TripleExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).PP1 / fitting_results.(best_model).R_normalization.(name_embryo).triple(1,1) ...
            * exp(-x/fitting_results.(best_model).TT1)+ fitting_results.(best_model).PP2 / fitting_results.(best_model).R_normalization.(name_embryo).triple(2,1) ...
            * exp(-x/fitting_results.(best_model).TT2)+ fitting_results.(best_model).PP3 / fitting_results.(best_model).R_normalization.(name_embryo).triple(3,1) ...
            * exp(-x/fitting_results.(best_model).TT3) ) ;
        
        fplot(f_TripleExpo,[0 max(binranges)], '-r');
        
        string7 = ['residency time = ' num2str(round2(fitting_results.(best_model).TT1,1e-2)) '+/-' ...
            num2str(round2(fitting_results.(best_model).TT1_se,1e-2)) ' s and ' ...
            num2str(round2(fitting_results.(best_model).PP1*100,1e-1)) '+/-' ...
            num2str(round2(fitting_results.(best_model).PP1_se*100,1e-1)) '%'];
        text(25,75,string7,'Units','pixels')
        
        string8 = ['residency time = ' num2str(round2(fitting_results.(best_model).TT2,1e-2)) '+/-' ...
            num2str(round2(fitting_results.(best_model).TT2_se,1e-2)) ' s and ' ...
            num2str(round2(fitting_results.(best_model).PP2*100,1e-1)) '+/-' ...
            num2str(round2(fitting_results.(best_model).PP2_se*100,1e-1)) '%'];
        text(25,50,string8,'Units','pixels')
        
        string9 = ['residency time = ' num2str(round2(fitting_results.(best_model).TT3,1e-2)) '+/-' ...
            num2str(round2(fitting_results.(best_model).TT3_se,1e-2)) ' s and ' ...
            num2str(round2(fitting_results.(best_model).PP3*100,1e-1)) '+/-' ...
            num2str(round2(fitting_results.(best_model).PP3_se*100,1e-1)) '%'];
        text(25,25,string9,'Units','pixels')
        
        string10 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string10,'Units','pixels')
        
    %-------------
    % best model = triple-expo with fixed T0        
    elseif strcmp(best_model,'TripleExpo_fixedT0')
        
        f_TripleExpo_fixedT0 = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).P0 / fitting_results.(best_model).R_normalization.(name_embryo).triple_(1,1) ...
            * exp(-x/fitting_results.(best_model).T0)+ fitting_results.(best_model).P1 / fitting_results.(best_model).R_normalization.(name_embryo).triple_(2,1) ...
            * exp(-x/fitting_results.(best_model).T1)+ fitting_results.(best_model).P2 / fitting_results.(best_model).R_normalization.(name_embryo).triple_(3,1) ...
            * exp(-x/fitting_results.(best_model).T2) ) ;
        
        fplot(f_TripleExpo_fixedT0,[0 max(binranges)], '-r');
        
        string7 = ['residency time = ' num2str(round2(fitting_results.(best_model).T0,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P0*100,1e-1)) ' % '];
        text(25,75,string7,'Units','pixels')
        
        string8 = ['residency time = ' num2str(round2(fitting_results.(best_model).T1,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P1*100,1e-1)) ' % '];
        text(25,50,string8,'Units','pixels')
        
        string9 = ['residency time = ' num2str(round2(fitting_results.(best_model).T2,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P2*100,1e-1)) ' % '];
        text(25,25,string9,'Units','pixels')
        
        string10 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string10,'Units','pixels')        
        
    %-------------
    % best model = triple-expo with fixed T0 and P0        
    elseif strcmp(best_model,'TripleExpo_fixedT0P0')
        
        f_TripleExpo_fixedT0P0 = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).P0 / fitting_results.(best_model).R_normalization.(name_embryo).triple__(1,1) ...
            * exp(-x/fitting_results.(best_model).T0)+ fitting_results.(best_model).P1 / fitting_results.(best_model).R_normalization.(name_embryo).triple__(2,1) ...
            * exp(-x/fitting_results.(best_model).T1)+ fitting_results.(best_model).P2 / fitting_results.(best_model).R_normalization.(name_embryo).triple__(3,1) ...
            * exp(-x/fitting_results.(best_model).T2) ) ;
        
        fplot(f_TripleExpo_fixedT0P0,[0 max(binranges)], '-r');
        
        string7 = ['residency time = ' num2str(round2(fitting_results.(best_model).T0,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P0*100,1e-1)) ' % '];
        text(25,75,string7,'Units','pixels')
        
        string8 = ['residency time = ' num2str(round2(fitting_results.(best_model).T1,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P1*100,1e-1)) ' % '];
        text(25,50,string8,'Units','pixels')
        
        string9 = ['residency time = ' num2str(round2(fitting_results.(best_model).T2,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).P2*100,1e-1)) ' % '];
        text(25,5,string9,'Units','pixels')
        
        string10 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string10,'Units','pixels')        

            %-------------
    % best model =quadro-expo        
    elseif strcmp(best_model,'QuadroExpo')
        
        f_QuadroExpo = @(x) size_population.(name_embryo).raw * ( fitting_results.(best_model).PPP1 / fitting_results.(best_model).R_normalization.(name_embryo).quadro(1,1) ...
            * exp(-x/fitting_results.(best_model).TTT1)+ fitting_results.(best_model).PPP2 / fitting_results.(best_model).R_normalization.(name_embryo).quadro(2,1) ...
            * exp(-x/fitting_results.(best_model).TTT2)+ fitting_results.(best_model).PPP3 / fitting_results.(best_model).R_normalization.(name_embryo).quadro(3,1) ...
            * exp(-x/fitting_results.(best_model).TTT3)+ fitting_results.(best_model).PPP4 / fitting_results.(best_model).R_normalization.(name_embryo).quadro(4,1) ...
            * exp(-x/fitting_results.(best_model).TTT4) ) ;
        
        fplot(f_QuadroExpo,[0 max(binranges)], '-r');
        
        string7 = ['residency time = ' num2str(round2(fitting_results.(best_model).TTT1,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).PPP1*100,1e-1)) ' % '];
        text(25,75,string7,'Units','pixels')
        
        string8 = ['residency time = ' num2str(round2(fitting_results.(best_model).TTT2,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).PPP2*100,1e-1)) ' % '];
        text(25,50,string8,'Units','pixels')
        
        string9 = ['residency time = ' num2str(round2(fitting_results.(best_model).TTT3,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).PPP3*100,1e-1)) ' % '];
        text(25,25,string9,'Units','pixels')
        
        string9 = ['residency time = ' num2str(round2(fitting_results.(best_model).TTT4,1e-2)) ' sec. and ' ...
            num2str(round2(fitting_results.(best_model).PPP4*100,1e-1)) ' % '];
        text(25,25,string9,'Units','pixels')
        
        string10 = ['flag= ' num2str(fitting_results.(best_model).flag)];
        text(25,100,string10,'Units','pixels')
 
    end
    
    %-------------
    % cosmetic stuff
    if log_scale == 1
        set(gca,'xscale','log');
    elseif log_scale ==2
        set(gca,'yscale','log');
    end
    
    xlabel('tracks duration in sec')
    
    if classification_performed == 0
            title(['Duration distributions : ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+bestfit-', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 1
        title(['Duration distributions : ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution_T0removal+bestfit-', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    elseif classification_performed == 2
        title(['Duration distributions : ' given_set, ', ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+bestfit_', given_set,'_', name0, '-', name_embryo, '-', name1 , '-', name2, '.fig');
    end
    saveas(gcf,[save_stem namePlot]);
    
    if classification_performed == 0
        title(['Duration distributions : ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+bestfit-', name0, '-', name_embryo, '-', name1 , '-', name2, '.tif');
    elseif classification_performed == 1
        title(['Duration distributions : ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution_T0removal+bestfit-', name0, '-', name_embryo, '-', name1 , '-', name2, '.tif');
    elseif classification_performed == 2
         title(['Duration distributions : ' given_set, ', ', name1, ' and ', 10, name2, '(', best_model, ' for ', name_embryo, ' )' ]);
        namePlot = strcat('tracks_duration_distribution+bestfit_', given_set,'_', name0, '-', name_embryo, '-', name1 , '-', name2, '.tif');
    end
    saveas(gcf,[save_stem namePlot]);
    
    close(gcf)
    
end


end

