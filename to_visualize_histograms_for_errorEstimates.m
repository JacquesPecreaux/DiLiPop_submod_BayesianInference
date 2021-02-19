function to_visualize_histograms_for_errorEstimates( timing_phase,blastomere,T_allSimulations,T1_allSimulations,T2_allSimulations,P1_allSimulations,P2_allSimulations, ...
    T2_2_allSimulations,P2_1_allSimulations,P2_2_allSimulations,TT1_allSimulations,TT2_allSimulations,TT3_allSimulations,PP1_allSimulations,PP2_allSimulations,PP3_allSimulations,...
    Ts_allSimulations,power_allSimulations,nbEmbryo_givenCondition,save_stem_error,best_model,nb_bins,display_plots)


if display_plots == 1
    
    if nbEmbryo_givenCondition == 1
        
        if strcmp(best_model,'MonoExpo')
            
            figure
            hist(T_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T : MonoExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T-MonoExpo-', blastomere,'-', timing_phase, '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'DoubleExpo')
            
            figure
            hist(T1_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T1 : DoubleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T1-DoubleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(T2_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T2 : DoubleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T2-DoubleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P1_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P1 : DoubleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P1-DoubleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P2_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P2 : DoubleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P2-DoubleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
 
        elseif strcmp(best_model,'DoubleExpo_fixedT0')

            figure
            hist(T2_2_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T2 : DoubleExpo_fixedT0 :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T2-DoubleExpo_fixedT0-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P2_1_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P1 : DoubleExpo_fixedT0 :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P1-DoubleExpo_fixedT0-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P2_2_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P2 : DoubleExpo_fixedT0 :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P2-DoubleExpo_fixedT0-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'TripleExpo')
            figure
            hist(TT1_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT1 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT1-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(TT2_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT2 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT2-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(TT3_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT3 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT3-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(PP1_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP1 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP1-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,fullfile(save_stem_error, namePlot));
            
            figure
            hist(PP2_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP2 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP2-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(PP3_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP3 : TripleExpo :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP3-TripleExpo-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'MonoExpo_stretched')
            
            figure
            hist(Ts_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T : MonoExpo_stretched :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T-MonoExpoStretched-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(power_allSimulations, nb_bins);
            xlabel('power (a.u.)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['power : MonoExpo stretched :' blastomere ' and ' timing_phase ]);
            namePlot = strcat('Histogram_ErrorsEstimate_power-MonoExpoStretched-', blastomere,'-', timing_phase,  '.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end
        
    else
        
        if strcmp(best_model,'MonoExpo')
            
            figure
            hist(T_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T : MonoExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T-MonoExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'DoubleExpo')
            
            figure
            hist(T1_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T1 : DoubleExpo :'  ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T1-DoubleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(T2_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T2 : DoubleExpo :'  ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T2-DoubleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P1_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P1 : DoubleExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P1-DoubleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(P2_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['P2 : DoubleExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_P2-DoubleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'TripleExpo')
            figure
            hist(TT1_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT1 : TripleExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT1-TripleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(TT2_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT2 : TripleExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT2-TripleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(TT3_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['TT3 : TripleExpo :'  ]);
            namePlot = strcat('Histogram_ErrorsEstimate_TT3-TripleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(PP1_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP1 : TripleExpo :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP1-TripleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(PP2_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP2 : TripleExpo :'  ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP2.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(PP3_allSimulations, nb_bins);
            xlabel('population fraction (%)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['PP3 : TripleExpo :'  ]);
            namePlot = strcat('Histogram_ErrorsEstimate_PP3-TripleExpo.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        elseif strcmp(best_model,'MonoExpo_stretched')
            
            figure
            hist(Ts_allSimulations, nb_bins);
            xlabel('tracks duration (s)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['T : MonoExpo_stretched :' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_T-MonoExpoStretched.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
            figure
            hist(power_allSimulations, nb_bins);
            xlabel('power (a.u.)','FontSize',10);
            ylabel('Counts','FontSize',10);
            title(['power : MonoExpo stretched' ]);
            namePlot = strcat('Histogram_ErrorsEstimate_power-MonoExpoStretched.fig');
            saveas(gcf,[save_stem_error namePlot]);
            
        end
        
        
    end
    
end

close all

end

