function [ fitting_results ] = to_calculate_parameters_using_min_loglikelihood_sumEmbryo...
    ( input,fitting_results,models,nbEmbryo_givenCondition,name1,name2,save_stem)

%% get raw size population for each embryos that will be necessary for normalization

size_population_raw = 0;

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    name_embryo = ['embryo' num2str(iEmbryo)];
    size_population_raw =  size_population_raw + nansum( input.(name_embryo).data(:,2) );
    size_population.(name_embryo).raw = nansum( input.(name_embryo).data(:,2) );
    
end

size_population.total.raw = size_population_raw;
fitting_results.size_population = size_population;
size_population_mean = size_population_raw /nbEmbryo_givenCondition;
size_population2 = size_population; % to avoid conflict later

clear size_population_raw size_population


if ~isempty(find(ismember(models,'MonoExpo'),1))
    
    np = 500; % nb of data that will be tested
    pm = 0.2; %min of parameter value
    pM = 10; % max of parameter value
    
    p_mono = linspace(pm, pM, np);
    
    LLt_mono = zeros(np,1);
    
    for j = 1:np
        par = p_mono(j);
        LLt_mono(j) = Loglikelihood2_total(par, input,'MonoExpo', @simple_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
    end
    
    min_LLt_mono = min(LLt_mono);
    fi_min = find(LLt_mono == min_LLt_mono);
    fitting_results.MonoExpo.T_lik = 1 / p_mono(fi_min);
    
    figure
    plot(1./p_mono(1,:),LLt_mono(:,1) )
    hold all
    plot(1/p_mono(fi_min),min_LLt_mono,'ro');
    xlabel('T (s)')
    ylabel('maximum likelihood')
    title(['maximum likelihood : ' name1 ' and ' name2 ]);
    string1 = ['T= ' num2str(fitting_results.MonoExpo.T_lik) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_MonoExpo-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure
    plot(p_mono(1,:),LLt_mono(:,1) )
    hold all
    plot(p_mono(fi_min),min_LLt_mono,'ro');
    xlabel('1/T (s-1)')
    ylabel('maximum likelihood')
    title(['maximum likelihood : ' name1 ' and ' name2 ]);
    string1 = ['T= ' num2str(fitting_results.MonoExpo.T_lik) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_MonoExpo2-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    clear np pm pM p_mono LLT_mono par min_LLt_mono fi_min namePlot
    close all
    
end


nb_classes = 25;

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    
    np = 100;
    p1m = 0;
    p1M = 1;
    p2m = 1;
    p2M = 10;
    p3m = 0.2;
    p3M = 2;
    
    p1_double = linspace(p1m, p1M, np);
    p2_double = linspace(p2m, p2M, np);
    p3_double = linspace(p3m, p3M, np);
    
    par_default = [p1m p2m p3m];
    min_LLt_double = Loglikelihood2_total(par_default, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
    
    LLt_double = zeros(np, np, np);
    for i = 1:np
        for j=1:np
            for k=1:np
                par = [p1_double(i) p2_double(j) p3_double(k) ];
                LLt_double(i,j,k) = Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
                if LLt_double(i,j,k) < min_LLt_double
                    min_LLt_double = LLt_double(i,j,k);
                    p_double = par;
                    index_double = [i j k];
                end
            end
        end
    end
    
    [min_LLt_double,idx_LLt_double] = min(LLt_double(:));
    [a, b, c]=ind2sub(size(LLt_double),idx_LLt_double);
    fitting_results.DoubleExpo.P1_lik = p1_double(a);
    fitting_results.DoubleExpo.P2_lik = 1 - p1_double(a);
    fitting_results.DoubleExpo.T1_lik = 1 / p2_double(b);
    fitting_results.DoubleExpo.T2_lik = 1 / p3_double(c);
    
    %%
    
    %--------------------
    % variation of parameters 2 and 3
    LLt_double_23 = zeros(np,np);
    ifixe_1 = index_double(1);
    p1_double_ref = p1_double(ifixe_1);
    
    for i=1:np
        for j = 1:np
            LLt_double_23(i,j) = LLt_double(ifixe_1, i,j);
        end
    end
    
    figure,
    mesh(p2_double,p3_double,LLt_double_23)
    hold all
    xlabel('1/T1 (s-1)')
    ylabel('1/T2 (s-1)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['P1 fixed to ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par23_3D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(p2_double,p3_double,LLt_double_23,nb_classes)
    hold all
    plot(1/fitting_results.DoubleExpo.T1_lik,1/fitting_results.DoubleExpo.T2_lik,'or');
    xlabel('1/T1 (s-1)')
    ylabel('1/T2 (s-1)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['P1 fixed to ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par23_2D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    mesh(1./p2_double,1./p3_double,LLt_double_23);
    hold all
    xlabel('T1 (s)')
    ylabel('T2 (s)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['P1 fixed to ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par23_3D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(1./p2_double,1./p3_double,LLt_double_23,nb_classes)
    hold all
    plot(fitting_results.DoubleExpo.T1_lik,fitting_results.DoubleExpo.T2_lik,'or');
    xlabel('T1 (s)')
    ylabel('T2 (s)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['P1 fixed to ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
    text(25,100,string1,'Units','pixels')
    string2 = ['T1 = ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s and T2 = ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,50,string2,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par23_2D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    
    %----------------------------
    % variation of parameters 1 and 2
    
    LLt_double_12 = zeros(np,np);
    kfixe_3 = index_double(3);
    p3_double_ref = p3_double(kfixe_3);
    
    for i=1:np
        for j = 1:np
            LLt_double_12(i,j) = LLt_double(i, j, kfixe_3);
        end
    end
    
    figure,
    mesh(p1_double,p2_double,LLt_double_12)
    hold all
    xlabel('P1 (a.u.)')
    ylabel('1/T1 (s-1)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['T2 fixed to ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par12_3D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(p1_double,p2_double,LLt_double_12,nb_classes)
    hold all
    plot(fitting_results.DoubleExpo.P1_lik,1/fitting_results.DoubleExpo.T1_lik,'or');
    xlabel('P1 (s-1)')
    ylabel('1/T1 (s-1)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['T2 fixed to ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par12_2D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    mesh(p1_double,1./p2_double,LLt_double_12)
    hold all
    xlabel('P1 (a.u.)')
    ylabel('T1 (s)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['T2 fixed to ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par12_3D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(p1_double,1./p2_double,LLt_double_12,nb_classes)
    hold all
    plot(fitting_results.DoubleExpo.P1_lik,fitting_results.DoubleExpo.T1_lik,'or');
    xlabel('P1 (a.u.)')
    ylabel('T1 (s)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['T2 fixed to ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    string2 = ['T1 = ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s and P1 = ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
    text(25,50,string2,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par12_2D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    %----------------
    % variation of parameters 1 and 3
    
    LLt_double_13 = zeros(np,np);
    jfixe_2 = index_double(2);
    p2_double_ref = p2_double(jfixe_2);
    
    for i=1:np
        for j = 1:np
            LLt_double_13(i,j) = LLt_double(i, jfixe_2,j);
        end
    end
    
    figure,
    mesh(p1_double,p3_double,LLt_double_13)
    hold all
    xlabel('P1 (a.u.)')
    ylabel('1/T2 (s-1)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['T1 fixed to ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par13_3D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(p1_double,p3_double,LLt_double_13,nb_classes)
    hold all
    plot(fitting_results.DoubleExpo.P1_lik,1/fitting_results.DoubleExpo.T2_lik,'or');
    xlabel('P1 (a.u.)')
    ylabel('1/T2 (s-1)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['T1 fixed to ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par13_2D-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    mesh(p1_double,1./p3_double,LLt_double_13)
    hold all
    xlabel('P1 (a.u.)')
    ylabel('T2 (s)')
    zlabel('maximum likelihood')
    title([ name1 ' and ' name2 ]);
    string1 = ['T1 fixed to ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par13_3D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
    figure,
    contourf(p1_double,1./p3_double,LLt_double_13,nb_classes)
    hold all
    plot(fitting_results.DoubleExpo.P1_lik,fitting_results.DoubleExpo.T2_lik,'or');
    xlabel('P1 (s)')
    ylabel('T2 (s)')
    zlabel('maximum likelihood')
    colormap(gca,'jet')
    title([ name1 ' and ' name2 ]);
    string1 = ['T1 fixed to ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' s'];
    text(25,100,string1,'Units','pixels')
    string2 = ['P1 = ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' % and T2 = ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
    text(25,50,string2,'Units','pixels')
    namePlot = strcat('maximum_likelihood_DoubleExpo_par13_2D_bis-', name1 , '-', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    
%     %%
%     % new calculations of the likelihood for a given parameter value
%     
%     %--------------------
%     % variation of parameters 2 and 3
%     LLt_double_23_bis = zeros(np,np);
%     
%     par_default = [p1m p2m p3m];
%     min_LLt_double_23 = Loglikelihood2_total(par_default, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%     
%     for i=1:np
%         for j = 1:np
%                 par = [fitting_results.DoubleExpo.P1_lik p2_double(i) p3_double(j) ];
%                 LLt_double_23_bis(i,j) = Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%                 if LLt_double_23_bis(i,j) < min_LLt_double_23
%                     min_LLt_double_23 = LLt_double_23_bis(i,j);
%                     p_double_23 = par;
%                     index_double_23 = [i j];
%                 end
%         end
%     end
% 
%     [min_LLt_double_23,idx_LLt_double_23] = min(LLt_double_23_bis(:));
%     [a_, b_]=ind2sub(size(LLt_double_23_bis),idx_LLt_double_23);
%     
%     figure,
%     contourf(1./p2_double,1./p3_double,LLt_double_23_bis,nb_classes)
%     hold all
%     plot(1/p2_double(a_),1/p2_double(b_),'or');
%     xlabel('T1 (s)')
%     ylabel('T2 (s)')
%     zlabel('maximum likelihood')
%     colormap(gca,'jet')
%     title([ name1 ' and ' name2 ]);
%     string1 = ['P1 fixed to ' num2str(round2(fitting_results.DoubleExpo.P1_lik*100,1e-2)) ' %'];
%     text(25,100,string1,'Units','pixels')
%     string2 = ['T1 = ' num2str(round2(1/p2_double(a_),1e-2)) ' s and T2 = ' num2str(round2(1/p2_double(b_),1e-2)) ' s'];
%     text(25,50,string2,'Units','pixels')
%     namePlot = strcat('maximum_likelihood_DoubleExpo_par23_2D_ter-', name1 , '-', name2, '.fig');
%     saveas(gcf,[save_stem namePlot]);
%     
%     
%     %----------------------------
%     % variation of parameters 1 and 2
%     LLt_double_12_bis = zeros(np,np);
%     
%     par_default = [p1m p2m p3m];
%     min_LLt_double_12 = Loglikelihood2_total(par_default, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%     
%     for i=1:np
%         for j = 1:np
%                 par = [p1_double(i) p2_double(j) 1/fitting_results.DoubleExpo.T2_lik ];
%                 LLt_double_12_bis(i,j) = Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%                 if LLt_double_12_bis(i,j) < min_LLt_double_12
%                     min_LLt_double_12 = LLt_double_12_bis(i,j);
%                     p_double_12 = par;
%                     index_double_12 = [i j];
%                 end
%         end
%     end
% 
%     [min_LLt_double_12,idx_LLt_double_12] = min(LLt_double_12_bis);
%     [a_, b_]=ind2sub(size(LLt_double_12_bis),idx_LLt_double_12);
%     
%     figure,
%     contourf(p1_double,1./p2_double,LLt_double_12_bis,nb_classes)
%     hold all
%     plot(p1_double(a_),1./p2_double(b_),'or');
%     xlabel('P1 (s)')
%     ylabel('T1 (s)')
%     zlabel('maximum likelihood')
%     colormap(gca,'jet')
%     title([ name1 ' and ' name2 ]);
%     string1 = ['T2 fixed to ' num2str(round2(fitting_results.DoubleExpo.T2_lik,1e-2)) ' s'];
%     text(25,100,string1,'Units','pixels')
%     string2 = ['P1 = ' num2str(round2(p1_double(a_)*100,1e-2)) ' s and T1 = ' num2str(round2(1/p2_double(b_),1e-2)) ' s'];
%     text(25,50,string2,'Units','pixels')
%     namePlot = strcat('maximum_likelihood_DoubleExpo_par12_2D_ter-', name1 , '-', name2, '.fig');
%     saveas(gcf,[save_stem namePlot]);
%     
%     %----------------
%     % variation of parameters 1 and 3
%     
%      LLt_double_13_bis = zeros(np,np);
%     
%     par_default = [p1m p2m p3m];
%     min_LLt_double_13 = Loglikelihood2_total(par_default, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%     
%     for i=1:np
%         for j = 1:np
%                 par = [ p1_double(i) 1/fitting_results.DoubleExpo.T1_lik p3_double(j) ];
%                 LLt_double_13_bis(i,j) = Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
%                 if LLt_double_13_bis(i,j) < min_LLt_double_13
%                     min_LLt_double_13 = LLt_double_13_bis(i,j);
%                     p_double_13 = par;
%                     index_double_13 = [i j];
%                 end
%         end
%     end
% 
%     [min_LLt_double_13,idx_LLt_double_13] = min(LLt_double_13_bis);
%     [a_, b_]=ind2sub(size(LLt_double_13_bis),idx_LLt_double_13);
%     
%     figure,
%     contourf(p1_double,1./p3_double,LLt_double_13_bis,nb_classes)
%     hold all
%     plot(p1_double(a_),1/p3_double(b_),'or');
%     xlabel('P1 (s)')
%     ylabel('T2 (s)')
%     zlabel('maximum likelihood')
%     colormap(gca,'jet')
%     title([ name1 ' and ' name2 ]);
%     string1 = ['T1 fixed to ' num2str(round2(fitting_results.DoubleExpo.T1_lik,1e-2)) ' %'];
%     text(25,100,string1,'Units','pixels')
%     string2 = ['P1 = ' num2str(round2(p1_double(a_)*100,1e-2)) ' s and T2 = ' num2str(round2(1/p3_double(b_),1e-2)) ' s'];
%     text(25,50,string2,'Units','pixels')
%     namePlot = strcat('maximum_likelihood_DoubleExpo_par13_2D_ter-', name1 , '-', name2, '.fig');
%     saveas(gcf,[save_stem namePlot]);
%         
%     clear np p1m p1M p2m p2M p3m p3M p1_double p2_double p3_double
%     clear par_default min_LLt_double LLt_double p_double index_double idx_LLt_double a b c
%     clear LLT_double_23 ifixe_1 p1_double_ref
%     clear LLT_double_12 kfixe_3 p3_double_ref
%     clear LLT_double_13 jfixe_2 p2_double_ref
%     close all
    
end



if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    
    np = 100;
    p1m = 0.2;
    p1M = 10;
    p2m = 0;
    p2M = 4;
    
    p1_monoS = linspace(p1m, p1M, np);
    p2_monoS = linspace(p2m, p2M, np);
    
    par_default = [p1m p2m];
    min_LLt_monoS = Loglikelihood2_total(par_default, input,'MonoExpo_stretched', @simple_exp2_stretched_beta_norm, nbEmbryo_givenCondition, size_population2);
    
    LLt_monoS = zeros(np, np);
    for i = 1:np
        for j=1:np
            par = [p1_double(i) p2_double(j) ];
            LLt_monoS(i,j) = Loglikelihood2_total(par, input,'MonoExpo_stretched', @simple_exp2_stretched_beta_norm, nbEmbryo_givenCondition, size_population2);
            if LLt_monoS(i,j) < min_LLt_monoS
                min_LLt_monoS = LLt_monoS(i,j);
                p_monoS = par;
                index_monoS = [i j];
            end
        end
    end
    
    [min_LLt_monoS,idx_LLt_monoS] = min(LLt_monoS(:));
    [a, b]=ind2sub(size(LLt_monoS),idx_LLt_monoS);
    fitting_results.MonoExpo_stretched.Ts_lik = 1/p1_monoS(a);
    fitting_results.MonoExpo_stretched.power_lik = p2_monoS(b);

end

if ~isempty(find(ismember(models,'TripleExpo'),1))
    
    np = 100;
    p1m = 0;
    p1M = 1;
    p2m = 1;
    p2M = 10;
    p3m = 0;
    p3M = 1;
    p4m = 0.2;
    p4M = 2;
    p5m = 0.1;
    p5M = 1;
    
    p1_triple = linspace(p1m, p1M, np);
    p2_triple = linspace(p2m, p2M, np);
    p3_triple = linspace(p3m, p3M, np);
    p4_triple = linspace(p4m, p4M, np);
    p5_triple = linspace(p5m, p5M, np);
    
    par_default = [p1m p2m p3m p4m p5m];
    min_LLt_triple = Loglikelihood2_total(par_default, input,'TripleExpo', @triple_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
    
    LLt_triple = zeros(np, np, np, np, np);
    for i = 1:np
        for j=1:np
            for k=1:np
                for l = 1 : np
                    for m = 1 : np
                        par = [p1_triple(i) p2_triple(j) p3_triple(k) p4_triple(l) p5_triple(m)];
                            LLt_triple(i,j,k,l,m) = Loglikelihood2_total(par, input,'TripleExpo', @triple_exp2_beta_norm, nbEmbryo_givenCondition, size_population2);
                            if LLt_triple(i,j,k,l,m) < min_LLt_triple
                                min_LLt_triple = LLt_triple(i,j,k,l,m);
                                p_triple = par;
                                index_triple = [i j k l m];
                            end
                    end
                end
            end
        end
    end
    
    [min_LLt_triple,idx_LLt_triple] = min(LLt_triple(:));
    [a, b, c, d, e]=ind2sub(size(LLt_triple),idx_LLt_triple);
    fitting_results.TripleExpo.P1_lik = p1_triple(a);
    fitting_results.TripleExpo.P2_lik = p3_triple(c);
    fitting_results.TripleExpo.P3_lik = 1 - p1_triple(a) - p3_triple(c);
    fitting_results.TripleExpo.T1_lik = 1 / p2_double(b);
    fitting_results.TripleExpo.T2_lik = 1 / p4_double(d);
    fitting_results.TripleExpo.T3_lik = 1 / p5_double(e);
    

end



end

