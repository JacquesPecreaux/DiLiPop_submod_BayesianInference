function [ fitting_results ] = to_calculate_errors_parameters_using_ProfileLikelihood...
    ( input,fitting_results,best_model,nb_embryos,name1,name2,save_stem,fval_mono,fval_double,fval_triple )

size_population2 = fitting_results.size_population ;

if strcmp(best_model,'MonoExpo')
    
    %lower range of parameter
    np = 100;  % nb of data that will be tested
    pm = 1 / (fitting_results.MonoExpo.T + 0.05);  % min of p = 1/T
    pM = 1 / ( fitting_results.MonoExpo.T - 0.05); % max of p = 1/T
    
    p_mono_ = linspace(pm, pM, np);
    
    LLt_fixedT = zeros(np,1);
    for i=1:np
        p_fixed = p_mono_(i);
        LLt_fixedT(i,1) = - Loglikelihood2_total(p_fixed, input,'MonoExpo', @simple_exp2_beta_norm, nb_embryos, size_population2);
    end
    
    % find interval of confidence
    confidence_level = 0.95;
    degree_of_freedom = 1;
    %log Lrest >= log Labs - 3.84/2
    %chi2inv(0.95,1) = 3.84
    max_LLt_mono_minus05 =  -fval_mono - chi2inv(confidence_level,degree_of_freedom); %max_LLt_mono-2;
    index_above = find(LLt_fixedT >= max_LLt_mono_minus05);
    if ~isempty(index_above)
        T_upperBound = 1/p_mono_(min(index_above));
        T_lowerBound = 1/p_mono_(max(index_above));
    else
        T_upperBound = NaN;
        T_lowerBound = NaN;
    end
    
    %display results
    figure
    plot(1./p_mono_(1,:),LLt_fixedT(:,1) )
    hold all
    plot(1./p_mono_(1,:),max_LLt_mono_minus05.*ones(1,np),'--k')
    plot(fitting_results.MonoExpo.T,-fval_mono,'bo');
    plot(T_upperBound,max_LLt_mono_minus05,'gx');
    plot(T_lowerBound,max_LLt_mono_minus05,'gx');
    xlabel('T (s)')
    ylabel('maximum likelihood')
    string = [ 'Confidence interval: ', num2str(round2(T_lowerBound,2e-2)), ' - ',  num2str(round2(T_upperBound,2e-2)) ,' s'];
    text(25,50,string,'Units','pixels')
    namePlot = strcat('maximum_likelihood_MonoExpo_T', name1, '_', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    namePlot = strcat('maximum_likelihood_MonoExpo_T', name1, '_', name2, '.tif');
    saveas(gcf,[save_stem namePlot]);
    
    fitting_results.MonoExpo.T_lowerBound = T_lowerBound;
    fitting_results.MonoExpo.T_upperBound = T_upperBound;
    fitting_results.MonoExpo.T_se = abs(T_upperBound-T_lowerBound)/2;
    
    
elseif strcmp(best_model,'DoubleExpo')
    
    % Set constant T1 in the model and calculate the max likelihood with
    % this new model.
    % repeat this for T1 value set in arange between -X and +X from optimal
    % value
    
    np = 100;  % nb of data that will be tested
    p2m = 1 / ( fitting_results.DoubleExpo.T1 + 0.1);  % min of p2 = 1/T1
    p2M = 1 / ( fitting_results.DoubleExpo.T1 - 0.1); % max of p2 = 1/T1
    p2_double_ = linspace(p2m, p2M, np);
    
    Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);
    
    LLt_fixedT1 = zeros(np,1);
    for i=1:np
        p2_fixed = p2_double_(i);
        [x2_, fval, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedT0', @double_exp2_beta_norm_fixedT0, ...
            nb_embryos, size_population2,1/p2_fixed),  [0.5 0.5], [], [], [], [],[0 0], [], [], Options);
        
        LLt_fixedT1(i,1) = -fval;
        p1_double__(i) = x2_(1);
        p3_double__(i) = x2_(2);
    end
    
    % find interval of confidence
    confidence_level = 0.95;
    degree_of_freedom = 3-1; % or 3??
    %log Lrest >= log Labs - 7.81/2
    %chi2inv(0.95,3) = 7.81
    %chi2inv(0.95,2) = 5.99
    max_LLt_double_minus05 =  -fval_double - chi2inv(confidence_level,degree_of_freedom); %max_LLt_double-3;
    index_above = find(LLt_fixedT1 >= max_LLt_double_minus05);
    if ~isempty(index_above)
        T1_upperBound = 1/p2_double_(min(index_above));
        T1_lowerBound = 1/p2_double_(max(index_above));
    else
        T1_upperBound = NaN;
        T1_lowerBound = NaN;
    end
    
    % display results
    figure
    plot(1./p2_double_(1,:),LLt_fixedT1(:,1) )
    hold all
    plot(1./p2_double_(1,:),max_LLt_double_minus05.*ones(1,np),'--k');
    plot(fitting_results.DoubleExpo.T1,-fval_double,'bo');
    plot(T1_upperBound,max_LLt_double_minus05,'gx');
    plot(T1_lowerBound,max_LLt_double_minus05,'gx');
    string = [ 'Confidence interval for T1: ', num2str(round2(T1_lowerBound,2e-2)), ' - ',  num2str(round2(T1_upperBound,2e-2)) ,' s'];
    text(25,50,string,'Units','pixels')
    xlabel('T1 (s)')
    ylabel('maximum likelihood')
    namePlot = strcat('maximum_likelihood_DoubleExpo_T1', name1, '_', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    namePlot = strcat('maximum_likelihood_DoubleExpo_T1', name1, '_', name2, '.tif');
    saveas(gcf,[save_stem namePlot]);
    
    fitting_results.DoubleExpo.T1_lowerBound = T1_lowerBound;
    fitting_results.DoubleExpo.T1_upperBound = T1_upperBound;
    fitting_results.DoubleExpo.T1_se = abs(T1_upperBound-T1_lowerBound)/2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set constant T2 in the model and calculate the max likelihood with
    % this new model.
    % repeat this for T2 value set in arange between -X and +X from optimal
    % value
    
    np = 100;  % nb of data that will be tested
    p3m = 1 / ( fitting_results.DoubleExpo.T2 + 0.25); % min of p2 = 1/T1
    p3M = 1 / ( fitting_results.DoubleExpo.T2 - 0.25);% max of p2 = 1/T1
    p3_double_ = linspace(p3m, p3M, np);
    
    Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);
    
    LLt_fixedT2 = zeros(np,1);
    for i=1:np
        p3_fixed = p3_double_(i);
        [x2_, fval, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedT0', @double_exp2_beta_norm_fixedT0, ...
            nb_embryos, size_population2,1/p3_fixed),  [0.6 0.4], [], [], [], [],[0 0], [], [], Options);
        LLt_fixedT2(i,1) = -fval;
        p1_double__(i) = x2_(1);
        p2_double__(i) = x2_(2);
    end
    
    % find interval of confidence
    index_above = find(LLt_fixedT2 >= max_LLt_double_minus05);
    if ~isempty(index_above)
        T2_upperBound = 1/p3_double_(min(index_above));
        T2_lowerBound = 1/p3_double_(max(index_above));
    else
        T2_upperBound = NaN;
        T2_lowerBound = NaN;
    end
    
    % display results
    figure
    plot(1./p3_double_(1,:),LLt_fixedT2(:,1) )
    hold all
    plot(1./p3_double_(1,:),max_LLt_double_minus05.*ones(1,np),'--k')
    plot(fitting_results.DoubleExpo.T2,-fval_double,'bo');
    plot(T2_upperBound,max_LLt_double_minus05,'gx');
    plot(T2_lowerBound,max_LLt_double_minus05,'gx');
    string = [ 'Confidence interval for T2: ', num2str(round2(T2_lowerBound,2e-2)), ' - ',  num2str(round2(T2_upperBound,2e-2)) ,' s'];
    text(25,50,string,'Units','pixels')
    xlabel('T2 (s)')
    ylabel('maximum likelihood')
    namePlot = strcat('maximum_likelihood_DoubleExpo_T2', name1, '_', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    namePlot = strcat('maximum_likelihood_DoubleExpo_T2', name1, '_', name2, '.tif');
    saveas(gcf,[save_stem namePlot]);
    
    fitting_results.DoubleExpo.T2_lowerBound = T2_lowerBound;
    fitting_results.DoubleExpo.T2_upperBound = T2_upperBound;
    fitting_results.DoubleExpo.T2_se = abs(T2_upperBound-T2_lowerBound)/2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set constant P1 in the model and calculate the max likelihood with
    % this new model.
    % repeat this for P1 value set in arange between -X and +X from optimal
    % value
    
    np = 100;  % nb of data that will be tested
    p1m = fitting_results.DoubleExpo.P1 + 0.1; % min of P1
    p1M = fitting_results.DoubleExpo.P1 - 0.1; % max of P1
    p1_double_ = linspace(p1m, p1M, np);
    
    Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);
    
    LLt_fixedP1 = zeros(np,1);
    for i=1:np
        p1_fixed = p1_double_(i);
        [x2_, fval, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedP0', @double_exp2_beta_norm_fixedP0, ...
            nb_embryos, size_population2,[],p1_fixed),  [0.6 0.4], [], [], [], [],[0 0], [], [], Options);
        
        LLt_fixedP1(i,1) = -fval;
        p2_double__(i) = x2_(1);
        p3_double__(i) = x2_(2);
    end
    
    % find interval of confidence
    index_above = find(LLt_fixedP1 >= max_LLt_double_minus05);
    if ~isempty(index_above)
        P1_lowerBound = p1_double_(min(index_above));
        P1_upperBound = p1_double_(max(index_above));
    else
        P1_upperBound = NaN;
        P1_lowerBound = NaN;
    end
    
    % display results
    figure
    plot(p1_double_(1,:),LLt_fixedP1(:,1) )
    hold all
    plot(p1_double_(1,:),max_LLt_double_minus05.*ones(1,np),'--k');
    plot(fitting_results.DoubleExpo.P1,-fval_double,'bo');
    plot(P1_lowerBound,max_LLt_double_minus05,'gx');
    plot(P1_upperBound,max_LLt_double_minus05,'gx');
    string = [ 'Confidence interval for P1: ', num2str(round2(P1_lowerBound,2e-2)*100), ' - ',  num2str(round2(P1_upperBound,2e-2)*100) ,'%'];
    text(25,50,string,'Units','pixels')
    xlabel('P1 (a.u.)')
    ylabel('maximum likelihood')
    namePlot = strcat('maximum_likelihood_DoubleExpo_P1', name1, '_', name2, '.fig');
    saveas(gcf,[save_stem namePlot]);
    namePlot = strcat('maximum_likelihood_DoubleExpo_P1', name1, '_', name2, '.tif');
    saveas(gcf,[save_stem namePlot]);
    
    fitting_results.DoubleExpo.P1_lowerBound = P1_lowerBound;
    fitting_results.DoubleExpo.P1_upperBound = P1_upperBound;
    fitting_results.DoubleExpo.P1_se = abs(P1_upperBound-P1_lowerBound)/2;
    fitting_results.DoubleExpo.P2_se = abs(P1_upperBound-P1_lowerBound)/2;
    
end

end

