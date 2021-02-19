function [ input,fitting_results,fval_mono,fval_double,fval_triple] = to_calculate_parameters_using_fmincon_sum_mle...
    ( input,models,nbEmbryo_givenCondition,simulation,fixed_short_lifetime,fixed_short_percent  )

% main script to get the parameters of each fitting model by looking at the
% sum of the likelihood of the different embryos, with the same model
% parameters kept to fit the different embryos

% for that, necessary to normalized the distribution + choice to write the
% model equation in the simplest way ( no denominator: exp(-beta*t) ) since
% we have to do derivate at least to get the errors on the fitting parameters



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


%% choice of the options for the fitting

Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);
fitting_results.options = Options;

%% fit with mono expo as models

if ~isempty(find(ismember(models,'MonoExpo'),1))
    
 %   [x, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'MonoExpo', @simple_exp2_beta_norm, nbEmbryo_givenCondition, size_population2), ...
  %      [1], [], [], [], [],[0.1], [1000],[],Options);
    [x, fval_mono, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'MonoExpo', @simple_exp2_beta_norm, nbEmbryo_givenCondition, size_population2), ...
        [1], [], [], [], [],[0], [],[],Options);
    
    fitting_results.MonoExpo.parameters = x;
    fitting_results.MonoExpo.nparam = 1;
    fitting_results.MonoExpo.T = 1/x(1);
    fitting_results.MonoExpo.flag = exitflag;
    
    %---------------
    % caluclation of R, the denominator in the simple_exp2_beta_norm
    % function due to the normalization
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        xdata=input.(name_embryo).data(:,1);
        R.(name_embryo).mono = sum( exp(-xdata.*x(1)) );
        clear xdata
    end
    fitting_results.MonoExpo.R_normalization = R;
    %----------------
    % calculation fo the size of the population after fitting for each
    % embryo and the total one
    size_population_tot = 0;
    size_population_norm = 0;
    size_population_normA = 0;
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        size_population_tot = size_population_tot + double( floor(vpa(sum(simple_exp2_beta_norm(x,input.(name_embryo).data,R.(name_embryo).mono, ...
            size_population2.(name_embryo).raw)),3)) );
        size_population.(name_embryo) = double( floor(vpa(sum(simple_exp2_beta_norm(x,input.(name_embryo).data, R.(name_embryo).mono,...
            size_population2.(name_embryo).raw)),3)) );
        if simulation == 0
            size_population_normalized.(name_embryo) = size_population.(name_embryo) / input.(name_embryo).duration_phase;
            size_population_norm = size_population_norm + size_population_normalized.(name_embryo) ;
            size_population_normalizedA.(name_embryo) = size_population_normalized.(name_embryo)*60 / input.(name_embryo).area;
            size_population_normA = size_population_normA + size_population_normalizedA.(name_embryo) ;
        end
    end
    size_population.total = double(size_population_tot);
    fitting_results.MonoExpo.size_population = size_population;
    if simulation == 0
        size_population_normalized.mean = size_population_norm / nbEmbryo_givenCondition;
        fitting_results.MonoExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
        size_population_normalizedA.mean = size_population_normA / nbEmbryo_givenCondition;
        fitting_results.MonoExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTS/min/um2
    end
    %-----------------
    % calculate the function value (mi) at each bins and the residues
    % between mi and ci
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        fitting_results.MonoExpo.(name_embryo).simple = simple_exp2_beta_norm(x,input.(name_embryo).data,R.(name_embryo).mono,size_population2.(name_embryo).raw);
        residues = (fitting_results.MonoExpo.(name_embryo).simple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.MonoExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
    end
    clear exitflag
    clear size_population
    clear residues
    clear R
    clear size_population_tot
    clear x
    clear size_population_normalized
    clear size_population_norm
    clear size_population_normalizedA
    clear size_population_normA
    
else
    fval_mono = NaN;
end


%% fit with double expo as model

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    
 %   [x2, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2), ...
 %       [0.5 2 0.5], [], [], [], [],[], [], [], Options);
    [x2, fval_double, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo', @double_exp2_beta_norm, nbEmbryo_givenCondition, size_population2), ...
        [0.5 2 0.5], [], [], [], [],[0 0 0], [], [], Options);
    
    fitting_results.DoubleExpo.parameters = x2;
    fitting_results.DoubleExpo.nparam = 3;
    if 1/x2(2) < 1/x2(3)
        fitting_results.DoubleExpo.P1 = x2(1);
        fitting_results.DoubleExpo.T1 = 1/x2(2);
        fitting_results.DoubleExpo.P2 = 1 - x2(1);
        fitting_results.DoubleExpo.T2 = 1/x2(3);
        fitting_results.DoubleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double(1,1) = sum( exp(-xdata.*x2(2)) );
            R.(name_embryo).double(2,1) = sum( exp(-xdata.*x2(3)) );
            clear xdata
        end
        fitting_results.DoubleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(double_exp2_beta_norm(x2,input.(name_embryo).data,R.(name_embryo).double(1,1),...
                R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],...
                input.(name_embryo).data,R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],...
                input.(name_embryo).data,R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm(x2,input.(name_embryo).data,...
                R.(name_embryo).double(1,1),R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],...
                input.(name_embryo).data,R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],...
                input.(name_embryo).data,R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo.(name_embryo).double = double_exp2_beta_norm(x2,input.(name_embryo).data,...
                R.(name_embryo).double(1,1),R.(name_embryo).double(2,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.DoubleExpo.(name_embryo).double(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
        
    elseif 1/x2(3) < 1/x2(2)
        
        fitting_results.DoubleExpo.P1 = 1-x2(1);
        fitting_results.DoubleExpo.T1 = 1/x2(3);
        fitting_results.DoubleExpo.P2 = x2(1);
        fitting_results.DoubleExpo.T2 = 1/x2(2);
        fitting_results.DoubleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double(1,1) = sum( exp(-xdata.*x2(3)) );
            R.(name_embryo).double(2,1) = sum( exp(-xdata.*x2(2)) );
            clear xdata
        end
        fitting_results.DoubleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(double_exp2_beta_norm(x2,input.(name_embryo).data,...
                R.(name_embryo).double(2,1),R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],...
                input.(name_embryo).data,R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],...
                input.(name_embryo).data,R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm(x2,input.(name_embryo).data,...
                R.(name_embryo).double(2,1),R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],input.(name_embryo).data,...
                R.(name_embryo).double(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],input.(name_embryo).data,...
                R.(name_embryo).double(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2        
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo.(name_embryo).double = double_exp2_beta_norm(x2,input.(name_embryo).data,...
                R.(name_embryo).double(2,1),R.(name_embryo).double(1,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.DoubleExpo.(name_embryo).double(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot1 size_population_tot2 percent_population_tot1 percent_population_tot2
    clear residues
    clear size_population percentage_population
    clear R
    clear x2
    clear size_population_norm1 size_population_norm2 size_population_normalized
    clear size_population_normA1 size_population_normA2 size_population_normalizedA size_population_normA_tot
    
else
    fval_double = NaN;
end


%% fit with double expo with fixed T0 as model

if ~isempty(find(ismember(models,'DoubleExpo_fixedT0'),1)) && exist('fixed_short_lifetime','var') == 1 && ~isempty(fixed_short_lifetime)
    
 %   [x2_, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedT0', @double_exp2_beta_norm_fixedT0, nbEmbryo_givenCondition, size_population2,fixed_short_lifetime), ...
  %      [0.5 0.5], [], [], [], [],[0 0.2], [1 1/(fixed_short_lifetime+0.1)], [], Options);
     [x2_, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedT0', @double_exp2_beta_norm_fixedT0, nbEmbryo_givenCondition, size_population2,fixed_short_lifetime), ...
        [0.5 0.5], [], [], [], [],[0 0], [], [], Options);
     
    fitting_results.DoubleExpo_fixedT0.parameters = x2_;
    fitting_results.DoubleExpo_fixedT0.nparam = 2;
    
    if fixed_short_lifetime < 1/x2_(2)
        
        fitting_results.DoubleExpo_fixedT0.P1 = x2_(1);
        fitting_results.DoubleExpo_fixedT0.T1 = fixed_short_lifetime;
        fitting_results.DoubleExpo_fixedT0.P2 = 1 - x2_(1);
        fitting_results.DoubleExpo_fixedT0.T2 = 1/x2_(2);
        fitting_results.DoubleExpo_fixedT0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm_fixedT0
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double_(1,1) = sum( exp(-xdata./fixed_short_lifetime) );
            R.(name_embryo).double_(2,1) = sum( exp(-xdata.*x2_(2)) );
            clear xdata
        end
        fitting_results.DoubleExpo_fixedT0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,R.(name_embryo).double_(1,1),...
                R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
            size_population_tot1 = size_population_tot1 + double( x2_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x2_(1))*floor(vpa(sum(simple_exp2_beta_norm([x2_(2)],...
                input.(name_embryo).data,R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,...
                R.(name_embryo).double_(1,1),R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
            size_population.(name_embryo).pop1 = double( x2_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( (1-x2_(1))*floor(vpa(sum(simple_exp2_beta_norm([x2_(2)],...
                input.(name_embryo).data,R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo_fixedT0.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo_fixedT0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo_fixedT0.(name_embryo).double_ = double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,...
                R.(name_embryo).double_(1,1),R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime);
            residues = (fitting_results.DoubleExpo_fixedT0.(name_embryo).double_(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo_fixedT0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
       
    elseif 1/x2_(2) < fixed_short_lifetime
        
        fitting_results.DoubleExpo_fixedT0.P1 = 1-x2_(1);
        fitting_results.DoubleExpo_fixedT0.T1 = 1/x2_(2);
        fitting_results.DoubleExpo_fixedT0.P2 = x2_(1);
        fitting_results.DoubleExpo_fixedT0.T2 = 1/fixed_short_lifetime;
        fitting_results.DoubleExpo_fixedT0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double_(1,1) = sum( exp(-xdata.*x2_(2)) );
            R.(name_embryo).double_(2,1) = sum( exp(-xdata./fixed_short_lifetime) );
            clear xdata
        end
        fitting_results.DoubleExpo_fixedT0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,...
                R.(name_embryo).double_(2,1),R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
            size_population_tot1 = size_population_tot1 + double( (1-x2_(1))*floor(vpa(sum(simple_exp2_beta_norm([x2_(2)],...
                input.(name_embryo).data,R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x2_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,...
                R.(name_embryo).double_(2,1),R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x2_(1))*floor(vpa(sum(simple_exp2_beta_norm([x2_(2)],input.(name_embryo).data,...
                R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x2_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double_(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo_fixedT0.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo_fixedT0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2        
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo_fixedT0.(name_embryo).double_ = double_exp2_beta_norm_fixedT0(x2_,input.(name_embryo).data,...
                R.(name_embryo).double_(2,1),R.(name_embryo).double_(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime);
            residues = (fitting_results.DoubleExpo_fixedT0.(name_embryo).double_(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo_fixedT0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot1 size_population_tot2 percent_population_tot1 percent_population_tot2
    clear residues
    clear size_population percentage_population
    clear R
    clear x2_
    clear size_population_norm1 size_population_norm2 size_population_normalized
    clear size_population_normA1 size_population_normA2 size_population_normalizedA size_population_normA_tot
    
end
   

%% fit with double expo with fixed T0P0 as model

if ~isempty(find(ismember(models,'DoubleExpo_fixedT0P0'),1)) && exist('fixed_short_lifetime','var') == 1 && ~isempty(fixed_short_lifetime) ...
        && exist('fixed_short_percent','var') == 1 && ~isempty(fixed_short_percent)
    
    [x2__, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'DoubleExpo_fixedT0P0', ...
        @double_exp2_beta_norm_fixedT0P0, nbEmbryo_givenCondition, size_population2,fixed_short_lifetime,fixed_short_percent), ...
        [0.5], [], [], [], [],[0], [], [], Options);
    
    fitting_results.DoubleExpo_fixedT0P0.parameters = x2__;
    fitting_results.DoubleExpo_fixedT0P0.nparam = 1;
    
    if fixed_short_lifetime < 1/x2__(1)
        
        fitting_results.DoubleExpo_fixedT0P0.P1 = fixed_short_percent;
        fitting_results.DoubleExpo_fixedT0P0.T1 = fixed_short_lifetime;
        fitting_results.DoubleExpo_fixedT0P0.P2 = 1 - fixed_short_percent;
        fitting_results.DoubleExpo_fixedT0P0.T2 = 1/x2__(1);
        fitting_results.DoubleExpo_fixedT0P0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm_fixedT0P0
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double__(1,1) = sum( exp(-xdata./fixed_short_lifetime) );
            R.(name_embryo).double__(2,1) = sum( exp(-xdata.*x2__(1)) );
            clear xdata
        end
        fitting_results.DoubleExpo_fixedT0P0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,R.(name_embryo).double__(1,1),...
                R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population_tot1 = size_population_tot1 + double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x2__(1)],...
                input.(name_embryo).data,R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,...
                R.(name_embryo).double__(1,1),R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population.(name_embryo).pop1 = double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( (1-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x2__(1)],...
                input.(name_embryo).data,R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo_fixedT0P0.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo_fixedT0P0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0P0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0P0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo_fixedT0P0.(name_embryo).double__ = double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,...
                R.(name_embryo).double__(1,1),R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent);
            residues = (fitting_results.DoubleExpo_fixedT0P0.(name_embryo).double__(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo_fixedT0P0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
       
    elseif 1/x2__(1) < fixed_short_lifetime
        
        fitting_results.DoubleExpo_fixedT0P0.P1 = 1-fixed_short_percent;
        fitting_results.DoubleExpo_fixedT0P0.T1 = 1/x2__(1);
        fitting_results.DoubleExpo_fixedT0P0.P2 = fixed_short_percent;
        fitting_results.DoubleExpo_fixedT0P0.T2 = 1/fixed_short_lifetime;
        fitting_results.DoubleExpo_fixedT0P0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).double__(1,1) = sum( exp(-xdata.*x2__(1)) );
            R.(name_embryo).double__(2,1) = sum( exp(-xdata./fixed_short_lifetime) );
            clear xdata
        end
        fitting_results.DoubleExpo_fixedT0P0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,...
                R.(name_embryo).double__(2,1),R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population_tot1 = size_population_tot1 + double( (1-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x2__(1)],...
                input.(name_embryo).data,R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,...
                R.(name_embryo).double__(2,1),R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population.(name_embryo).pop1 = double( (1-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x2__(1)],input.(name_embryo).data,...
                R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).double__(2,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.DoubleExpo_fixedT0P0.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.DoubleExpo_fixedT0P0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0P0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.DoubleExpo_fixedT0P0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2        
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.DoubleExpo_fixedT0P0.(name_embryo).double__ = double_exp2_beta_norm_fixedT0P0(x2__,input.(name_embryo).data,...
                R.(name_embryo).double__(2,1),R.(name_embryo).double__(1,1),size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent);
            residues = (fitting_results.DoubleExpo_fixedT0P0.(name_embryo).double__(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.DoubleExpo_fixedT0P0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot1 size_population_tot2 percent_population_tot1 percent_population_tot2
    clear residues
    clear size_population percentage_population
    clear R
    clear x2__
    clear size_population_norm1 size_population_norm2 size_population_normalized
    clear size_population_normA1 size_population_normA2 size_population_normalizedA size_population_normA_tot
    
end
 

%% fit with mono expo stretched as model

if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    
 %   [xS, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'MonoExpo_stretched', @simple_exp2_stretched_beta_norm, nbEmbryo_givenCondition,size_population2), ...
 %       [1 1],[], [], [], [],[0.1 0], [1000 10], [],Options);
        [xS, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'MonoExpo_stretched', @simple_exp2_stretched_beta_norm, nbEmbryo_givenCondition,size_population2), ...
        [1 1],[], [], [], [],[], [], [],Options);
    
    fitting_results.MonoExpo_stretched.parameters = xS;
    fitting_results.MonoExpo_stretched.nparam = 2;
    fitting_results.MonoExpo_stretched.Ts = 1/xS(1);
    fitting_results.MonoExpo_stretched.power = xS(2);
    fitting_results.MonoExpo_stretched.flag = exitflag;
    
    %---------------
    % caluclation of R, the denominator in the simple_exp2_beta_norm
    % function due to the normalization
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        xdata=input.(name_embryo).data(:,1);
        R.(name_embryo).monoS = sum( exp( -(xdata.*xS(1)).^(1/xS(2)) ) );
        clear xdata
    end
    fitting_results.MonoExpo_stretched.R_normalization = R;
    %----------------
    % calculation fo the size of the population after fitting for each
    % embryo and the total one
    size_population_tot = 0;
    size_population_norm = 0;
    size_population_normA = 0;
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        size_population_tot = size_population_tot + double( floor(vpa(sum(simple_exp2_stretched_beta_norm(xS,input.(name_embryo).data,R.(name_embryo).monoS,...
            size_population2.(name_embryo).raw)),3)) );
        size_population.(name_embryo) = double( floor(vpa(sum(simple_exp2_stretched_beta_norm(xS,input.(name_embryo).data, R.(name_embryo).monoS,...
            size_population2.(name_embryo).raw)),3)) );
        if simulation == 0
            size_population_normalized.(name_embryo) = size_population.(name_embryo) / input.(name_embryo).duration_phase;
            size_population_norm = size_population_norm + size_population_normalized.(name_embryo);
            size_population_normalizedA.(name_embryo) = size_population_normalized.(name_embryo)*60 / input.(name_embryo).area;
            size_population_normA = size_population_normA + size_population_normalizedA.(name_embryo);            
        end
    end
    size_population.total = double(size_population_tot);
    fitting_results.MonoExpo_stretched.size_population = size_population;
    if simulation == 0
        size_population_normalized.mean = size_population_norm / nbEmbryo_givenCondition;
        fitting_results.MonoExpo_stretched.size_population_normalized_time = size_population_normalized; % nbMTs / second
        size_population_normalizedA.mean = size_population_normA / nbEmbryo_givenCondition;
        fitting_results.MonoExpo_stretched.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2        
    end
    %-----------------
    % calculate the function value (mi) at each bins and the residues
    % between mi and ci
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        fitting_results.MonoExpo_stretched.(name_embryo).simpleS = simple_exp2_stretched_beta_norm(xS,input.(name_embryo).data,R.(name_embryo).monoS,...
            size_population2.(name_embryo).raw);
        residues = (fitting_results.MonoExpo_stretched.(name_embryo).simpleS(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.MonoExpo_stretched.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
    end
    clear exitflag
    clear size_population
    clear residues
    clear R
    clear size_population_tot
    clear xS
    clear size_population_norm size_population_normalized
    clear size_population_normA size_population_normalizedA
    
end


%% fit with triple expo as model

if ~isempty(find(ismember(models,'TripleExpo'),1))
    
 %   [x3, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'TripleExpo', @triple_exp2_beta_norm, nbEmbryo_givenCondition,size_population2),...
 %       [0.4 5 0.4 1.4 0.5 ], [1 0 1 0 0], [1], [], [],[0 2.1 0 0.5 0.1], [1 1000 1 4 0.67], [],Options);
    [x3, fval_triple, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'TripleExpo', @triple_exp2_beta_norm, nbEmbryo_givenCondition,size_population2),...
        [0.4 5 0.4 1.4 0.5 ], [1 0 1 0 0], [1], [], [],[0 0 0 0 0], [], [],Options);
    
    fitting_results.TripleExpo.parameters = x3;
    fitting_results.TripleExpo.nparam = 5;
    
    if 1/x3(2) < 1/x3(4) < 1/x3(5)
        
        fitting_results.TripleExpo.PP1 = x3(1);
        fitting_results.TripleExpo.TT1 = 1/x3(2);
        fitting_results.TripleExpo.PP2 = x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(4);
        fitting_results.TripleExpo.PP3 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(5);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(2)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(4)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(5)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;               
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
        
    elseif 1/x3(2) < 1/x3(5) < 1/x3(4)
        
        fitting_results.TripleExpo.PP1 = x3(1);
        fitting_results.TripleExpo.TT1 = 1/x3(2);
        fitting_results.TripleExpo.PP2 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(5);
        fitting_results.TripleExpo.PP3 = x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(4);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(2)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(5)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(4)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(3,1),R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    elseif 1/x3(4) < 1/x3(2) < 1/x3(5)
        
        fitting_results.TripleExpo.PP1 = x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(4);
        fitting_results.TripleExpo.PP2 = x3(1);
        fitting_results.TripleExpo.TT2 = 1/x3(2);
        fitting_results.TripleExpo.PP3 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(5);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(4)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(2)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(5)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;       
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2          
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    elseif 1/x3(4) < 1/x3(5) < 1/x3(2)
        
        fitting_results.TripleExpo.PP1 = x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(4);
        fitting_results.TripleExpo.PP2 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(5);
        fitting_results.TripleExpo.PP3 = x3(1);
        fitting_results.TripleExpo.TT3 = 1/x3(2);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(4)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(5)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(2)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2     
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    elseif 1/x3(5) < 1/x3(2) < 1/x3(4)
        
        fitting_results.TripleExpo.PP1 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(5);
        fitting_results.TripleExpo.PP2 = x3(1);
        fitting_results.TripleExpo.TT2 = 1/x3(2);
        fitting_results.TripleExpo.PP3 = x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(4);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(5)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(2)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(4)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min / um2        
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    elseif 1/x3(5) < 1/x3(4) < 1/x3(2)
        
        fitting_results.TripleExpo.PP1 = 1-x3(1)-x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(5);
        fitting_results.TripleExpo.PP2 = x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(4);
        fitting_results.TripleExpo.PP3 = x3(1);
        fitting_results.TripleExpo.TT3 = 1/x3(2);
        fitting_results.TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple(1,1) = sum( exp(-xdata.*x3(5)) );
            R.(name_embryo).triple(2,1) = sum( exp(-xdata.*x3(4)) );
            R.(name_embryo).triple(3,1) = sum( exp(-xdata.*x3(2)) );
            clear xdata
        end
        fitting_results.TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,R.(name_embryo).triple(1,1),...
                R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],input.(name_embryo).data,...
                R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],input.(name_embryo).data,...
                R.(name_embryo).triple(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        fitting_results.TripleExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        fitting_results.TripleExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2      
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo.(name_embryo).triple = triple_exp2_beta_norm(x3,input.(name_embryo).data,...
                R.(name_embryo).triple(3,1),R.(name_embryo).triple(2,1),R.(name_embryo).triple(1,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.TripleExpo.(name_embryo).triple(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot1 size_population_tot2 size_population_tot3 percent_population_tot1 percent_population_tot2 percent_population_tot3
    clear residues
    clear size_population percentage_population
    clear x3
    clear size_population_norm1 size_population_norm2 size_population_norm3 size_population_normalized
    clear size_population_normA1 size_population_normA2 size_population_normA3 size_population_normalizedA
    
else
    fval_triple = NaN;
end


%% fit with triple expo with fixed short lifetime as model

if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1)) && exist('fixed_short_lifetime','var') == 1 && ~isempty(fixed_short_lifetime)
    
    [x3_, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'TripleExpo_fixedT0', @triple_exp2_beta_norm_fixedT0, nbEmbryo_givenCondition, size_population2, fixed_short_lifetime), ...
         [0.4 0.4 1.5 0.5], [1 1 0 0], [1], [], [],[0 0 0 0], [], [], Options);
    
    fitting_results.TripleExpo_fixedT0.parameters = x3_;
    fitting_results.TripleExpo_fixedT0.nparam = 4;
    
    if 1/x3_(3) < 1/x3_(4)
        
        fitting_results.TripleExpo_fixedT0.P0 = x3_(1);
        fitting_results.TripleExpo_fixedT0.T0 = fixed_short_lifetime;        
        fitting_results.TripleExpo_fixedT0.P1 = x3_(2);
        fitting_results.TripleExpo_fixedT0.T1 = 1/x3_(3);
        fitting_results.TripleExpo_fixedT0.P2 = 1 - x3_(1) - x3_(2);
        fitting_results.TripleExpo_fixedT0.T2 = 1/x3_(4);
        fitting_results.TripleExpo_fixedT0.flag = exitflag;        
              
        %---------------
        % caluclation of R, the denominator in the
        % double_exp3_beta_norm_fixedT0
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple_(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
            R.(name_embryo).triple_(2,1) = sum( exp(-xdata.*x3_(3)) );
            R.(name_embryo).triple_(3,1) = sum( exp(-xdata.*x3_(4)) );
            clear xdata
        end
        fitting_results.TripleExpo_fixedT0.R_normalization = R;
        
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot0 = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm0 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA0 = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0( x3_,input.(name_embryo).data,R.(name_embryo).triple_(1,1),...
                R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
             size_population_tot0 = size_population_tot0 + double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),size_population2.(name_embryo).raw)),3)) );           
            size_population_tot1 = size_population_tot1 + double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],...
                input.(name_embryo).data,R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0(x3_,input.(name_embryo).data,...
                R.(name_embryo).triple_(1,1),R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw,...
                fixed_short_lifetime)),3)) );
            size_population.(name_embryo).pop0 = double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(2,1),size_population2.(name_embryo).raw)),3)) );            
            size_population.(name_embryo).pop2 = double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],...
                input.(name_embryo).data,R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop0 = size_population.(name_embryo).pop0 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop0 = size_population.(name_embryo).pop0 / input.(name_embryo).duration_phase;
                size_population_norm0 = size_population_norm0 + size_population_normalized.(name_embryo).pop0 ;                
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
             
                size_population_normalizedA.(name_embryo).pop0 = size_population_normalized.(name_embryo).pop0 *60 / input.(name_embryo).area;
                size_population_normA0 = size_population_normA0 + size_population_normalizedA.(name_embryo).pop0 ;                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total0 = double(size_population_tot0);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.TripleExpo_fixedT0.size_population = size_population;
        percent_population_tot0 =  size_population_tot0 /  size_population_tot *100;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total0 = double(percent_population_tot0);
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.TripleExpo_fixedT0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean0 = size_population_norm0 / nbEmbryo_givenCondition;
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean0 = size_population_normA0 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo_fixedT0.(name_embryo).triple_ = triple_exp2_beta_norm_fixedT0(x3_,input.(name_embryo).data,...
                R.(name_embryo).triple_(1,1),R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw,...
                fixed_short_lifetime);
            residues = (fitting_results.TripleExpo_fixedT0.(name_embryo).triple_(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo_fixedT0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
        
    elseif 1/x3_(3) > 1/x3_(4)
        
        fitting_results.TripleExpo_fixedT0.P0 = x3_(1);
        fitting_results.TripleExpo_fixedT0.T0 = fixed_short_lifetime;        
        fitting_results.TripleExpo_fixedT0.P1 = 1 - x3_(1) - x3_(2);
        fitting_results.TripleExpo_fixedT0.T1 = 1/x3_(4);
        fitting_results.TripleExpo_fixedT0.P2 = x3_(2);
        fitting_results.TripleExpo_fixedT0.T2 = 1/x3_(3);
        fitting_results.TripleExpo_fixedT0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple_(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
            R.(name_embryo).triple_(2,1) = sum( exp(-xdata.*x3_(4)) );
            R.(name_embryo).triple_(3,1) = sum( exp(-xdata.*x3_(3)) );
            clear xdata
        end
        fitting_results.TripleExpo_fixedT0.R_normalization = R;
 
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot0 = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm0 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA0 = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0([x3_(1) (1-x3_(1)-x3_(2)) x3_(4) x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
             size_population_tot0 = size_population_tot0 + double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),size_population2.(name_embryo).raw)),3)) );           
            size_population_tot1 = size_population_tot1 + double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],...
                input.(name_embryo).data,R.(name_embryo).triple_(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0([x3_(1) (1-x3_(1)-x3_(2)) x3_(4) x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime)),3)) );
            size_population.(name_embryo).pop0 = double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],...
                input.(name_embryo).data,R.(name_embryo).triple_(2,1),size_population2.(name_embryo).raw)),3)) );            
            size_population.(name_embryo).pop2 = double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop0 = size_population.(name_embryo).pop0 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop0 = size_population.(name_embryo).pop0 / input.(name_embryo).duration_phase;
                size_population_norm0 = size_population_norm0 + size_population_normalized.(name_embryo).pop0 ;                
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
             
                size_population_normalizedA.(name_embryo).pop0 = size_population_normalized.(name_embryo).pop0 *60 / input.(name_embryo).area;
                size_population_normA0 = size_population_normA0 + size_population_normalizedA.(name_embryo).pop0 ;                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total0 = double(size_population_tot0);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.TripleExpo_fixedT0.size_population = size_population;
        percent_population_tot0 =  size_population_tot0 /  size_population_tot *100;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total0 = double(percent_population_tot0);
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.TripleExpo_fixedT0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean0 = size_population_norm0 / nbEmbryo_givenCondition;
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean0 = size_population_normA0 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo_fixedT0.(name_embryo).triple_ = triple_exp2_beta_norm_fixedT0([x3_(1) (1-x3_(1)-x3_(2)) x3_(4) x3_(3)],...
                input.(name_embryo).data,R.(name_embryo).triple_(1,1),R.(name_embryo).triple_(2,1),R.(name_embryo).triple_(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime);
            residues = (fitting_results.TripleExpo_fixedT0.(name_embryo).triple_(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo_fixedT0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
     
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot0 size_population_tot1 size_population_tot2 
    clear residues
    clear size_population percentage_population percent_population_tot0 percent_population_tot1 percent_population_tot2
    clear R
    clear x3_
    clear size_population_norm0 size_population_norm1 size_population_norm2 size_population_normalized
    clear size_population_normA0 size_population_normA1 size_population_normA2 size_population_normalizedA size_population_normA_tot
    clear triple_
    
end


%% Triple expo model with fixed short lifetime and percent

if ~isempty(find(ismember(models,'TripleExpo_fixedT0P0'),1)) && exist('fixed_short_lifetime','var') == 1  && ~isempty(fixed_short_lifetime)...
        && exist('fixed_short_percent','var') == 1 && ~isempty(fixed_short_percent)

   [x3__, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'TripleExpo_fixedT0P0', @triple_exp2_beta_norm_fixedT0P0, nbEmbryo_givenCondition, ...
       size_population2, fixed_short_lifetime,fixed_short_percent ), [0.4 1.5 0.5], [], [], [], [],[0 0 0], [], [], Options);
    
    fitting_results.TripleExpo_fixedT0P0.parameters = x3__;
    fitting_results.TripleExpo_fixedT0P0.nparam = 3;
    
    if 1/x3__(2) < 1/x3__(3)
        
        fitting_results.TripleExpo_fixedT0P0.P0 = fixed_short_percent;
        fitting_results.TripleExpo_fixedT0P0.T0 = fixed_short_lifetime;        
        fitting_results.TripleExpo_fixedT0P0.P1 = x3__(1);
        fitting_results.TripleExpo_fixedT0P0.T1 = 1/x3__(2);
        fitting_results.TripleExpo_fixedT0P0.P2 = 1 - x3__(1) - fixed_short_percent;
        fitting_results.TripleExpo_fixedT0P0.T2 = 1/x3__(3);
        fitting_results.TripleExpo_fixedT0P0.flag = exitflag;        
              
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple__(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
            R.(name_embryo).triple__(2,1) = sum( exp(-xdata.*x3__(2)) );
            R.(name_embryo).triple__(3,1) = sum( exp(-xdata.*x3__(3)) );
            clear xdata
        end
        fitting_results.TripleExpo_fixedT0P0.R_normalization = R;
        
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot0 = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm0 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA0 = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0P0( x3__,input.(name_embryo).data,R.(name_embryo).triple__(1,1),...
                R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw,fixed_short_lifetime,...
                fixed_short_percent)),3)) );
             size_population_tot0 = size_population_tot0 + double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),size_population2.(name_embryo).raw)),3)) );           
            size_population_tot1 = size_population_tot1 + double( x3__(1)*floor(vpa(sum(simple_exp2_beta_norm([x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( (1-x3__(1)-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x3__(3)],...
                input.(name_embryo).data,R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0P0(x3__,input.(name_embryo).data,...
                R.(name_embryo).triple__(1,1),R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw,...
                fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population.(name_embryo).pop0 = double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x3__(1)*floor(vpa(sum(simple_exp2_beta_norm([x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(2,1),size_population2.(name_embryo).raw)),3)) );            
            size_population.(name_embryo).pop2 = double( (1-x3__(1)-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x3__(3)],...
                input.(name_embryo).data,R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop0 = size_population.(name_embryo).pop0 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop0 = size_population.(name_embryo).pop0 / input.(name_embryo).duration_phase;
                size_population_norm0 = size_population_norm0 + size_population_normalized.(name_embryo).pop0 ;                
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
             
                size_population_normalizedA.(name_embryo).pop0 = size_population_normalized.(name_embryo).pop0 *60 / input.(name_embryo).area;
                size_population_normA0 = size_population_normA0 + size_population_normalizedA.(name_embryo).pop0 ;                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total0 = double(size_population_tot0);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.TripleExpo_fixedT0P0.size_population = size_population;
        percent_population_tot0 =  size_population_tot0 /  size_population_tot *100;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total0 = double(percent_population_tot0);
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.TripleExpo_fixedT0P0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean0 = size_population_norm0 / nbEmbryo_givenCondition;
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0P0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean0 = size_population_normA0 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0P0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo_fixedT0P0.(name_embryo).triple__ = triple_exp2_beta_norm_fixedT0P0(x3__,input.(name_embryo).data,...
                R.(name_embryo).triple__(1,1),R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw,...
                fixed_short_lifetime,fixed_short_percent);
            residues = (fitting_results.TripleExpo_fixedT0P0.(name_embryo).triple__(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo_fixedT0P0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
        
        
    elseif 1/x3__(2) > 1/x3__(3)
        
        fitting_results.TripleExpo_fixedT0P0.P0 = fixed_short_percent;
        fitting_results.TripleExpo_fixedT0P0.T0 = fixed_short_lifetime;        
        fitting_results.TripleExpo_fixedT0P0.P1 = 1 - x3__(1) - fixed_short_percent;
        fitting_results.TripleExpo_fixedT0P0.T1 = 1/x3__(3);
        fitting_results.TripleExpo_fixedT0P0.P2 = x3__(1);
        fitting_results.TripleExpo_fixedT0P0.T2 = 1/x3__(2);
        fitting_results.TripleExpo_fixedT0P0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).triple__(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
            R.(name_embryo).triple__(2,1) = sum( exp(-xdata.*x3__(3)) );
            R.(name_embryo).triple__(3,1) = sum( exp(-xdata.*x3__(2)) );
            clear xdata
        end
        fitting_results.TripleExpo_fixedT0P0.R_normalization = R;
 
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot0 = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_norm0 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm_tot = 0;
        size_population_normA0 = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA_tot = 0;        
        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            
            size_population_tot = size_population_tot + ...
                double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0P0([fixed_short_percent (1-x3__(1)-fixed_short_percent) x3__(3) x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
             size_population_tot0 = size_population_tot0 + double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),size_population2.(name_embryo).raw)),3)) );           
            size_population_tot1 = size_population_tot1 + double( (1-x3__(1)-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x3__(3)],...
                input.(name_embryo).data,R.(name_embryo).triple__(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x3__(1)*floor(vpa(sum(simple_exp2_beta_norm([x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0P0([x3__(1) (1-x3__(1)-fixed_short_percent) x3__(3) x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent)),3)) );
            size_population.(name_embryo).pop0 = double( fixed_short_percent*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( (1-x3__(1)-fixed_short_percent)*floor(vpa(sum(simple_exp2_beta_norm([x3__(3)],...
                input.(name_embryo).data,R.(name_embryo).triple__(2,1),size_population2.(name_embryo).raw)),3)) );            
            size_population.(name_embryo).pop2 = double( x3__(1)*floor(vpa(sum(simple_exp2_beta_norm([x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(3,1),size_population2.(name_embryo).raw)),3)) );
            percentage_population.(name_embryo).pop0 = size_population.(name_embryo).pop0 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop0 = size_population.(name_embryo).pop0 / input.(name_embryo).duration_phase;
                size_population_norm0 = size_population_norm0 + size_population_normalized.(name_embryo).pop0 ;                
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
             
                size_population_normalizedA.(name_embryo).pop0 = size_population_normalized.(name_embryo).pop0 *60 / input.(name_embryo).area;
                size_population_normA0 = size_population_normA0 + size_population_normalizedA.(name_embryo).pop0 ;                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot= size_population_normA_tot + size_population_normalizedA.(name_embryo).total;                
            end
        end
        size_population.total = double(size_population_tot);
        size_population.total0 = double(size_population_tot0);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        fitting_results.TripleExpo_fixedT0P0.size_population = size_population;
        percent_population_tot0 =  size_population_tot0 /  size_population_tot *100;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percentage_population.total0 = double(percent_population_tot0);
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        fitting_results.TripleExpo_fixedT0P0.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean0 = size_population_norm0 / nbEmbryo_givenCondition;
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0P0.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean0 = size_population_normA0 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.TripleExpo_fixedT0P0.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.TripleExpo_fixedT0P0.(name_embryo).triple__ = triple_exp2_beta_norm_fixedT0P0([fixed_short_percent (1-x3__(1)-fixed_short_percent) x3__(3) x3__(2)],...
                input.(name_embryo).data,R.(name_embryo).triple__(1,1),R.(name_embryo).triple__(2,1),R.(name_embryo).triple__(3,1),...
                size_population2.(name_embryo).raw,fixed_short_lifetime,fixed_short_percent);
            residues = (fitting_results.TripleExpo_fixedT0P0.(name_embryo).triple__(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.TripleExpo_fixedT0P0.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
     
        
    end
    
    clear exitflag
    clear size_population_tot size_population_tot0 size_population_tot1 size_population_tot2 
    clear residues
    clear size_population percentage_population percent_population_tot0 percent_population_tot1 percent_population_tot2
    clear R
    clear x3__
    clear size_population_norm0 size_population_norm1 size_population_norm2 size_population_normalized
    clear size_population_normA0 size_population_normA1 size_population_normA2 size_population_normalizedA size_population_normA_tot
    clear triple__
    
end





% notice before residues are ponderated ones!!


%% fit with quadro expo as model

if ~isempty(find(ismember(models,'QuadroExpo'),1))
       
%    [x4, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'QuadroExpo', @quadro_exp2_beta_norm,nbEmbryo_givenCondition,size_population2), ...
%        [0.25 5 0.25 1.4 0.25 0.6 0.3], [1 0 1 0 1 0 0], [1], [], [],[0 1 0 0.3 0 0.2 0.1], [1 1000 1 2 1 1 0.3], [], Options);
    [x4, ~, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'QuadroExpo', @quadro_exp2_beta_norm,nbEmbryo_givenCondition,size_population2), ...
        [0.25 5 0.25 1.4 0.25 0.6 0.3], [1 0 1 0 1 0 0], [1], [], [],[0 0 0 0 0 0 0], [], [], Options);
    
    fitting_results.QuadroExpo.parameters = x4;
    fitting_results.QuadroExpo.nparam = 7;
    
    fitting_results.QuadroExpo.PPP1 = x4(1);
    fitting_results.QuadroExpo.TTT1 = 1/x4(2);
    fitting_results.QuadroExpo.PPP2 = x4(3);
    fitting_results.QuadroExpo.TTT2 = 1/x4(4);
    fitting_results.QuadroExpo.PPP3 = x4(5);
    fitting_results.QuadroExpo.TTT3 = 1/x4(6);
    fitting_results.QuadroExpo.PPP4 = 1-x4(1)-x4(3)-x4(5);
    fitting_results.QuadroExpo.TTT4 = 1/x4(7);    
    fitting_results.QuadroExpo.flag = exitflag;
    
    %if 1/x4(2) < 1/x4(4) < 1/x4(6) < 1/x4(7)
        
        %---------------
        % caluclation of R, the denominator in the quadro_exp2_beta_norm
        % function due to the normalization
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            xdata=input.(name_embryo).data(:,1);
            R.(name_embryo).quadro(1,1) = sum( exp(-xdata.*x4(2)) );
            R.(name_embryo).quadro(2,1) = sum( exp(-xdata.*x4(4)) );
            R.(name_embryo).quadro(3,1) = sum( exp(-xdata.*x4(6)) );
            R.(name_embryo).quadro(4,1) = sum( exp(-xdata.*x4(7)) );
            clear xdata
        end
        fitting_results.QuadroExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population_tot = 0;
        size_population_tot1 = 0;
        size_population_tot2 = 0;
        size_population_tot3 = 0;
        size_population_tot4 = 0;
        size_population_norm1 = 0;
        size_population_norm2 = 0;
        size_population_norm3 = 0;
        size_population_norm4 = 0;
        size_population_norm_tot = 0;
        size_population_normA1 = 0;
        size_population_normA2 = 0;
        size_population_normA3 = 0;
        size_population_normA4 = 0;
        size_population_normA_tot = 0;        
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            size_population_tot = size_population_tot + double( floor(vpa(sum(quadro_exp2_beta_norm(x4,input.(name_embryo).data,R.(name_embryo).quadro(1,1),...
                R.(name_embryo).quadro(2,1),R.(name_embryo).quadro(3,1),R.(name_embryo).quadro(4,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot1 = size_population_tot1 + double( x4(1)*floor(vpa(sum(simple_exp2_beta_norm([x4(2)],input.(name_embryo).data,...
                R.(name_embryo).quadro(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot2 = size_population_tot2 + double( x4(3)*floor(vpa(sum(simple_exp2_beta_norm([x4(4)],input.(name_embryo).data,...
                R.(name_embryo).quadro(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot3 = size_population_tot3 + double( (x4(5))*floor(vpa(sum(simple_exp2_beta_norm([x4(6)],input.(name_embryo).data,...
                R.(name_embryo).quadro(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population_tot4 = size_population_tot4 + double( (1-x4(1)-x4(3)-x4(5))*floor(vpa(sum(simple_exp2_beta_norm([x4(7)],input.(name_embryo).data,...
                R.(name_embryo).quadro(4,1),size_population2.(name_embryo).raw)),3)) );
            
            size_population.(name_embryo).total = double( floor(vpa(sum(quadro_exp2_beta_norm(x4,input.(name_embryo).data,...
                R.(name_embryo).quadro(1,1),R.(name_embryo).quadro(2,1),R.(name_embryo).quadro(3,1),R.(name_embryo).quadro(4,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop1 = double( x4(1)*floor(vpa(sum(simple_exp2_beta_norm([x4(2)],input.(name_embryo).data,...
                R.(name_embryo).quadro(1,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop2 = double( x4(3)*floor(vpa(sum(simple_exp2_beta_norm([x4(4)],input.(name_embryo).data,...
                R.(name_embryo).quadro(2,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop3 = double( (x4(5))*floor(vpa(sum(simple_exp2_beta_norm([x4(6)],input.(name_embryo).data,...
                R.(name_embryo).quadro(3,1),size_population2.(name_embryo).raw)),3)) );
            size_population.(name_embryo).pop4 = double( (1-x4(1)-x4(3)-x4(5))*floor(vpa(sum(simple_exp2_beta_norm([x4(7)],input.(name_embryo).data,...
                R.(name_embryo).quadro(4,1),size_population2.(name_embryo).raw)),3)) );            
            percentage_population.(name_embryo).pop1 = size_population.(name_embryo).pop1 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop2 = size_population.(name_embryo).pop2 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop3 = size_population.(name_embryo).pop3 / size_population.(name_embryo).total *100;
            percentage_population.(name_embryo).pop4 = size_population.(name_embryo).pop4 / size_population.(name_embryo).total *100;
            
            if simulation == 0
                size_population_normalized.(name_embryo).pop1 = size_population.(name_embryo).pop1 / input.(name_embryo).duration_phase;
                size_population_norm1 = size_population_norm1 + size_population_normalized.(name_embryo).pop1 ;
                size_population_normalized.(name_embryo).pop2 = size_population.(name_embryo).pop2 / input.(name_embryo).duration_phase;
                size_population_norm2 = size_population_norm2 + size_population_normalized.(name_embryo).pop2 ;
                size_population_normalized.(name_embryo).pop3 = size_population.(name_embryo).pop3 / input.(name_embryo).duration_phase;
                size_population_norm3 = size_population_norm3 + size_population_normalized.(name_embryo).pop3 ;
                size_population_normalized.(name_embryo).pop4 = size_population.(name_embryo).pop4 / input.(name_embryo).duration_phase;
                size_population_norm4 = size_population_norm4 + size_population_normalized.(name_embryo).pop4 ;                
                size_population_normalized.(name_embryo).total = size_population.(name_embryo).total / input.(name_embryo).duration_phase;
                size_population_norm_tot = size_population_norm_tot + size_population_normalized.(name_embryo).total;
                
                size_population_normalizedA.(name_embryo).pop1 = size_population_normalized.(name_embryo).pop1 *60 / input.(name_embryo).area;
                size_population_normA1 = size_population_normA1 + size_population_normalizedA.(name_embryo).pop1 ;
                size_population_normalizedA.(name_embryo).pop2 = size_population_normalized.(name_embryo).pop2 *60 / input.(name_embryo).area;
                size_population_normA2 = size_population_normA2 + size_population_normalizedA.(name_embryo).pop2 ;
                size_population_normalizedA.(name_embryo).pop3 = size_population_normalized.(name_embryo).pop3 *60 / input.(name_embryo).area;
                size_population_normA3 = size_population_normA3 + size_population_normalizedA.(name_embryo).pop3 ;
                 size_population_normalizedA.(name_embryo).pop4 = size_population_normalized.(name_embryo).pop4 *60 / input.(name_embryo).area;
                size_population_normA4 = size_population_normA4 + size_population_normalizedA.(name_embryo).pop4 ;               
                size_population_normalizedA.(name_embryo).total = size_population_normalized.(name_embryo).total *60 / input.(name_embryo).area;
                size_population_normA_tot = size_population_normA_tot + size_population_normalizedA.(name_embryo).total;               
            end
            
        end
        size_population.total = double(size_population_tot);
        size_population.total1 = double(size_population_tot1);
        size_population.total2 = double(size_population_tot2);
        size_population.total3 = double(size_population_tot3);
        size_population.total4 = double(size_population_tot4);
        fitting_results.QuadroExpo.size_population = size_population;
        percent_population_tot1 =  size_population_tot1 /  size_population_tot *100;
        percent_population_tot2 =  size_population_tot2 /  size_population_tot *100;
        percent_population_tot3 =  size_population_tot3 /  size_population_tot *100;
        percent_population_tot4 =  size_population_tot4 /  size_population_tot *100;
        percentage_population.total1 = double(percent_population_tot1);
        percentage_population.total2 = double(percent_population_tot2);
        percentage_population.total3 = double(percent_population_tot3);
        percentage_population.total4 = double(percent_population_tot4);
        fitting_results.QuadroExpo.percentage_population = percentage_population;
        if simulation == 0
            size_population_normalized.mean1 = size_population_norm1 / nbEmbryo_givenCondition;
            size_population_normalized.mean2 = size_population_norm2 / nbEmbryo_givenCondition;
            size_population_normalized.mean3 = size_population_norm3 / nbEmbryo_givenCondition;
            size_population_normalized.mean4 = size_population_norm4 / nbEmbryo_givenCondition;
            size_population_normalized.mean = size_population_norm_tot / nbEmbryo_givenCondition;
            fitting_results.QuadroExpo.size_population_normalized_time = size_population_normalized; % nbMTs / second
            
            size_population_normalizedA.mean1 = size_population_normA1 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean2 = size_population_normA2 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean3 = size_population_normA3 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean4 = size_population_normA4 / nbEmbryo_givenCondition;
            size_population_normalizedA.mean = size_population_normA_tot / nbEmbryo_givenCondition;
            fitting_results.QuadroExpo.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTs / min /um2            
        end
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        for iEmbryo = 1 : nbEmbryo_givenCondition
            name_embryo = ['embryo' num2str(iEmbryo)];
            fitting_results.QuadroExpo.(name_embryo).quadro = quadro_exp2_beta_norm(x4,input.(name_embryo).data,...
                R.(name_embryo).quadro(1,1),R.(name_embryo).quadro(2,1),R.(name_embryo).quadro(3,1),R.(name_embryo).quadro(4,1),size_population2.(name_embryo).raw);
            residues = (fitting_results.QuadroExpo.(name_embryo).quadro(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
            residues(~isfinite(residues))=NaN;
            fitting_results.QuadroExpo.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
        end
   
    clear exitflag
    clear size_population_tot size_population_tot1 size_population_tot2 size_population_tot3 size_population_tot4
    clear percent_population_tot1 percent_population_tot2 percent_population_tot3 percent_population_tot4
    clear residues
    clear size_population percentage_population
    clear x4
    clear size_population_norm1 size_population_norm2 size_population_norm3 size_population_norm4 size_population_normalized
    clear size_population_normA1 size_population_normA2 size_population_normA3 size_population_normA4 size_population_normalizedA
    
end


%% fit with drift and diffusion model

if ~isempty(find(ismember(models,'Drift_diffusion'),1))
    
    [x_d, fval_dd, exitflag]=fmincon(@(par) Loglikelihood2_total(par, input,'Drift_diffusion', @diffusion_drift_model_norm, nbEmbryo_givenCondition, size_population2), ...
        [1], [], [], [], [],[0], [],[],Options);
    
    fitting_results.Drift_diffusion.parameters = x_d;
    fitting_results.Drift_diffusion.nparam = 1;
    fitting_results.Drift_diffusion.T = 1/x_d(1);
    fitting_results.Drift_diffusion.flag = exitflag;
    
    %---------------
    % caluclation of R, the denominator in the diffusion_drift_model_norm
    % function due to the normalization
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        xdata=input.(name_embryo).data(:,1);
        R.(name_embryo).dd = sum( xdata.^(-3/2) .* exp(-xdata.*x_d(1)) );
        clear xdata
    end
    fitting_results.Drift_diffusion.R_normalization = R;
    %----------------
    % calculation fo the size of the population after fitting for each
    % embryo and the total one
    size_population_tot = 0;
    size_population_norm = 0;
    size_population_normA = 0;
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        size_population_tot = size_population_tot + double( floor(vpa(sum(diffusion_drift_model_norm(x_d,input.(name_embryo).data,R.(name_embryo).dd, ...
            size_population2.(name_embryo).raw)),3)) );
        size_population.(name_embryo) = double( floor(vpa(sum(diffusion_drift_model_norm(x_d,input.(name_embryo).data, R.(name_embryo).dd,...
            size_population2.(name_embryo).raw)),3)) );
        if simulation == 0
            size_population_normalized.(name_embryo) = size_population.(name_embryo) / input.(name_embryo).duration_phase;
            size_population_norm = size_population_norm + size_population_normalized.(name_embryo) ;
            size_population_normalizedA.(name_embryo) = size_population_normalized.(name_embryo)*60 / input.(name_embryo).area;
            size_population_normA = size_population_normA + size_population_normalizedA.(name_embryo) ;
        end
    end
    size_population.total = double(size_population_tot);
    fitting_results.Drift_diffusion.size_population = size_population;
    if simulation == 0
        size_population_normalized.mean = size_population_norm / nbEmbryo_givenCondition;
        fitting_results.Drift_diffusion.size_population_normalized_time = size_population_normalized; % nbMTs / second
        size_population_normalizedA.mean = size_population_normA / nbEmbryo_givenCondition;
        fitting_results.Drift_diffusion.size_population_normalized_timeAndArea = size_population_normalizedA; % nbMTS/min/um2
    end
    %-----------------
    % calculate the function value (mi) at each bins and the residues
    % between mi and ci
    for iEmbryo = 1 : nbEmbryo_givenCondition
        name_embryo = ['embryo' num2str(iEmbryo)];
        fitting_results.Drift_diffusion.(name_embryo).d_d = diffusion_drift_model_norm(x_d,input.(name_embryo).data,R.(name_embryo).dd,size_population2.(name_embryo).raw);
        residues = (fitting_results.Drift_diffusion.(name_embryo).d_d(:) - input.(name_embryo).data(:,2)) ./input.(name_embryo).data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.Drift_diffusion.(name_embryo).residuals = cat(2,input.(name_embryo).data(:,1),residues);
    end
    clear exitflag
    clear size_population
    clear residues
    clear R
    clear size_population_tot
    clear x_d
    clear size_population_normalized
    clear size_population_norm
    clear size_population_normalizedA
    clear size_population_normA
    
else
    fval_dd = NaN;
end


end


