function [ fitting_results ] = to_calculate_parameters_using_fmincon_mle_(data,models,timing_phase,blastomere,fitting_results,fixed_short_lifetime )

% fitting of a single embryo using fmincon function

% check for any inf in data
[ ri] = find(~isfinite(data(:,2)));
data(ri,:)= [];

size_population_raw = nansum(data(:,2));
fitting_results.(blastomere).(timing_phase).size_population_raw = size_population_raw;
fitting_results.(blastomere).(timing_phase).data = data;


%% choice of the options for the fitting

Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);

fitting_results.(blastomere).(timing_phase).options = Options;


%% fit with mono expo as models

if ~isempty(find(ismember(models,'MonoExpo'),1))
    
    [x, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'MonoExpo', @simple_exp2_beta_norm, size_population_raw), [1],...
        [], [], [], [],[0.05], [30],[],Options);
    
    fitting_results.(blastomere).(timing_phase).MonoExpo.parameters = x;
    fitting_results.(blastomere).(timing_phase).MonoExpo.nparam = 1;
    fitting_results.(blastomere).(timing_phase).MonoExpo.T = 1/x(1);
    fitting_results.(blastomere).(timing_phase).MonoExpo.flag = exitflag;
    
    %---------------
    % caluclation of R, the denominator in the simple_exp2_beta_norm
    % function due to the normalization
    xdata=data(:,1);
    R = sum( exp(-xdata.*x(1)) );
    clear xdata
    fitting_results.(blastomere).(timing_phase).MonoExpo.R_normalization = R;
    %----------------
    % calculation fo the size of the population after fitting for each
    % embryo and the total one
    size_population = double( floor(vpa(sum(simple_exp2_beta_norm(x,data,R,size_population_raw)),3)) );
    fitting_results.(blastomere).(timing_phase).MonoExpo.size_population = size_population;
    %-----------------
    % calculate the function value (mi) at each bins and the residues
    % between mi and ci
    fitting_results.(blastomere).(timing_phase).MonoExpo.simple = simple_exp2_beta_norm(x,data,R,size_population_raw);
    residues = (fitting_results.(blastomere).(timing_phase).MonoExpo.simple(:) - data(:,2)) ./data(:,2);
    residues(~isfinite(residues))=NaN;
    fitting_results.(blastomere).(timing_phase).MonoExpo.residuals = cat(2,data(:,1),residues);
    
    clear exitflag
    clear size_population
    clear residues
    clear R
    clear x
    
end

%% fit with double expo as model

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    
    [x2, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'DoubleExpo', @double_exp2_beta_norm, size_population_raw), ...
        [0.5 0.5 2], [], [], [], [],[0 0.1 0.1], [1 10  10], [], Options);
    
    fitting_results.(blastomere).(timing_phase).DoubleExpo.parameters = x2;
    fitting_results.(blastomere).(timing_phase).DoubleExpo.nparam = 3;
    
    if 1/x2(2) < 1/x2(3)
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P1 = x2(1);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T1 = 1/x2(2);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P2 = 1 - x2(1);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T2 = 1/x2(3);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.double(1,1) = sum( exp(-xdata.*x2(2)) );
        R.double(2,1) = sum( exp(-xdata.*x2(3)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).DoubleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        
        size_population.total = double( floor(vpa(sum(double_exp2_beta_norm(x2,data,R.double(1,1),R.double(2,1),size_population_raw)),3)) );
        size_population.pop1 = double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],data,R.double(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],data,R.double(2,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.percentage_population = percent_population;
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).DoubleExpo.double = double_exp2_beta_norm(x2,data,R.double(1,1),R.double(2,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).DoubleExpo.double(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x2(3) < 1/x2(2)
        
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P1 = 1-x2(1);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T1 = 1/x2(3);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.P2 = x2(1);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.T2 = 1/x2(2);
        fitting_results.(blastomere).(timing_phase).DoubleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.double(1,1) = sum( exp(-xdata.*x2(3)) );
        R.double(2,1) = sum( exp(-xdata.*x2(2)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).DoubleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(double_exp2_beta_norm([(1-x2(1)) x2(3) x2(2)],data,R.double(1,1),R.double(2,1),size_population_raw)),3)) );
        size_population.pop1 = double( (1-x2(1))*floor(vpa(sum(simple_exp2_beta_norm([x2(3)],data,R.double(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( x2(1)*floor(vpa(sum(simple_exp2_beta_norm([x2(2)],data,R.double(2,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).DoubleExpo.double = double_exp2_beta_norm([(1-x2(1)) x2(3) x2(2)],data,R.double(1,1),R.double(2,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).DoubleExpo.double(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).DoubleExpo.residuals = cat(2,data(:,1),residues);
        
    end
    clear exitflag
    clear residues
    clear size_population percentage_population
    clear R
    clear x2
    
end

%% fit with mono expo stretched as model

if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    
    [xS, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'MonoExpo_stretched', @simple_exp2_stretched_beta_norm, size_population_raw), ...
        [0.5 1 1],[], [], [], [],[0 0.1 0], [1 10 10], [],Options);
    
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.parameters = xS;
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.nparam = 3;
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.Ps = xS(1);
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.Ts = 1/xS(2);
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.power = xS(3);
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.flag = exitflag;
    
    %---------------
    % caluclation of R, the denominator in the simple_exp2_beta_norm
    % function due to the normalization
    xdata=data(:,1);
    R.monoS = sum( exp( -(xdata.*xS(2)).^(1/xS(3)) ) );
    clear xdata
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.R_normalization = R;
    %----------------
    % calculation fo the size of the population after fitting for each
    % embryo and the total one
    size_population = double( floor(vpa(sum(simple_exp2_stretched_beta_norm(xS,data, R.monoS,size_population_raw)),3)) );
       % *fitting_results.(blastomere).(timing_phase).size_population_raw),3)) );
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.size_population = size_population;
    %-----------------
    % calculate the function value (mi) at each bins and the residues
    % between mi and ci
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.simpleS = simple_exp2_stretched_beta_norm(xS,data,R.monoS,size_population_raw);
    residues = (fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.simpleS(:) - data(:,2)) ./data(:,2);
    residues(~isfinite(residues))=NaN;
    fitting_results.(blastomere).(timing_phase).MonoExpo_stretched.residuals = cat(2,data(:,1),residues);
    
    clear exitflag
    clear size_population
    clear residues
    clear R
    clear xS
    
end

%% fit with triple expo as model

if ~isempty(find(ismember(models,'TripleExpo'),1))
    
    [x3, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'TripleExpo', @triple_exp2_beta_norm, size_population_raw),...
        [0.5 0.5 0.5 2 5 ], [], [], [], [],[0 0.1 0 0.1 0.1], [1 100 1 100 100], [],Options);
    
    fitting_results.(blastomere).(timing_phase).TripleExpo.parameters = x3;
    fitting_results.(blastomere).(timing_phase).TripleExpo.nparam = 5;
    
    if 1/x3(2) < 1/x3(4) < 1/x3(5)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(2)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(4)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(5)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(1,1),R.triple(2,1),R.triple(3,1),size_population_raw)),3)) );
        size_population.pop1 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(1,1),size_population_raw)),3)) );    
        size_population.pop2 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(1,1),R.triple(2,1),R.triple(3,1), size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3(2)  < 1/x3(5) < 1/x3(4)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(2)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(5)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(4)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(1,1),R.triple(3,1),R.triple(2,1),size_population_raw)),3)) );
        size_population.pop1 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(1,1),R.triple(3,1),R.triple(2,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3(4)  < 1/x3(2) < 1/x3(5)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(4)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(2)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(5)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(2,1),R.triple(1,1),R.triple(3,1),size_population_raw)),3)) );
        size_population.pop1 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(2,1),R.triple(1,1),R.triple(3,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3(4)  < 1/x3(5) < 1/x3(2)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(4)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(5)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(2)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(2,1),R.triple(3,1),R.triple(1,1),size_population_raw)),3)) );
        size_population.pop1 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(2,1),R.triple(3,1),R.triple(1,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3(5)  < 1/x3(2) < 1/x3(4)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(5)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(2)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(4)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(3,1),R.triple(1,1),R.triple(2,1),size_population_raw)),3)) );
        size_population.pop1 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(3,1),R.triple(1,1),R.triple(2,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3(5)  < 1/x3(4) < 1/x3(2)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP1 = 1 - x3(1) - x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT1 = 1/x3(5);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP2 = x3(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT2 = 1/x3(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo.PP3 = x3(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo.TT3 = 1/x3(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the triple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple(1,1) = sum( exp(-xdata.*x3(6)) );
        R.triple(2,1) = sum( exp(-xdata.*x3(4)) );
        R.triple(3,1) = sum( exp(-xdata.*x3(2)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm(x3,data,R.triple(3,1),R.triple(2,1),R.triple(1,1),size_population_raw)),3)) );
        size_population.pop1 = double( (1-x3(1)-x3(3))*floor(vpa(sum(simple_exp2_beta_norm([x3(5)],data,R.triple(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( x3(3)*floor(vpa(sum(simple_exp2_beta_norm([x3(4)],data,R.triple(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( x3(1)*floor(vpa(sum(simple_exp2_beta_norm([x3(2)],data,R.triple(3,1),size_population_raw)),3)) );
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo.triple = triple_exp2_beta_norm(x3,data,R.triple(3,1),R.triple(2,1),R.triple(1,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo.triple(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo.residuals = cat(2,data(:,1),residues);
        
    end
    clear exitflag
    clear residues
    clear size_population percentage_population
    clear x3
    
end


%% fit with triple expo with fixed short lifetime as model

if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1)) && exist('fixed_short_lifetime','var') == 1
    
    [x3_, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'TripleExpo_fixedT0', @triple_exp2_beta_norm_fixedT0, size_population_raw, fixed_short_lifetime), ...
        [0.4 0.4 1.5 0.5], [], [], [], [],[0 0 0.1 0.1], [1 1 5 5], [], Options);
    
    fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.parameters = x3_;
    fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.nparam = 3;
    
    if 1/x3_(3) < 1/x3_(4)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P0 = x3_(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T0 = fixed_short_lifetime;        
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P1 = x3_(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T1 = 1/x3_(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P2 = 1 - x3_(1) - x3_(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T2 = 1/x3_(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple_(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
        R.triple_(2,1) = sum( exp(-xdata.*x3_(3)) );
        R.triple_(3,1) = sum( exp(-xdata.*x3_(4)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0(x3_,data,R.triple_(1,1),R.triple_(2,1),R.triple_(3,1),...
            size_population_raw, fixed_short_lifetime)),3)) );
        size_population.pop0 = double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],data,R.triple_(1,1),size_population_raw)),3)) );
        size_population.pop1 = double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],data,R.triple_(2,1),size_population_raw)),3)) );
        size_population.pop2 = double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],data,R.triple_(3,1),size_population_raw)),3)) );
        percent_population.pop0 = size_population.pop0 / size_population.total *100;
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.percentage_population = percent_population;
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.triple_ = triple_exp2_beta_norm_fixedT0(x3_,data,R.triple_(1,1),R.triple_(2,1),...
            R.triple_(3,1),size_population_raw, fixed_short_lifetime);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.residuals = cat(2,data(:,1),residues);
        
    elseif 1/x3_(3) > 1/x3_(4)
        
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P0 = x3_(1);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T0 = fixed_short_lifetime;        
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P1 = 1 - x3_(1) - x3_(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T1 = 1/x3_(4);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P2 = x3_(2);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T2 = 1/x3_(3);
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the simple_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.triple_(1,1) = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
        R.triple_(2,1) = sum( exp(-xdata.*x3_(4)) );
        R.triple_(3,1) = sum( exp(-xdata.*x3_(3)) );
        clear xdata
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        size_population.total = double( floor(vpa(sum(triple_exp2_beta_norm_fixedT0([x3_(1) (1-x3_(1)-x3_(2)) x3_(4) x3_(3)],data,...
            R.triple_(1,1),R.triple_(2,1),R.triple_(3,1),size_population_raw, fixed_short_lifetime)),3)) );
        size_population.pop0 = double( x3_(1)*floor(vpa(sum(simple_exp2_beta_norm([1/fixed_short_lifetime],data,R.triple_(1,1),size_population_raw)),3)) );
        size_population.pop1 = double( (1-x3_(1)-x3_(2))*floor(vpa(sum(simple_exp2_beta_norm([x3_(4)],data,R.triple_(2,1),size_population_raw)),3)) );
        size_population.pop2 = double( x3_(2)*floor(vpa(sum(simple_exp2_beta_norm([x3_(3)],data,R.triple_(3,1),size_population_raw)),3)) );
        percent_population.pop0 = size_population.pop0 / size_population.total *100;
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.percentage_population = percent_population;
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.triple_ = triple_exp2_beta_norm_fixedT0([(x3_(1)) (1-x3_(1)-x3_(2)) x3_(4)  x3_(3)],...
            data,R.triple_(1,1),R.triple_(2,1),R.triple_(3,1),size_population_raw, fixed_short_lifetime);
        residues = (fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.residuals = cat(2,data(:,1),residues);
        
    end
    clear exitflag
    clear residues
    clear size_population percentage_population
    clear R
    clear x3_
    
end

%% fit with quadro exponential model

if ~isempty(find(ismember(models,'QuadroExpo'),1))
    
    [x4, ~, exitflag]=fmincon(@(par) Loglikelihood2_norm(par, data,'QuadroExpo', @quadro_exp2_beta_norm, size_population_raw), ...
        [0.25 0.3 0.25 0.6 0.25 1 2], [], [], [], [],[0 0.05 0 0.05 0 0.05 0.05], [1 10 1 10 1 10 10], [], Options);
    
    fitting_results.(blastomere).(timing_phase).QuadroExpo.parameters = x4;
    fitting_results.(blastomere).(timing_phase).QuadroExpo.nparam = 7;
    
    %if 1/x4(2) < 1/x4(4) < 1/x4(6) < 1/x4(7)
        
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP1 = x4(1);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT1 = 1/x4(2);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP2 = x4(3);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT2 = 1/x4(4);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP3 = x4(5);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT3 = 1/x4(6);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP4 = 1 - x4(1) - x4(3) - x4(5);
        fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT4 = 1/x4(7);        
        fitting_results.(blastomere).(timing_phase).QuadroExpo.flag = exitflag;
        
        %---------------
        % caluclation of R, the denominator in the double_exp2_beta_norm
        % function due to the normalization
        xdata=data(:,1);
        R.quadro(1,1) = sum( exp(-xdata.*x4(2)) );
        R.quadro(2,1) = sum( exp(-xdata.*x4(4)) );
        R.quadro(3,1) = sum( exp(-xdata.*x4(6)) );
        R.quadro(4,1) = sum( exp(-xdata.*x4(7)) );        
        clear xdata
        fitting_results.(blastomere).(timing_phase).QuadroExpo.R_normalization = R;
        %----------------
        % calculation fo the size of the population after fitting for each
        % embryo and the total one
        
        size_population.total = double( floor(vpa(sum(quadro_exp2_beta_norm(x4,data,R.quadro(1,1),R.quadro(2,1),R.quadro(3,1),R.quadro(4,1),...
            size_population_raw)),3)) );
        size_population.pop1 = double( x4(1)*floor(vpa(sum(simple_exp2_beta_norm([x4(2)],data,R.quadro(1,1),size_population_raw)),3)) );
        size_population.pop2 = double( x4(3)*floor(vpa(sum(simple_exp2_beta_norm([x4(4)],data,R.quadro(2,1),size_population_raw)),3)) );
        size_population.pop3 = double( x4(5)*floor(vpa(sum(simple_exp2_beta_norm([x4(6)],data,R.quadro(1,1),size_population_raw)),3)) );
        size_population.pop4 = double( (1-x4(1)-x4(3)-x4(5))*floor(vpa(sum(simple_exp2_beta_norm([x4(7)],data,R.quadro(2,1),size_population_raw)),3)) );      
        percent_population.pop1 = size_population.pop1 / size_population.total *100;
        percent_population.pop2 = size_population.pop2 / size_population.total *100;
        percent_population.pop3 = size_population.pop3 / size_population.total *100;
        percent_population.pop4 = size_population.pop4 / size_population.total *100;        
        fitting_results.(blastomere).(timing_phase).QuadroExpo.size_population = size_population;
        fitting_results.(blastomere).(timing_phase).QuadroExpo.percentage_population = percent_population;
        
        %-----------------
        % calculate the function value (mi) at each bins and the residues
        % between mi and ci
        fitting_results.(blastomere).(timing_phase).QuadroExpo.quadro = quadro_exp2_beta_norm(x4,data,R.quadro(1,1),R.quadro(2,1),...
            R.quadro(3,1),R.quadro(4,1),size_population_raw);
        residues = (fitting_results.(blastomere).(timing_phase).QuadroExpo.quadro(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.(blastomere).(timing_phase).QuadroExpo.residuals = cat(2,data(:,1),residues);
        
        
    %end
    
    clear x4
    clear size_population
    clear percent_population
    clear residues 
    clear exitflag
    
end

clear data_norm size_population_raw


end

