function [fitting_results,size_population_raw] = tracks_duration_fitting_poisson...
    (bincounts,binranges,models,save_stem,name1,name2,name_choiceConditions,tracks_duration_histo)


%% fitting

data(:,1) = binranges;
data(:,2) = bincounts;

 ordinate = data(1,2);
 power = floor(log(abs(ordinate))./log(10));
 
 nrec = 2;

%Options=optimset('TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',100000, 'MaxIter' , 100000);
Options=optimset('TolFun',1e-15,'MaxFunEvals',10000);
fitting_results.options = Options;

% NB : Dividing par(2) by 1000 to get comparable numbers and better
% optimization
if ~isempty(find(ismember(models,'MonoExpo'),1))
    [x,resnorm,residual,exitflag] = lsqnonlin(@(par) simple_exp_opt(par,data,power), [1*10^(power) 1*10^(power)/2], [0 0 ],...
        [1*10^(power+2) 1*10^(power+2)], Options);
    fitting_results.MonoExpo.parameters = x;
    fitting_results.MonoExpo.A = x(1);
    fitting_results.MonoExpo.T = 1/ (x(2)*1*10^(-power));
    fitting_results.MonoExpo.residuals_weighted = residual;
    fitting_results.MonoExpo.flag = exitflag;
    fitting_results.MonoExpo.resnorm = resnorm; 
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    [x2,resnorm,residual,exitflag] = lsqnonlin(@(par) double_exp_opt(par,data,power), [1*10^(power)/2 1*10^(power) 1*10^(power)/2 1*10^(power)],...
        [0 0 0 0], [1*10^(power+2) 1*10^(power+2) 1*10^(power+2) 1*10^(power+2)], Options);
    fitting_results.DoubleExpo.parameters = x2;
    fitting_results.DoubleExpo.B1 = x2(1);
    fitting_results.DoubleExpo.T1 = 1/ ( x2(2)*1*10^(-power));
    fitting_results.DoubleExpo.B2 = x2(3);
    fitting_results.DoubleExpo.T2 = 1/ (x2(4)*1*10^(-power));
    fitting_results.DoubleExpo.residuals_weighted = residual;
    fitting_results.DoubleExpo.flag = exitflag;
    fitting_results.DoubleExpo.resnorm = resnorm; 
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    [xN1,resnorm,residual,exitflag] = lsqnonlin(@(par) simple_exp_optN1(par,data,power), [1*10^(power) 1*10^(power)/2], ...
        [0 0 ], [1*10^(power+2) 1*10^(power+2)], Options);
    fitting_results.MonoExpo_Neyman1.parameters = xN1;
    fitting_results.MonoExpo_Neyman1.A = xN1(1);
    fitting_results.MonoExpo_Neyman1.T = 1/ (xN1(2)*1*10^(-power));
    fitting_results.MonoExpo_Neyman1.residuals_weighted = residual;
    fitting_results.MonoExpo_Neyman1.flag = exitflag;
    fitting_results.MonoExpo_Neyman1.resnorm = resnorm; 
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    [xN2,resnorm,residual,exitflag] = lsqnonlin(@(par) simple_exp_optN2(par,data,power), [1*10^(power) 1*10^(power)/2], ...
        [0 0 ], [1*10^(power+2) 1*10^(power+2)], Options);
    fitting_results.MonoExpo_Neyman2.parameters = xN2;
    fitting_results.MonoExpo_Neyman2.A = xN2(1);
    fitting_results.MonoExpo_Neyman2.T = 1/ (xN2(2)*1*10^(-power));
    fitting_results.MonoExpo_Neyman2.residuals_weighted = residual;
    fitting_results.MonoExpo_Neyman2.flag = exitflag;
    fitting_results.MonoExpo_Neyman2.resnorm = resnorm; 
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    [xP,resnorm,residual,exitflag] = lsqnonlin(@(par) simple_exp_optP(par,data,power), [1*10^(power) 1*10^(power)/2], ...
        [0 0 ], [1*10^(power+2) 1*10^(power+2)], Options);
    fitting_results.MonoExpo_Pearson.parameters = xP;
    fitting_results.MonoExpo_Pearson.A = xP(1);
    fitting_results.MonoExpo_Pearson.T = 1/ (xP(2)*1*10^(-power));
    fitting_results.MonoExpo_Pearson.residuals_weighted = residual;
    fitting_results.MonoExpo_Pearson.flag = exitflag;
    fitting_results.MonoExpo_Pearson.resnorm = resnorm; 
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    [xFW,resnorm,residual,exitflag] = opti_simple_FW(data, [1*10^(power) 1*10^(power)/2], [0 0 ], ...
        [1*10^(power+3) 1*10^(power+3)], power, nrec, Options);
    fitting_results.MonoExpo_fitted_weighting.parameters = xFW;
    fitting_results.MonoExpo_fitted_weighting.A = xFW(1);
    fitting_results.MonoExpo_fitted_weighting.T = 1/( xFW(2)*1*10^(-power));
    fitting_results.MonoExpo_fitted_weighting.residuals_weighted = residual;
    fitting_results.MonoExpo_fitted_weighting.flag = exitflag;
    fitting_results.MonoExpo_fitted_weighting.resnorm = resnorm; 
    fitting_results.MonoExpo_fitted_weighting.nb_recurrences = nrec;
    clear resnorm
    clear residual
    clear exitflag
end

if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
     [x2FW,resnorm,residual,exitflag] = opti_double_FW(data,[1*10^(power)/2 1*10^(power) 1*10^(power)/2 1*10^(power)], ...
         [0 0 0 0], [1*10^(power+3) 1*10^(power+3) 1*10^(power+3) 1*10^(power+3)],power, nrec, Options);
    fitting_results.DoubleExpo_fitted_weighting.parameters = x2FW;
    fitting_results.DoubleExpo_fitted_weighting.B1 = x2FW(1);
    fitting_results.DoubleExpo_fitted_weighting.T1 = 1/ ( x2FW(2)*1*10^(-power));
    fitting_results.DoubleExpo_fitted_weighting.B2 = x2FW(3);
    fitting_results.DoubleExpo_fitted_weighting.T2 = 1/ (x2FW(4)*1*10^(-power));
    fitting_results.DoubleExpo_fitted_weighting.residuals_weighted = residual;
    fitting_results.DoubleExpo_fitted_weighting.flag = exitflag;
    fitting_results.DoubleExpo_fitted_weighting.resnorm = resnorm; 
    fitting_results.DoubleExpo_fitted_weighting.nb_recurrences = nrec;
    clear resnorm
    clear residual
    clear exitflag
end

% notice before residues are ponderated ones!!

%% Checking the population size

% why doing this checking? to se that not so much variation in total nb of
% counts

size_population_raw = sum(data(:,2));

if ~isempty(find(ismember(models,'MonoExpo'),1))
    fitting_results.MonoExpo.size_population=floor(vpa(sum(simple_exp(x,data,power)),3));
else
    fitting_results.MonoExpo.size_population = NaN;
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    fitting_results.MonoExpo_Neyman1.size_population=floor(vpa(sum(simple_exp(xN1,data,power)),3));
else
    fitting_results.MonoExpo_Neyman1.size_population = NaN;
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    fitting_results.MonoExpo_Neyman2.size_population=floor(vpa(sum(simple_exp(xN2,data,power)),3));
else
    fitting_results.MonoExpo_Neyman2.size_population = NaN;
end
if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    fitting_results.MonoExpo_Pearson.size_population=floor(vpa(sum(simple_exp(xP,data,power)),3));
else
    fitting_results.MonoExpo_Pearson.size_population = NaN;
end
if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    fitting_results.MonoExpo_fitted_weighting.size_population=floor(vpa(sum(simple_exp(xFW,data,power)),3)); 
else
    fitting_results.MonoExpo_fitted_weighting.size_population = NaN;
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    fitting_results.DoubleExpo.size_population=floor(vpa(sum(double_exp(x2,data,power)),3));
else
    fitting_results.DoubleExpo.size_population = NaN;
end
if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    fitting_results.DoubleExpo_fitted_weighting.size_population=floor(vpa(sum(double_exp(x2FW,data,power)),3));
else
    fitting_results.DoubleExpo_fitted_weighting.size_population = NaN;
end



%% Calculating the goodness-of-fit

weight_sum = 0; %normalization factor for model probabilities
AIC = [];

if ~isempty(find(ismember(models,'MonoExpo'),1))
    simple=simple_exp(x,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    simpleN1=simple_exp(xN1,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    simpleN2=simple_exp(xN2,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    simpleP=simple_exp(xP,data,power);    
end
if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    simpleFW=simple_exp(xFW,data,power);     
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    double=double_exp(x2,data,power);
end
if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    doubleFW=double_exp(x2FW,data,power);
end


%---------------------

if ~isempty(find(ismember(models,'MonoExpo'),1))
    MLE=0;
        for i=1:size(binranges)
            if data(i,2) > 0
            MLE= MLE + (data(i,2)-simple(i))-data(i,2)*log(simple(i)/data(i,2));
            else
                MLE = MLE + (data(i,2)-simple(i)) ;
            end
        end
    MLE=2*MLE;
    AIC_E = MLE +2*2; 
    AIC = [ [AIC] AIC_E];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    MLEN1=0;
        for i=1:size(binranges)
            if data(i,2) > 0
            MLEN1= MLEN1 + (data(i,2)-simpleN1(i))-data(i,2)*log(simpleN1(i)/data(i,2));
            else
                MLEN1 = MLEN1 + (data(i,2)-simpleN1(i)) ;
            end
        end
    AIC_EN1=2*MLEN1 +2*2;
    AIC = [ [AIC] AIC_EN1];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    MLEN2=0;
    for i=1:size(binranges)
        if data(i,2) > 0
        MLEN2= MLEN2 + (data(i,2)-simpleN2(i))-data(i,2)*log(simpleN2(i)/data(i,2));
        else
            MLEN2 = MLEN2 + (data(i,2)-simpleN2(i)) ;
        end
    end
    AIC_EN2=2*MLEN2 +2*2;
    AIC = [ [AIC] AIC_EN2];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    MLEP=0;
    for i=1:size(binranges)
        if data(i,2) > 0
        MLEP= MLEP + (data(i,2)-simpleP(i))-data(i,2)*log(simpleP(i)/data(i,2));
        else
            MLEP = MLEP + (data(i,2)-simpleP(i)) ;
        end
    end
    AIC_EP=2*MLEP +2*2;
    AIC = [ [AIC] AIC_EP];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    MLEFW=0;
    for i=1:size(binranges)
        if data(i,2) > 0
        MLEFW= MLEFW + (data(i,2)-simpleFW(i))-data(i,2)*log(simpleFW(i)/data(i,2));
        else
            MLEFW = MLEFW + (data(i,2)-simpleFW(i)) ;
        end
    end
    AIC_EFW=2*MLEFW +2*2;
    AIC = [ [AIC] AIC_EFW];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    MLE2=0;
        for i=1:size(binranges)
            if data(i,2) > 0
            MLE2= MLE2 + (data(i,2)-double(i))-data(i,2)*log(double(i)/data(i,2));
            else
                MLE2 = MLE2 + (data(i,2)-double(i)) ;
            end
        end
    AIC_E2=2*MLE2 +2*4;
    AIC = [ [AIC] AIC_E2];
else
    AIC = [ [AIC] NaN];
end

if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    MLE2FW=0;
        for i=1:size(binranges)
            if data(i,2) > 0
            MLE2FW= MLE2FW + (data(i,2)-doubleFW(i))-data(i,2)*log(doubleFW(i)/data(i,2));
            else
                MLE2FW = MLE2FW + (data(i,2)-doubleFW(i)) ;
            end
        end
    AIC_E2FW=2*MLE2FW +2*4;
    AIC = [ [AIC] AIC_E2FW];
else
    AIC = [ [AIC] NaN];
end

[AIC_min,Index_min]=min(AIC);

if Index_min == 1
    fitting_results.best_model = 'MonoExpo';
elseif Index_min ==2
    fitting_results.best_model = 'MonoExpo_Neyman1';
elseif Index_min ==3
    fitting_results.best_model = 'MonoExpo_Neyman2';
elseif Index_min ==4    
    fitting_results.best_model = 'MonoExpo_Pearson';
elseif Index_min ==5
    fitting_results.best_model = 'MonoExpo_fitted_weighting';
elseif Index_min ==6
    fitting_results.best_model = 'DoubleExpo';
elseif Index_min ==7    
    fitting_results.best_model = 'DoubleExpo_fitted_weighting';
end
    
%----------------------------
% why doing this calculation? to see difference when no ponderation
% it is the sum of residues squared

if ~isempty(find(ismember(models,'MonoExpo'),1)) 
    fitting_results.MonoExpo.sum_residues = sum(simple_exp_opt(x,data,power).^2); 
    fitting_results.MonoExpo.residuals = simple_exp_opt(x,data,power);
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    fitting_results.DoubleExpo.sum_residues = sum(double_exp_opt(x2,data,power).^2); 
    fitting_results.DoubleExpo.residuals = double_exp_opt(x2,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    fitting_results.MonoExpo_Neyman1.sum_residues = sum(simple_exp_opt(xN1,data,power).^2); 
    fitting_results.MonoExpo_Neyman1.residuals = simple_exp_opt(xN1,data,power); 
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    fitting_results.MonoExpo_Neyman2.sum_residues = sum(simple_exp_opt(xN2,data,power).^2);
    fitting_results.MonoExpo_Neyman2.residuals = simple_exp_opt(xN2,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    fitting_results.MonoExpo_Pearson.sum_residues = sum(simple_exp_opt(xP,data,power).^2);
    fitting_results.MonoExpo_Pearson.residuals = simple_exp_opt(xP,data,power);
end
if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    fitting_results.MonoExpo_fitted_weighting.sum_residues= sum(simple_exp_opt(xFW,data,power).^2);
    fitting_results.MonoExpo_fitted_weighting.residuals = simple_exp_opt(xFW,data,power);
end
if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    fitting_results.DoubleExpo_fitted_weighting.sum_residues = sum(double_exp_opt(x2FW,data,power).^2);
    fitting_results.DoubleExpo_fitted_weighting.residuals = double_exp_opt(x2FW,data,power);
end

%----------------------------
% Calculate the normalized probablity of each model 

if ~isempty(find(ismember(models,'MonoExpo'),1)) 
    weight_E = exp(-0.5*(AIC_E-AIC_min)); 
    weight_sum = weight_sum + weight_E;
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    weight_E2 = exp(-0.5*(AIC_E2-AIC_min)); 
    weight_sum = weight_sum + weight_E2;
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    weight_EN1 = exp(-0.5*(AIC_EN1-AIC_min)); 
    weight_sum = weight_sum + weight_EN1;
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    weight_EN2 = exp(-0.5*(AIC_EN2-AIC_min)); 
    weight_sum = weight_sum + weight_EN2;
end
if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    weight_EP = exp(-0.5*(AIC_EP-AIC_min)); 
    weight_sum = weight_sum + weight_EP;
end
if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    weight_EFW = exp(-0.5*(AIC_EFW-AIC_min)); 
    weight_sum = weight_sum + weight_EFW;
end
if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    weight_E2FW = exp(-0.5*(AIC_E2FW-AIC_min)); 
    weight_sum = weight_sum + weight_E2FW;
end

%-------------------------------
% Calculate the normalized probablity of each model 

if ~isempty(find(ismember(models,'MonoExpo'),1)) 
    fitting_results.MonoExpo.AIC = AIC_E; 
    fitting_results.MonoExpo.PrM = weight_E/weight_sum; 
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    fitting_results.DoubleExpo.AIC = AIC_E2; 
    fitting_results.DoubleExpo.PrM = weight_E2/weight_sum;
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    fitting_results.MonoExpo_Neyman1.AIC = AIC_EN1; 
    fitting_results.MonoExpo_Neyman1.PrM = weight_EN1/weight_sum; 
end
if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    fitting_results.MonoExpo_Neyman2.AIC = AIC_EN2;
    fitting_results.MonoExpo_Neyman2.PrM = weight_EN2/weight_sum;
end
if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    fitting_results.MonoExpo_Pearson.AIC = AIC_EP;
    fitting_results.MonoExpo_Pearson.PrM = weight_EP/weight_sum;
end
if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    fitting_results.MonoExpo_fitted_weighting.AIC = AIC_EFW;
    fitting_results.MonoExpo_fitted_weighting.PrM = weight_EFW/weight_sum;
end
if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    fitting_results.DoubleExpo_fitted_weighting.AIC = AIC_E2FW;
    fitting_results.DoubleExpo_fitted_weighting.PrM = weight_E2FW/weight_sum;
end

% %% 1 possibility to evaluate the errors on the fitting parameters
% 
% % P(A,T/Ci) = PI mi ^ ci exp(-mi) /ci!;
% % -2ln(L) = 2 sum_i (mi-ci*ln(mi)) +cste
% 
% cste = 0;
% for i=1:length(data)
%     for j=1:data(i,2)
%         cste = cste + log(j);
%     end
% end
% 
% size_population_lik = [];
% raw_size = sum(data(:,2));
% size_population_lik = [ [size_population_lik ] raw_size];
% 
% 
% if ~isempty(find(ismember(models,'MonoExpo'),1)) || ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1)) || ... 
%         ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1)) || ~isempty(find(ismember(models,'MonoExpo_Pearson'),1)) || ...
%         ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
%     
%     nA = 500;
%     nT= 500 ;
%     Amin = 0;
%     Amax = 5*10^(power);
%     Tmin = 0;
%     Tmax = 10000;
%     A=linspace(Amin,Amax,nA);
%     T=linspace(Tmin, Tmax, nT);
%     llik=zeros(nT,nA);
%     for i=1:nA
%         Ac=A(i);
%         for j=1:nT
%             Tc=T(j);
%             mod=simple_exp([Ac Tc],data,power);
%             llik(i,j)= sum (mod - data(:,2).*log(mod));
%         end
%     end
% 
%     figure;
%     contourf(A,T,llik,100);
%     xlabel('A')
%     ylabel('T (sec-1)')
%     title(['maximum likelihood : ' name1 ' and ' name2 ]);
%     namePlot = strcat('maximum_likelihood_MonoExpo-', name1 , '-', name2, '.fig');
%     saveas(gcf,[save_stem namePlot]);
% 
%     [i_ j_]=find(llik == min(min(llik)));
%     xlik=[A(i_) T(j_)];
%     fitting_results.MonoExpo.A_lik = xlik(1);
%     fitting_results.MonoExpo.T_lik = 1/ (xlik(2)/1*10^(-power));
%     size_ME_lik=floor(vpa(sum(simple_exp(xlik,data,power)),3));
%     size_population_lik = [ [size_population_lik ] size_ME_lik];   
% else
%     size_population_lik = [ [size_population_lik ] NaN];
% end
% 
% fitting_results.size_population = size_population;
% fitting_results.size_population_lik = size_population_lik;

%% 2nd possibility to evaluate errors : do boostrap on all the tracks (possible to mix different embryos??)


[ errors,errors2 ] = to_get_errors_using_bootstrap( tracks_duration_histo,models,save_stem,name1,name2,name_choiceConditions );

if ~isempty(find(ismember(models,'MonoExpo'),1))
    fitting_results.MonoExpo.A_se = errors2.MonoExpo.A_se;
    fitting_results.MonoExpo.T_se = errors2.MonoExpo.T_se;
else
    fitting_results.MonoExpo.A_se = NaN;
    fitting_results.MonoExpo.T_se = NaN;
end

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    fitting_results.DoubleExpo.B1_se = errors2.DoubleExpo.B1_se;
    fitting_results.DoubleExpo.T1_se = errors2.DoubleExpo.T1_se;
    fitting_results.DoubleExpo.B2_se = errors2.DoubleExpo.B2_se;
    fitting_results.DoubleExpo.T2_se = errors2.DoubleExpo.T2_se;
else
    fitting_results.DoubleExpo.B1_se = NaN;
    fitting_results.DoubleExpo.T1_se = NaN;
    fitting_results.DoubleExpo.B2_se = NaN;
    fitting_results.DoubleExpo.T2_se = NaN;    
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman1'),1))
    fitting_results.MonoExpo_Neyman1.A_se = errors2.MonoExpo_Neyman1.A_se;
    fitting_results.MonoExpo_Neyman1.T_se = errors2.MonoExpo_Neyman1.T_se;
else
    fitting_results.MonoExpo_Neyman1.A_se = NaN;
    fitting_results.MonoExpo_Neyman1.T_se = NaN;    
end

if ~isempty(find(ismember(models,'MonoExpo_Neyman2'),1))
    fitting_results.MonoExpo_Neyman2.A_se = errors2.MonoExpo_Neyman2.A_se;
    fitting_results.MonoExpo_Neyman2.T_se = errors2.MonoExpo_Neyman2.T_se;
else
    fitting_results.MonoExpo_Neyman2.A_se = NaN;
    fitting_results.MonoExpo_Neyman2.T_se = NaN;    
end

if ~isempty(find(ismember(models,'MonoExpo_Pearson'),1))
    fitting_results.MonoExpo_Pearson.A_se = errors2.MonoExpo_Pearson.A_se;
    fitting_results.MonoExpo_Pearson.T_se = errors2.MonoExpo_Pearson.T_se;
else
    fitting_results.MonoExpo_Pearson.A_se = NaN;
    fitting_results.MonoExpo_Pearson.T_se = NaN;    
end

if ~isempty(find(ismember(models,'MonoExpo_fitted_weighting'),1))
    fitting_results.MonoExpo_fitted_weighting.A_se = errors2.MonoExpo_fitted_weighting.A_se;
    fitting_results.MonoExpo_fitted_weighting.T_se = errors2.MonoExpo_fitted_weighting.T_se;
else
    fitting_results.MonoExpo_fitted_weighting.A_se = NaN;
    fitting_results.MonoExpo_fitted_weighting.T_se = NaN;    
end

if ~isempty(find(ismember(models,'DoubleExpo_fitted_weighting'),1))
    fitting_results.DoubleExpo_fitted_weighting.B1_se = errors2.DoubleExpo_fitted_weighting.B1_se;
    fitting_results.DoubleExpo_fitted_weighting.T1_se = errors2.DoubleExpo_fitted_weighting.T1_se;
    fitting_results.DoubleExpo_fitted_weighting.B2_se = errors2.DoubleExpo_fitted_weighting.B2_se;
    fitting_results.DoubleExpo_fitted_weighting.T2_se = errors2.DoubleExpo_fitted_weighting.T2_se;
else
    fitting_results.DoubleExpo_fitted_weighting.B1_se = NaN;
    fitting_results.DoubleExpo_fitted_weighting.T1_se = NaN;
    fitting_results.DoubleExpo_fitted_weighting.B2_se = NaN;
    fitting_results.DoubleExpo_fitted_weighting.T2_se = NaN;    
end


%% clear

clear error
clear errors2
clear binranges;
clear bincounts;
clear data;
clear power;

end

