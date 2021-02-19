function [ V_sum ] = Loglikelihood2_total( par, input, model, f, nbEmbryo_givenCondition, size_population, fixed_short_lifetime, fixed_short_percent )

% to calculate the sum of the likelihood of the different embryos that will
% be necessary to find the fitting parameters that are identical for all
% the embryos
    
V_sum = 0;

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    
    name_embryo = ['embryo' num2str(iEmbryo)];
    
    data =  input.(name_embryo).data;
    size_population_raw = size_population.(name_embryo).raw;
    
    if strcmp(model,'TripleExpo_fixedT0P0') && exist('fixed_short_lifetime','var') == 1 &&  exist('fixed_short_percent','var') == 1
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, fixed_short_lifetime, fixed_short_percent );
    elseif strcmp(model,'DoubleExpo_fixedT0P0') && exist('fixed_short_lifetime','var') == 1 &&  exist('fixed_short_percent','var') == 1
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, fixed_short_lifetime, fixed_short_percent );
    elseif strcmp(model,'TripleExpo_fixedT0') && exist('fixed_short_lifetime','var') == 1 
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, fixed_short_lifetime );
    elseif strcmp(model,'DoubleExpo_fixedT0') && exist('fixed_short_lifetime','var') == 1 
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, fixed_short_lifetime );    
    elseif strcmp(model,'DoubleExpo_fixedP0') && exist('fixed_short_percent','var') == 1 
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, [], fixed_short_percent );  
    elseif strcmp(model,'TripleExpo_fixedP0') && exist('fixed_short_percent','var') == 1 
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, [], fixed_short_percent );         
    else
        [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw );    
    end
    
    V_sum = V_sum + V ;

end


end

