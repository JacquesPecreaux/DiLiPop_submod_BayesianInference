function y = double_exp2_beta_norm_fixedT0(par, data, R1, R2,size_population_raw,fixed_short_lifetime)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( par(1)/R1*exp(-xdata./fixed_short_lifetime)  + (1 -par(1))/R2*exp(-xdata.*par(2)) ) ;

end