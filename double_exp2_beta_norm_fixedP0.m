function y = double_exp2_beta_norm_fixedP0(par, data, R1, R2,size_population_raw,fixed_short_percent)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( fixed_short_percent/R1*exp(-xdata.*par(1))  + (1 -fixed_short_percent)/R2*exp(-xdata.*par(2)) ) ;

end