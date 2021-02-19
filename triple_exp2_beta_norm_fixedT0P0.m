function y = triple_exp2_beta_norm_fixedT0P0(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_short_lifetime,fixed_short_percent)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( fixed_short_percent/R3_1*exp(-xdata./fixed_short_lifetime)  + par(1)/R3_2*exp(-xdata.*par(2)) + ...
    ( 1 - par(1) - fixed_short_percent ) /R3_3*exp(-xdata.*par(3)) );

end