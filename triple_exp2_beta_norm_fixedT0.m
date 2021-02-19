function y = triple_exp2_beta_norm_fixedT0(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_short_lifetime)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( par(1)/R3_1*exp(-xdata./fixed_short_lifetime)  + par(2)/R3_2*exp(-xdata.*par(3)) + ...
    ( 1 - par(1) - par(2) ) /R3_3*exp(-xdata.*par(4)) );

end