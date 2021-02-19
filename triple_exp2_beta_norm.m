function y = triple_exp2_beta_norm(par, data, R3_1, R3_2, R3_3, size_population_raw)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( par(1)/R3_1*exp(-xdata.*par(2))  + par(3)/R3_2*exp(-xdata.*par(4)) + ( 1 - par(1) - par(3) ) /R3_3*exp(-xdata.*par(5)) );

end