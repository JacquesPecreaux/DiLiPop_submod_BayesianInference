function y = double_exp2_beta_norm(par, data, R1, R2,size_population_raw)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( par(1)/R1*exp(-xdata.*par(2))  + (1 -par(1))/R2*exp(-xdata.*par(3)) ) ;

end