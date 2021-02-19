function y = simple_exp2_stretched_beta_norm(par,data,RS,size_population_raw)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw / RS *exp( -(xdata.*par(1)).^(1/par(2)) );

end