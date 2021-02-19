function y = simple_exp2_beta_norm(par,data,R,size_population_raw)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw/R*exp(-xdata.*par(1)) ;

end