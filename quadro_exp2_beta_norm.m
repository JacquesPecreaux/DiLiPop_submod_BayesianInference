function y = quadro_exp2_beta_norm(par, data, R1, R2, R3, R4, size_population_raw)

% Function parametrized with exp(-t*beta) and normalized by sum(
% exp(-t*beta) )

xdata=data(:,1);

y = size_population_raw * ( par(1)/R1*exp(-xdata.*par(2))  + par(3)/R2*exp(-xdata.*par(4)) + ...
    par(5)/R3*exp(-xdata.*par(6))  + ( 1 - par(1) - par(3) - par(5) )/R4*exp(-xdata.*par(7)) ) ;

end