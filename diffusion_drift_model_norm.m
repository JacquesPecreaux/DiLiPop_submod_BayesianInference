function y = diffusion_drift_model_norm(par,data,Rd,size_population_raw)

% Function parametrized with t^(-3/2) * exp(-t*beta) and normalized by sum(
% t^(-3/2) * exp(-t*beta) )
% see paper from Bicout 1997

xdata=data(:,1);

y = size_population_raw/Rd* xdata.^(-3/2) .* exp(-xdata.*par(1)) ;

end