function [ V ] = Loglikelihood2_norm( par, data, model, f, size_population_raw, fixed_short_lifetime,fixed_short_percent)

xdata=data(:,1);

if strcmp(model,'MonoExpo')
    
    R = sum( exp(-xdata.*par(1)) );
    
    mod = f(par, data, R, size_population_raw);
    
elseif strcmp(model,'DoubleExpo')
    
    R1 = sum( exp(-xdata.*par(2)) );
    R2 = sum( exp(-xdata.*par(3)) );
    
    mod = f(par, data, R1, R2, size_population_raw);   
    
elseif strcmp(model,'DoubleExpo_fixedT0') && exist('fixed_short_lifetime','var') == 1
    
    R2_1 = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
    R2_2 = sum( exp(-xdata.*par(2)) );
    
    mod = f(par, data, R2_1, R2_2, size_population_raw,fixed_short_lifetime);
    
elseif strcmp(model,'DoubleExpo_fixedP0') && exist('fixed_short_percent','var') == 1
    
    R2_1 = sum( exp(-xdata.*par(1)) );
    R2_2 = sum( exp(-xdata.*par(2)) );
    
    mod = f(par, data, R2_1, R2_2, size_population_raw,fixed_short_percent);
    
elseif strcmp(model,'DoubleExpo_fixedT0P0') && exist('fixed_short_lifetime','var') == 1 && exist('fixed_short_percent','var') == 1
    
    R2_1_ = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
    R2_2_ = sum( exp(-xdata.*par(1)) );
    
    mod = f(par, data, R2_1_, R2_2_, size_population_raw,fixed_short_lifetime,fixed_short_percent);    
    
elseif strcmp(model,'MonoExpo_stretched')
    
    RS = sum( exp( -(xdata.*par(1)).^(1/par(2)) ) );
    
    mod = f(par, data, RS, size_population_raw);       
    
elseif strcmp(model,'TripleExpo')
    
    R3_1 = sum( exp(-xdata.*par(2)) );
    R3_2 = sum( exp(-xdata.*par(4)) );
    R3_3 = sum( exp(-xdata.*par(5)) );
    
    mod = f(par, data, R3_1, R3_2, R3_3, size_population_raw);  
    
elseif strcmp(model,'TripleExpo_fixedT0') && exist('fixed_short_lifetime','var') == 1
    
    R3_1 = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
    R3_2 = sum( exp(-xdata.*par(3)) );
    R3_3 = sum( exp(-xdata.*par(4)) );
    
    mod = f(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_short_lifetime); 
        
elseif strcmp(model,'TripleExpo_fixedP0') && exist('fixed_short_percent','var') == 1
    
    R3_1 = sum( exp(-xdata.*par(1)) );
    R3_2 = sum( exp(-xdata.*par(3)) );
    R3_3 = sum( exp(-xdata.*par(4)) );
    
    mod = f(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_short_percent); 
    
elseif strcmp(model,'TripleExpo_fixedT0P0') && exist('fixed_short_lifetime','var') == 1 && exist('fixed_short_percent','var') == 1
    
    R3_1 = sum( exp(-xdata.*(1/fixed_short_lifetime)) );
    R3_2 = sum( exp(-xdata.*par(2)) );
    R3_3 = sum( exp(-xdata.*par(3)) );
    
    mod = f(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_short_lifetime,fixed_short_percent); 
    
elseif strcmp(model,'TripleExpo_fixedT1T2') && exist('fixed_lifetime1','var') == 1 && exist('fixed_lifetime2','var') == 1
    
    R3_1 = sum( exp(-xdata.*par(2)) );
    R3_2 = sum( exp(-xdata.*(1/fixed_lifetime1)) );
    R3_3 = sum( exp(-xdata.*(1/fixed_lifetime2)) );
    
    mod = f(par, data, R3_1, R3_2, R3_3, size_population_raw,fixed_lifetime1,fixed_lifetime2); 
 
elseif strcmp(model,'QuadroExpo')
    
    R4_1 = sum( exp(-xdata.*par(2)) );
    R4_2 = sum( exp(-xdata.*par(4)) );
    R4_3 = sum( exp(-xdata.*par(6)) );
    R4_4 = sum( exp(-xdata.*par(7)) );    
    
    mod = f(par, data, R4_1, R4_2, R4_3, R4_4, size_population_raw);    
    
elseif strcmp(model,'Drift_diffusion')
    
    Rd = sum( xdata.^(-3/2) .* exp(-xdata.*par(1)) );
    
    mod = f(par, data, Rd, size_population_raw);
    
end

V = sum (mod - data(:,2).*log(mod));


end

