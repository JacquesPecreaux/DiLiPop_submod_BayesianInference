function data_simulation = to_get_inSilicoSets_of_tracksDurationHisto...
    ( fitting_results,models,nbEmbryo_givenCondition,tmin,fixed_short_lifetime,fixed_short_percent)

global param
global general_param

if isempty(param)
    frequency = 10;
else
    frequency = param.sp6; %Hz
end
tmax = 10; % in sec
if tmin >= 1
    tmin = general_param.cortex_analysis.minLength / param.sp6;
end

tau1 = []; % double expo
tau2 = []; % double expo
tau = []; % mono expo
tau1_ = []; % triple expo
tau2_ = []; % triple expo
tau3_ = []; % triple expo
prop = []; % double expo
prop1_= []; % triple expo
prop2_ = []; % triple expo
tau1__ = []; % triple expo with fixed T0
tau2__ = []; % triple expo with fixed T0
prop0__ = []; % triple expo with fixed T0
prop1__ = []; % triple expo with fixed T0
tau1___ = []; % triple expo with fixed T0P0
tau2___ = []; % triple expo with fixed T0P0
prop1___ = []; % triple expo with fixed T0P0
tau2_2 = []; % double expo with fixed T0
prop_ = []; % double expo with fixed T0
tau2_2_ = []; % double expo with fixed T0P0
% quadro
tau4_1 = [];
tau4_2 = [];
tau4_3 = [];
tau4_4 = [];
prop4_1 = [];
prop4_2 = [];
prop4_3 = [];

if ~exist('fixed_short_lifetime','var')
    fixed_short_lifetime = [];
end
if ~exist('fixed_short_percent','var')
    fixed_short_percent = [];
end

%%  to generate duraation distribution

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    name_embryo = ['embryo' num2str(iEmbryo)];
    
   % Parameters of simulation
    np = round(fitting_results.size_population.(name_embryo).raw);
    
    if strcmp(models,'MonoExpo')
        tau = fitting_results.MonoExpo.T;
    elseif strcmp(models,'DoubleExpo')
        prop = fitting_results.DoubleExpo.P1;
        tau1 = fitting_results.DoubleExpo.T1;
        tau2 = fitting_results.DoubleExpo.T2;
    elseif strcmp(models,'DoubleExpo_fixedT0')
        prop_ = fitting_results.DoubleExpo_fixedT0.P1;
        tau2_2 = fitting_results.DoubleExpo_fixedT0.T2;       
    elseif strcmp(models,'DoubleExpo_fixedT0P0')
        tau2_2_ = fitting_results.DoubleExpo_fixedT0P0.T2;         
    elseif strcmp(models,'TripleExpo')
        prop1_ = fitting_results.TripleExpo.PP1;
        prop2_ = fitting_results.TripleExpo.PP2;
        tau1_ = fitting_results.TripleExpo.TT1;
        tau2_ = fitting_results.TripleExpo.TT2;
        tau3_ = fitting_results.TripleExpo.TT3;
    elseif strcmp(models,'TripleExpo_fixedT0')
        prop0__ = fitting_results.TripleExpo_fixedT0.P0;
        prop1__ = fitting_results.TripleExpo_fixedT0.P1;
        tau1__ = fitting_results.TripleExpo_fixedT0.T1;
        tau2__ = fitting_results.TripleExpo_fixedT0.T2;
    elseif strcmp(models,'TripleExpo_fixedT0P0')
        prop1___ = fitting_results.TripleExpo_fixedT0P0.P1;
        tau1___ = fitting_results.TripleExpo_fixedT0P0.T1;
        tau2___ = fitting_results.TripleExpo_fixedT0P0.T2;
        
    elseif strcmp(models,'QuadroExpo')        % quadro
        tau4_1 = fitting_results.QuadroExpo.TTT1;
        tau4_2 = fitting_results.QuadroExpo.TTT2;
        tau4_3 = fitting_results.QuadroExpo.TTT3;
        tau4_4 = fitting_results.QuadroExpo.TTT4;
        prop4_1 = fitting_results.QuadroExpo.PPP1;
        prop4_2 = fitting_results.QuadroExpo.PPP2;
        prop4_3 = fitting_results.QuadroExpo.PPP3;
    end
    
    [ data,bincounts,binranges ] = to_generate_inSilico_durationSet...
        ( models,np,frequency,tmin,tmax,tau,prop,tau1,tau2,prop_,tau2_2,tau2_2_,prop1_,prop2_,tau1_,tau2_,tau3_,prop0__,prop1__,fixed_short_lifetime,...
        tau1__,tau2__,fixed_short_percent,prop1___,tau1___,tau2___,prop4_1,prop4_2,prop4_3,tau4_1,tau4_2,tau4_3,tau4_4);
    
    data_simulation.(name_embryo).data = double(data);
    clear data
      
end


end

