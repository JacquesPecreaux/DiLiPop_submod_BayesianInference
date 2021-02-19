function [ familyAssignement,probability_long_duration,probability_short_duration,fitting_results ] = ...
    to_assign_MTs_toGivenPopulation_usingProba_mle_doubleExpo__...
    ( dataTracks,fitting_results,name1,name2,familyAssignement,proba_threshold,iEmbryo,frequency  )

% function to assigne the 2 sub-populations when duration distribution
% fitted by DoubleExpo

% to always have T1 < T2
backup_T2 = fitting_results.DoubleExpo.T2;
backup_P2 = fitting_results.DoubleExpo.P2;
if iEmbryo == 0
    backup_R2 = fitting_results.DoubleExpo.R_normalization.double(2,1);
else
    name0 = ['embryo', num2str(iEmbryo)];
    backup_R2 = fitting_results.DoubleExpo.R_normalization.(name0).double(2,1);
end

if fitting_results.DoubleExpo.T1 > fitting_results.DoubleExpo.T2
    fitting_results.DoubleExpo.T2 = fitting_results.DoubleExpo.T1;
    fitting_results.DoubleExpo.T1 = backup_T2;
    fitting_results.DoubleExpo.P2 = fitting_results.DoubleExpo.P1;
    fitting_results.DoubleExpo.P1 = backup_P2;
    if iEmbryo == 0
        fitting_results.DoubleExpo.R_normalization.double(2,1) = fitting_results.DoubleExpo.R_normalization.double(1,1);
        fitting_results.DoubleExpo.R_normalization.double(1,1) = backup_R2;
    else
        fitting_results.DoubleExpo.R_normalization.(name0).double(2,1) = fitting_results.DoubleExpo.R_normalization.(name0).double(1,1);
        fitting_results.DoubleExpo.R_normalization.(name0).double(1,1) = backup_R2;
    end
end


%to assigne each MT to the family displaying very short, short or long residency time
%at the cortex

% mi = duration of track i
% Prob(prop1|mi) = probability to choose pop1 knowing value of mi
% Prob(prop1|mi) = Prob(pop1) * Prob(mi|pop1) = proba of pop1 * proba that mi is in Pop 1

% calculate probability function of each population: labbook #9, page 64
%ex: Prob(mi|pop1) = (P1 * size_tot / R1 ) * exp(-mi/T1)
if iEmbryo == 0
    given_size_population = dataTracks.(name1).(name2).numTracks; %fitting_results.size_population_raw ;
    probability_short_duration = @(x) ( given_size_population* fitting_results.DoubleExpo.P1 ) ...
        / ( fitting_results.DoubleExpo.R_normalization.double(1,1)  ) * exp( - x / fitting_results.DoubleExpo.T1 );
    probability_long_duration = @(x) ( given_size_population * fitting_results.DoubleExpo.P2 ) ...
        / ( fitting_results.DoubleExpo.R_normalization.double(2,1)) * exp( - x / fitting_results.DoubleExpo.T2 );
else
    given_size_population = dataTracks.(name1).(name2).numTracks; %fitting_results.size_population.(name0).raw;
    probability_short_duration = @(x) ( given_size_population* fitting_results.DoubleExpo.P1 ) ...
        / ( fitting_results.DoubleExpo.R_normalization.(name0).double(1,1)  ) * exp( - x / fitting_results.DoubleExpo.T1 );
    probability_long_duration = @(x) ( given_size_population * fitting_results.DoubleExpo.P2 ) ...
        / ( fitting_results.DoubleExpo.R_normalization.(name0).double(2,1)) * exp( - x / fitting_results.DoubleExpo.T2 );
end

for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
    len = dataTracks.(name1).(name2).lengthTracks(yourNumberTrack)/frequency; % in sec.
    H_short = probability_short_duration(len);
    H_long = probability_long_duration(len);
    proba_short_percent(yourNumberTrack) = fitting_results.DoubleExpo.P1*H_short /...
        (fitting_results.DoubleExpo.P2*H_long + fitting_results.DoubleExpo.P1*H_short)*100;
    proba_long_percent(yourNumberTrack) = fitting_results.DoubleExpo.P2*H_long /...
        (fitting_results.DoubleExpo.P2*H_long + fitting_results.DoubleExpo.P1*H_short)*100;
end

familyAssignement.probability.short = proba_short_percent;
familyAssignement.probability.long = proba_long_percent;


% allocate variables to zeros
nb_MTs_short_duration = 0;
nb_MTs_long_duration = 0;
index_MTs_short_duration = [];
index_MTs_long_duration = [];
probability_MTs_long_duration = [];
probability_MTs_short_duration = [];

% assignement of each individual MT
alpha_confidence = proba_threshold; % the risk we offer to have

for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
    if  proba_long_percent(yourNumberTrack) > proba_short_percent(yourNumberTrack) % more probability track belongs to family with long duration
        nb_MTs_long_duration = nb_MTs_long_duration +1;
        index_MTs_long_duration(nb_MTs_long_duration) = yourNumberTrack;
        probability_MTs_long_duration(nb_MTs_long_duration) = proba_long_percent(yourNumberTrack);
    else
        nb_MTs_short_duration = nb_MTs_short_duration +1;
        index_MTs_short_duration(nb_MTs_short_duration) = yourNumberTrack;
        probability_MTs_short_duration(nb_MTs_short_duration) = proba_short_percent(yourNumberTrack);
    end
end


familyAssignement.short_duration.number = nb_MTs_short_duration;
familyAssignement.short_duration.index = index_MTs_short_duration;
familyAssignement.short_duration.probability = probability_MTs_short_duration;
familyAssignement.long_duration.number = nb_MTs_long_duration;
familyAssignement.long_duration.index = index_MTs_long_duration;
familyAssignement.long_duration.probability = probability_MTs_long_duration;


end

