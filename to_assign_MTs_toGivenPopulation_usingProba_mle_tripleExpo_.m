function [ familyAssignement,probability_long_duration,probability_short_duration,probability_veryshort_duration,fitting_results ] = ...
    to_assign_MTs_toGivenPopulation_usingProba_mle_tripleExpo_...
    ( dataTracks,fitting_results,name1,name2,familyAssignement,proba_threshold,iEmbryo,frequency )

% function to assigne the 3 sub-populations when duration distribution
% fitted by TripleExpo

% to always have T1 < T3
backup_TT3 = fitting_results.TripleExpo.TT3;
backup_PP3 = fitting_results.TripleExpo.PP3;
if iEmbryo == 0
    backup_R3 = fitting_results.TripleExpo.R_normalization.triple(3,1);
else
    name0 = ['embryo', num2str(iEmbryo)];
    backup_R3 = fitting_results.TripleExpo.R_normalization.(name0).triple(3,1);
end

if fitting_results.TripleExpo.TT1 > fitting_results.TripleExpo.TT3
    fitting_results.TripleExpo.TT3 = fitting_results.TripleExpo.TT1;
    fitting_results.TripleExpo.TT1 = backup_TT3;
    fitting_results.TripleExpo.PP3 = fitting_results.TripleExpo.PP1;
    fitting_results.TripleExpo.PP1 = backup_PP3;
    if iEmbryo == 0
        fitting_results.TripleExpo.R_normalization.triple(3,1) = fitting_results.TripleExpo.R_normalization.triple(1,1);
        fitting_results.TripleExpo.R_normalization.triple(1,1) = backup_R3;
    else
        fitting_results.TripleExpo.R_normalization.(name0).triple(3,1) = fitting_results.TripleExpo.R_normalization.(name0).triple(1,1);
        fitting_results.TripleExpo.R_normalization.(name0).triple(1,1) = backup_R3;
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
    given_size_population = fitting_results.size_population_raw ;
    probability_long_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP3 ) ...
    / ( fitting_results.TripleExpo.R_normalization.triple(3,1)  ) ...
    * exp( - x / fitting_results.TripleExpo.TT3 );
probability_short_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP2 ) ...
    / ( fitting_results.TripleExpo.R_normalization.triple(2,1)  ) ...
    * exp( - x / fitting_results.TripleExpo.TT2 );
probability_veryshort_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP1 ) ...
    / ( fitting_results.TripleExpo.R_normalization.triple(1,1)) ...
    * exp( - x / fitting_results.TripleExpo.TT1 );
else
    given_size_population = fitting_results.size_population.(name0).raw;
    probability_long_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP3 ) ...
    / ( fitting_results.TripleExpo.R_normalization.(name0).triple(3,1)  ) ...
    * exp( - x / fitting_results.TripleExpo.TT3 );
probability_short_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP2 ) ...
    / ( fitting_results.TripleExpo.R_normalization.(name0).triple(2,1)  ) ...
    * exp( - x / fitting_results.TripleExpo.TT2 );
probability_veryshort_duration = @(x) ( given_size_population * fitting_results.TripleExpo.PP1 ) ...
    / ( fitting_results.TripleExpo.R_normalization.(name0).triple(1,1)) ...
    * exp( - x / fitting_results.TripleExpo.TT1 );
end


% Prob(P1|mi) = P1 * Prob(mi|pop1) / ( P1 * Prob(mi|pop1) + P2 * Prob(mi|pop2) + P3 * Prob(mi|pop2) )
for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
    len = dataTracks.(name1).(name2).lengthTracks(yourNumberTrack)/frequency; % in sec.
    H_veryshort = probability_veryshort_duration(len);
    H_short = probability_short_duration(len);
    H_long = probability_long_duration(len);
    proba_veryshort_percent(yourNumberTrack) = fitting_results.TripleExpo.PP1*H_veryshort /...
        ( fitting_results.TripleExpo.PP3*H_long + fitting_results.TripleExpo.PP2*H_short + ...
        fitting_results.TripleExpo.PP1*H_veryshort )*100;
    proba_short_percent(yourNumberTrack) = fitting_results.TripleExpo.PP2*H_short /...
        ( fitting_results.TripleExpo.PP3*H_long + fitting_results.TripleExpo.PP2*H_short + ...
        fitting_results.TripleExpo.PP1*H_veryshort )*100;
    proba_long_percent(yourNumberTrack) = fitting_results.TripleExpo.PP3*H_long /...
        ( fitting_results.TripleExpo.PP3*H_long + fitting_results.TripleExpo.PP2*H_short + ...
        fitting_results.TripleExpo.PP1*H_veryshort )*100;
end

familyAssignement.probability.veryshort = proba_veryshort_percent;
familyAssignement.probability.short = proba_short_percent;
familyAssignement.probability.long = proba_long_percent;


% allocate variables to zeros

nb_MTs_veryshort_duration = 0;
nb_MTs_short_duration = 0;
nb_MTs_long_duration = 0;
nb_MTs_not_classified = 0;
index_MTs_veryshort_duration = [];
index_MTs_short_duration = [];
index_MTs_long_duration = [];
index_MTs_not_classified = [];
probability_MTs_long_duration = [];
probability_MTs_short_duration = [];
probability_MTs_veryshort_duration = [];
probability_MTs_not_classified_1 = [];
probability_MTs_not_classified_2 = [];

% assignement of each individual MT

alpha_confidence = proba_threshold; % the risk we offer to have

for yourNumberTrack = 1 : dataTracks.(name1).(name2).numTracks
    if  proba_long_percent(yourNumberTrack) > alpha_confidence * 100 % more probability track belongs to family with long duration
        nb_MTs_long_duration = nb_MTs_long_duration +1;
        index_MTs_long_duration(nb_MTs_long_duration) = yourNumberTrack;
        probability_MTs_long_duration(nb_MTs_long_duration) = proba_long_percent(yourNumberTrack);
    elseif proba_short_percent(yourNumberTrack) > alpha_confidence * 100 % more probability track belongs to family with short duration
        nb_MTs_short_duration = nb_MTs_short_duration +1;
        index_MTs_short_duration(nb_MTs_short_duration) = yourNumberTrack;
        probability_MTs_short_duration(nb_MTs_short_duration) = proba_short_percent(yourNumberTrack);
    elseif proba_veryshort_percent(yourNumberTrack) > alpha_confidence * 100 % more probability track belongs to family with short duration
        nb_MTs_veryshort_duration = nb_MTs_veryshort_duration +1;
        index_MTs_veryshort_duration(nb_MTs_veryshort_duration) = yourNumberTrack;
        probability_MTs_veryshort_duration(nb_MTs_veryshort_duration) = proba_veryshort_percent(yourNumberTrack);
    else
        nb_MTs_not_classified = nb_MTs_not_classified +1;
        index_MTs_not_classified(nb_MTs_not_classified) = yourNumberTrack;
        probability_MTs_not_classified_0(nb_MTs_not_classified) = proba_veryshort_percent(yourNumberTrack);
        probability_MTs_not_classified_1(nb_MTs_not_classified) = proba_short_percent(yourNumberTrack);
        probability_MTs_not_classified_2(nb_MTs_not_classified) = proba_long_percent(yourNumberTrack);
    end
end


familyAssignement.veryshort_duration.number = nb_MTs_veryshort_duration;
familyAssignement.veryshort_duration.index = index_MTs_veryshort_duration;
familyAssignement.veryshort_duration.probability = probability_MTs_veryshort_duration;
familyAssignement.short_duration.number = nb_MTs_short_duration;
familyAssignement.short_duration.index = index_MTs_short_duration;
familyAssignement.short_duration.probability = probability_MTs_short_duration;
familyAssignement.long_duration.number = nb_MTs_long_duration;
familyAssignement.long_duration.index = index_MTs_long_duration;
familyAssignement.long_duration.probability = probability_MTs_long_duration;
familyAssignement.not_classified.number = nb_MTs_not_classified;
familyAssignement.not_classified.index = index_MTs_not_classified;
if nb_MTs_not_classified > 0
    familyAssignement.not_classified.probability_veryshort = probability_MTs_not_classified_0;
    familyAssignement.not_classified.probability_short = probability_MTs_not_classified_1;
    familyAssignement.not_classified.probability_long = probability_MTs_not_classified_2;
end


end

