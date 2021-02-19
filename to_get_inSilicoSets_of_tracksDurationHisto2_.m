function data_simulation = to_get_inSilicoSets_of_tracksDurationHisto2_...
    ( fitting_results,models,blastomere,timing_phase,nbEmbryo_givenCondition,tmin,fixed_short_lifetime )

frequency = 10; %Hz
tmax = 10; % Maximal duration of a track in sec
%tmin = 0.3;

%% to generate duraation distribution

for iEmbryo = 1 : nbEmbryo_givenCondition
    
    %binranges = [tmin+0.05 : 0.1 : tmax-0.05].*frequency;
    binranges = [tmin : 0.1 : tmax].*frequency;
    name_embryo = ['embryo' num2str(iEmbryo)];
    
    if ~isempty(find(ismember(models,'MonoExpo'),1))
        
        % Parameters of simulation
        np = round(fitting_results.(blastomere).(timing_phase).size_population_raw);
        tau = fitting_results.(blastomere).(timing_phase).MonoExpo.T; %tau = 1;
        
        
        data_p=zeros(1,np);
        
        % notice: round to lowest closed integer to mimick tracks durations
        for i=1:np
            tmp = floor( exprnd(tau*frequency) );
            while tmp < tmin*frequency || tmp > tmax*frequency
                tmp = floor( exprnd(tau*frequency) );
            end
            data_p(i) = tmp;
        end
        
        [bincounts] = hist(data_p, binranges);
        %[bincounts] = histc(data_p, binranges);
        clear data_p
        
        
    elseif ~isempty(find(ismember(models,'DoubleExpo'),1))
        
        % parameters of simulation
        np = round(fitting_results.(blastomere).(timing_phase).size_population_raw);
        prop = fitting_results.(blastomere).(timing_phase).DoubleExpo.P1;
        tau1 = fitting_results.(blastomere).(timing_phase).DoubleExpo.T1;
        tau2 = fitting_results.(blastomere).(timing_phase).DoubleExpo.T2;
        
        data_p=zeros(1,np);
        
        for i=1:np
            r=rand(1);
            if (r<prop)
                tmp = floor( exprnd(tau1*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau1*frequency) );
                end
                data_p(i) = tmp;
            else
                tmp = floor( exprnd(tau2*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau2*frequency) );
                end
                data_p(i) = tmp;
            end
        end
        
        [bincounts] = hist(data_p, binranges);
        %[bincounts] = histc(data_p, binranges);
        clear data_p
        
        
    elseif ~isempty(find(ismember(models,'TripleExpo'),1))
        
        % parameters of simulation
        np = round(fitting_results.(blastomere).(timing_phase).size_population_raw);
        prop1 = fitting_results.(blastomere).(timing_phase).TripleExpo.PP1;
        prop2 = fitting_results.(blastomere).(timing_phase).TripleExpo.PP2;
        tau1 = fitting_results.(blastomere).(timing_phase).TripleExpo.TT1;
        tau2 = fitting_results.(blastomere).(timing_phase).TripleExpo.TT2;
        tau3 = fitting_results.(blastomere).(timing_phase).TripleExpo.TT3;
        
        data_p=zeros(1,np);
        
        for i=1:np
            r1=rand(1);
            if (r1<prop1+prop2)
                r2=rand(1)*(prop1+prop2);
                if (r2<prop1)
                    tmp = floor( exprnd(tau1*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau1*frequency) );
                    end
                    data_p(i) = tmp;
                else
                    tmp = floor( exprnd(tau2*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau2*frequency) );
                    end
                    data_p(i) = tmp;
                end
                
            else
                tmp = floor( exprnd(tau3*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau3*frequency) );
                end
                data_p(i) = tmp;
            end
        end
        
        [bincounts] = hist(data_p, binranges);
        %[bincounts] = histc(data_p, binranges);
        clear data_p
        
        
    elseif ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1)) && exist('fixed_short_lifetime','var') == 1
        
        
        % parameters of simulation
        np = round(fitting_results.(blastomere).(timing_phase).size_population_raw);
        prop0__ = fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P0;
        prop1__ = fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.P1;
        %tau0__ = fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T0;
        tau1__ = fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T1;
        tau2__ = fitting_results.(blastomere).(timing_phase).TripleExpo_fixedT0.T2;
        
        data_p=zeros(1,np);
        
        for i=1:np
            r1=rand(1);
            if (r1<prop0__+prop1__)
                r2=rand(1)*(prop0__+prop1__);
                if (r2<prop0__)
                    tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                    end
                    data_p(i) = tmp;
                else
                    tmp = floor( exprnd(tau1__*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau1__*frequency) );
                    end
                    data_p(i) = tmp;
                end
                
            else
                tmp = floor( exprnd(tau2__*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau2__*frequency) );
                end
                data_p(i) = tmp;
            end
        end
        
        [bincounts] = hist(data_p, binranges);
        %[bincounts] = histc(data_p, binranges);
        clear data_p
    
           
    elseif ~isempty(find(ismember(models,'TripleExpo'),1))
        
        % parameters of simulation
        np = round(fitting_results.(blastomere).(timing_phase).size_population_raw);
        prop1 = fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP1;
        prop2 = fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP2;
        prop3 = fitting_results.(blastomere).(timing_phase).QuadroExpo.PPP3;
        tau1 = fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT1;
        tau2 = fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT2;
        tau3 = fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT3;
        tau4 = fitting_results.(blastomere).(timing_phase).QuadroExpo.TTT4;
        
        data_p=zeros(1,np);
        
        for i=1:np
            r1=rand(1);
            if (r1<prop1+prop2+prop3)
                r2=rand(1)*(prop1+prop2+prop3);
                if (r2<prop1+prop2)
                    r3=rand(1)*(prop1+prop2);
                    if r3 < prop1
                        tmp = floor( exprnd(tau1*frequency) );
                        while tmp < tmin*frequency || tmp > tmax*frequency
                            tmp = floor( exprnd(tau1*frequency) );
                        end
                        data_p(i) = tmp;
                    else
                        tmp = floor( exprnd(tau2*frequency) );
                        while tmp < tmin*frequency || tmp > tmax*frequency
                            tmp = floor( exprnd(tau2*frequency) );
                        end
                        data_p(i) = tmp;
                    end
                else
                    tmp = floor( exprnd(tau3*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau3*frequency) );
                    end
                    data_p(i) = tmp;
                end
            else
                tmp = floor( exprnd(tau4*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau3*frequency) );
                end
                data_p(i) = tmp;
            end
        end
        
        [bincounts] = hist(data_p, binranges);
        %[bincounts] = histc(data_p, binranges);
        clear data_p
        
    end
    
    
    %----------------------
    nb_zero = 0;
    index_zero = [];
    for i = 1 : length(bincounts)
        if bincounts(i) ==0
            nb_zero = nb_zero+1;
            index_zero = [[index_zero] i ];
        end
        if nb_zero == 10
            max_index = i;
            break
        end
    end
    
    if nb_zero < 10
        max_index = max(index_zero);
    end
    
    max_max_index = length(bincounts);
    index_to_remove2 = [(max_index):(max_max_index)];
    binranges(index_to_remove2) = [];
    bincounts(index_to_remove2) = [];
    
    nb_one = 0;
    index_one = [];
    for i = 1 : length(bincounts)
        if bincounts(i) ==1
            nb_one = nb_one+1;
            index_one = [[index_one] i ];
        end
        if nb_one == 3
            max_index_one = i;
            break
        end
    end
    
    if nb_one < 3
        max_index_one = max(index_one);
    end
    
    max_max_index_one = length(bincounts);
    
    bincounts_above_max_one = 0;
    for index = max_index_one+1 : max_max_index_one
        bincounts_above_max_one = bincounts_above_max_one + bincounts(index);
    end
    
    if max_max_index_one >= (max_index_one + 1)
        bincounts_above_max_one = bincounts_above_max_one / (max_max_index_one - max_index_one - 1);
        
        if ~isnan(bincounts_above_max_one)
            index_to_remove3 = [(max_index_one):(max_max_index_one)];
            binranges(index_to_remove3) = [];
            bincounts(index_to_remove3) = [];
            binranges(max_index_one) = max_index_one + ( max_max_index_one - max_index_one)/2;
            bincounts(max_index_one) = bincounts_above_max_one;
        end
    end
    
    if isinf(bincounts(:))
        disp('problem')
    end
    index_to_remove_inf = find(isinf(bincounts));
    binranges(index_to_remove_inf) = [];
    bincounts(index_to_remove_inf) = [];
    
    %--------------------
    data = [binranges(:)./frequency bincounts(:)];
    data_simulation.(name_embryo).data = data;
    clear data
    
end


end





