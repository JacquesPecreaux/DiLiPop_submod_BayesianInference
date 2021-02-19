    function [ data,bincounts,binranges,data_p ] = to_generate_inSilico_durationSet...
        ( model_init,np,frequency,tmin,tmax,tau,prop,tau1,tau2,prop_,tau2_2,tau2_2_,prop1_,prop2_,tau1_,tau2_,tau3_,prop0__,prop1__,fixed_short_lifetime,...
        tau1__,tau2__,fixed_short_percent,prop1___,tau1___,tau2___,prop4_1,prop4_2,prop4_3,tau4_1,tau4_2,tau4_3,tau4_4,...
        tauS,powerS,use_threshold_count_fit,use_smoothing_fit,binCount_threshold,binCount_smoothW_size,extension)

    % script to create tracks duration distribution according to the model
    % chosen
    % output data_p argument in frame number, while binranges in seconds

    % create the tracks distribution
    % exprnd(mu) generates random numbers from the exponential distribution with mean parameter mu
    % exp(-x/mu)

    global pathMainDirectory

    if nargin < 36
        use_threshold_count_fit = 0;
    end
    if nargin < 37
        use_smoothing_fit = 0;
    end

    binranges = [tmin : 0.1 : tmax].*frequency;

    %% monoExpo

    if ~isempty(find(ismember(model_init,'MonoExpo'),1))

        data_p=zeros(1,np);

        for i=1:np
            tmp = floor( exprnd(tau*frequency) );
            while tmp < tmin*frequency || tmp > tmax*frequency
                tmp = floor( exprnd(tau*frequency) );
            end
            data_p(i) = tmp;
        end

        [bincounts] = hist(data_p, binranges);
    %     clear data_p

    
    %% mono expo stretched
    
    elseif ~isempty(find(ismember(model_init,'MonoExpo_stretched'),1))
                
        binranges = [tmin-0.1 : 0.1 : tmax].*frequency; % minus 0.1 to account for first bin being oversampled
            
        napprox = length(binranges);
        tauS_ = tauS*frequency; % in frame number       
        U = rand(np,1);
        tvec = linspace((tmin-0.1)*frequency,tmax*frequency, napprox);
        f = exp(-(tvec/tauS_).^(1/powerS)) ./ sum(exp(-(tvec/tauS_).^(1/powerS)));
        F = cumsum(f);
        if (exp(-(tmax/tauS_).^(1/powerS))>0.001)
            disp('Increase tmax, F did not reach 1');
        end
        
        [value, indx] = min( abs(F - U),  [], 2 );
        data_p = tvec(indx);
        [bincounts] = hist(data_p, binranges);
        clear value indx F U f tvec
        bincounts = bincounts(2:end); % use only data from tmin
        binranges = binranges(2:end);

    
    %% Dble Expo

    elseif ~isempty(find(ismember(model_init,'DoubleExpo'),1))

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
    %     clear data_p

    %% Dble Expo with fixed T0

    elseif ~isempty(find(ismember(model_init,'DoubleExpo_fixedT0'),1))

        data_p=zeros(1,np);

        for i=1:np
            r=rand(1);
            if (r<prop_)
                tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                end
                data_p(i) = tmp;
            else
                tmp = floor( exprnd(tau2_2*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau2_2*frequency) );
                end
                data_p(i) = tmp;
            end
        end

        [bincounts] = hist(data_p, binranges);

    %% Dble Expo with fixed T0P0

    elseif ~isempty(find(ismember(model_init,'DoubleExpo_fixedT0P0'),1))

        data_p=zeros(1,np);

        for i=1:np
            r=rand(1);
            if (r<fixed_short_percent)
                tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                end
                data_p(i) = tmp;
            else
                tmp = floor( exprnd(tau2_2_*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau2_2_*frequency) );
                end
                data_p(i) = tmp;
            end
        end

        [bincounts] = hist(data_p, binranges);

    %% triple expo

    elseif ~isempty(find(ismember(model_init,'TripleExpo'),1))

        data_p=zeros(1,np);

        for i=1:np
            r1=rand(1);
            if (r1<prop1_+prop2_)
                r2=rand(1)*(prop1_+prop2_);
                if (r2<prop1_)
                    tmp = floor( exprnd(tau1_*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau1_*frequency) );
                    end
                    data_p(i) = tmp;
                else
                    tmp = floor( exprnd(tau2_*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau2_*frequency) );
                    end
                    data_p(i) = tmp;
                end

            else
                tmp = floor( exprnd(tau3_*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau3_*frequency) );
                end
                data_p(i) = tmp;
            end
        end

        [bincounts] = hist(data_p, binranges);
    %     clear data_p


    %% triple expo with fixed T0

    elseif ~isempty(find(ismember(model_init,'TripleExpo_fixedT0'),1))

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
    %     clear data_p


    %% triple expo with fixed T0

    elseif ~isempty(find(ismember(model_init,'TripleExpo_fixedT0P0'),1))

        data_p=zeros(1,np);

        for i=1:np
            r1=rand(1);
            if (r1<fixed_short_percent+prop1___)
                r2=rand(1)*(fixed_short_percent+prop1___);
                if (r2<fixed_short_percent)
                    tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(fixed_short_lifetime*frequency) );
                    end
                    data_p(i) = tmp;
                else
                    tmp = floor( exprnd(tau1___*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau1___*frequency) );
                    end
                    data_p(i) = tmp;
                end

            else
                tmp = floor( exprnd(tau2___*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau2___*frequency) );
                end
                data_p(i) = tmp;
            end
        end

        [bincounts] = hist(data_p, binranges);
    %     clear data_p


    %% triple expo

    elseif ~isempty(find(ismember(model_init,'QuadroExpo'),1))

        data_p=zeros(1,np);

        for i=1:np
            r1=rand(1);
            if (r1<prop4_1+prop4_2+prop4_3)
                r2=rand(1)*(prop4_1+prop4_2+prop4_3);
                if (r2<prop4_1+prop4_2)
                    r3 = rand(1)*(prop4_1+prop4_2);
                    if (r3<prop4_1)
                        tmp = floor( exprnd(tau4_1*frequency) );
                        while tmp < tmin*frequency || tmp > tmax*frequency
                            tmp = floor( exprnd(tau4_1*frequency) );
                        end
                        data_p(i) = tmp;
                    else
                        tmp = floor( exprnd(tau4_2*frequency) );
                        while tmp < tmin*frequency || tmp > tmax*frequency
                            tmp = floor( exprnd(tau4_2*frequency) );
                        end
                        data_p(i) = tmp;
                    end                    
                else
                    tmp = floor( exprnd(tau4_3*frequency) );
                    while tmp < tmin*frequency || tmp > tmax*frequency
                        tmp = floor( exprnd(tau4_3*frequency) );
                    end
                    data_p(i) = tmp;
                end

            else
                tmp = floor( exprnd(tau4_4*frequency) );
                while tmp < tmin*frequency || tmp > tmax*frequency
                    tmp = floor( exprnd(tau4_4*frequency) );
                end
                data_p(i) = tmp;
            end
        end

        [bincounts] = hist(data_p, binranges);
    %     clear data_p

    end

    %% "cleaning"/preparation of the distribution to be fitted later

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

    if use_threshold_count_fit == 1
        nb_belowTreshold = 0;
        index_belowThreshold = [];
        for i = 1 : length(bincounts)
            if bincounts(i) <= binCount_threshold
                nb_belowTreshold = nb_belowTreshold+1;
                index_belowThreshold = [[index_belowThreshold] i ];
            end
        end
        if nb_belowTreshold > 0
            binranges(index_belowThreshold) = [];
            bincounts(index_belowThreshold) = [];
        end
    end

    % to perform smoothing of the distribution
    if use_smoothing_fit == 1
        bincounts_smoothed = smooth(bincounts,binCount_smoothW_size);
        figure,
        plot(binranges,bincounts,'xk');
        set(gca,'yscale','log');
        hold all
        plot(binranges,bincounts_smoothed,'or');
        xlabel ('Duration (image nb.)');
        ylabel ('Count (a.u.)');
        figureName = strcat('SmoothedDistribution-', extension, '.fig');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        figureName = strcat('SmoothedDistribution-', extension, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(gcf)
        clear bincounts
        bincounts(1,:) = bincounts_smoothed;
        clear bincounts_smoothed
    end

    % check for any inf in data
    [ ri] = find(~isfinite(bincounts));
    bincounts(ri)= [];
    binranges(ri)= [];

    if isinf(bincounts(:))
        disp('problem')
    end
    index_to_remove_inf = find(isinf(bincounts));
    binranges(index_to_remove_inf) = [];
    bincounts(index_to_remove_inf) = [];

    %% final results

    data = [binranges(:)./frequency bincounts(:)];
    binranges = binranges ./ frequency;


    end

