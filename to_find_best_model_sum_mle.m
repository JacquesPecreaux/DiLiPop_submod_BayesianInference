function [ model_choice,best_model ] = to_find_best_model_sum_mle( input,fitting_results,nbEmbryo_givenCondition,models )

% the purpose is to find the best model that will be unique/same for all
% the embryos

% needs to discuss with Yann if the way to calculate the global critera is
% the correct one (the one I think is not correct has been commented)

% what about the cste?? does it affect the choice of the best model? is sit
% correctly calculated?


%% calculate AIC, AICc and BIC of each model

AIC = [];
AICc = [];
BIC = [];

%-----------------
% model = mono-expo

if ~isempty(find(ismember(models,'MonoExpo'),1))
    
    binranges_length = [];
    AIC_mono = 0;
    BIC_mono = 0;
    AICc_mono = 0;
    sum_MLE = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            simple_ = fitting_results.MonoExpo.(name_embryo).simple;
            binranges_length = [ [binranges_length ] length(binranges) ];
            MLE = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE = MLE + simple_(i)-data(i,2)*log(simple_(i)) + cste;
                    %MLE = MLE + simple_(i)-data(i,2)*log(simple_(i));
                    MLEv(i) = simple_(i)-data(i,2)*log(simple_(i));
                else
                    MLE = MLE + simple_(i) ;
                    MLEv(i) = simple_(i) ;
                end
            end
            
            sum_MLE = sum_MLE + MLE;
            model_choice.MonoExpo.(name_embryo).logMLv = MLEv;
            model_choice.MonoExpo.(name_embryo).logML = MLE;
            model_choice.MonoExpo.(name_embryo).AIC = 2*MLE +2*fitting_results.MonoExpo.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.MonoExpo.(name_embryo).AICc = 2*MLE +2*fitting_results.MonoExpo.nparam + 2*fitting_results.MonoExpo.nparam*(fitting_results.MonoExpo.nparam+1)/...
                (length(binranges) - fitting_results.MonoExpo.nparam -1); %wagenmakers & farrell 2004
            model_choice.MonoExpo.(name_embryo).BIC = 2*MLE +fitting_results.MonoExpo.nparam * log(length(binranges) );
            AIC_mono = AIC_mono + model_choice.MonoExpo.(name_embryo).AIC;
            BIC_mono = BIC_mono + model_choice.MonoExpo.(name_embryo).BIC;
            AICc_mono = AICc_mono + model_choice.MonoExpo.(name_embryo).AICc;
            
            clear simple_
            clear data
            clear binranges
            clear cste
            clear Se
            clear MLEv
        end
        
    end
    
    % commented lines below since I think this is not correct
    %     model_choice.MonoExpo.AIC = 2*sum_MLE +2*fitting_results.MonoExpo.nparam; % AIC = -2 log(MLE) + 2 * nb_param
    %     model_choice.MonoExpo.AICc = 2*sum_MLE +2*fitting_results.MonoExpo.nparam + 2*fitting_results.MonoExpo.nparam*(fitting_results.MonoExpo.nparam+1)/...
    %        (mean(binranges_length) - fitting_results.MonoExpo.nparam -1); %wagenmakers & farrell 2004
    %     model_choice.MonoExpo.BIC = 2*sum_MLE +fitting_results.MonoExpo.nparam * log(mean(binranges_length) );
    model_choice.MonoExpo.logML = sum_MLE;
    model_choice.MonoExpo.AIC = AIC_mono;
    model_choice.MonoExpo.AICc = AICc_mono;
    model_choice.MonoExpo.BIC = BIC_mono;
    AIC = [ [AIC] model_choice.MonoExpo.AIC];
    AICc = [ [AICc] model_choice.MonoExpo.AICc];
    BIC = [ [BIC] model_choice.MonoExpo.BIC];
    
    clear MLE
    clear binranges_length
    clear AIC_mono AICc_mono BIC_mono
    
else
    AIC = [ [AIC] NaN];
    AICc = [ [AICc] NaN];
    BIC = [ [BIC]  NaN];
end


%-----------------
% model = mono-expo stretched

if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    
    binranges_length = [];
    AIC_monoS = 0;
    BIC_monoS = 0;
    AICc_monoS = 0;
    sum_MLES =0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            simple_stretched_ = fitting_results.MonoExpo_stretched.(name_embryo).simpleS;
            MLES = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLES = MLES + simple_stretched_(i)-data(i,2)*log(simple_stretched_(i)) + cste;
                    %MLES = MLES + simple_stretched_(i)-data(i,2)*log(simple_stretched_(i));
                    MLESv(i) = simple_stretched_(i)-data(i,2)*log(simple_stretched_(i));
                else
                    MLES = MLES + simple_stretched_(i) ;
                    MLESv(i) = simple_stretched_(i) ;
                end
            end
            
            sum_MLES = sum_MLES + MLES;
            model_choice.MonoExpo_stretched.(name_embryo).logMLv = MLESv;
            model_choice.MonoExpo_stretched.(name_embryo).logML = MLES;
            model_choice.MonoExpo_stretched.(name_embryo).AIC = 2*MLES +2*fitting_results.MonoExpo_stretched.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.MonoExpo_stretched.(name_embryo).AICc = 2*MLES +2*fitting_results.MonoExpo_stretched.nparam + ...
                2*fitting_results.MonoExpo_stretched.nparam*(fitting_results.MonoExpo_stretched.nparam+1)/(length(binranges) - fitting_results.MonoExpo_stretched.nparam -1); %wagenmakers & farrell 2004
            model_choice.MonoExpo_stretched.(name_embryo).BIC = 2*MLES +fitting_results.MonoExpo_stretched.nparam * log(length(binranges) );
            AIC_monoS = AIC_monoS + model_choice.MonoExpo_stretched.(name_embryo).AIC;
            BIC_monoS = BIC_monoS + model_choice.MonoExpo_stretched.(name_embryo).BIC;
            AICc_monoS = AICc_monoS + model_choice.MonoExpo_stretched.(name_embryo).AICc;
            
            clear simple_stretched_
            clear data
            clear binranges
            clear cste
            clear Se
            clear MLESv
        end
        
    end
    
    % commented lines below since I think this is not correct
    %     model_choice.MonoExpo_stretched.AIC = 2*sum_MLES +2*fitting_results.MonoExpo_stretched.nparam;
    %     model_choice.MonoExpo_stretched.AICc = 2*sum_MLES +2*fitting_results.MonoExpo_stretched.nparam + ...
    %         2*fitting_results.MonoExpo_stretched.nparam*(fitting_results.MonoExpo_stretched.nparam+1)/...
    %         (mean(binranges_length)  - fitting_results.MonoExpo_stretched.nparam -1); %wagenmakers & farrell 2004
    %     model_choice.MonoExpo_stretched.BIC = 2*sum_MLES +fitting_results.MonoExpo_stretched.nparam * log(mean(binranges_length) );
    model_choice.MonoExpo_stretched.logML = sum_MLES;
    model_choice.MonoExpo_stretched.AIC = AIC_monoS;
    model_choice.MonoExpo_stretched.AICc = AICc_monoS;
    model_choice.MonoExpo_stretched.BIC = BIC_monoS;
    AIC = [ [AIC] model_choice.MonoExpo_stretched.AIC];
    AICc = [ [AICc] model_choice.MonoExpo_stretched.AICc];
    BIC = [ [BIC] model_choice.MonoExpo_stretched.BIC];
    
    clear MLES
    clear binranges_length
    clear AIC_monoS AICc_monoS BIC_monoS
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end


%-----------------
% model = double-expo

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    
    binranges_length = [];
    AIC_double = 0;
    BIC_double = 0;
    AICc_double = 0;
    sum_MLE2 = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            double_ = fitting_results.DoubleExpo.(name_embryo).double;
            MLE2 = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE2 = MLE2 + double_(i)-data(i,2)*log(double_(i)) + cste;
                    %MLE2 = MLE2 + double_(i)-data(i,2)*log(double_(i));
                    MLE2v(i) = double_(i)-data(i,2)*log(double_(i));
                else
                    MLE2 = MLE2 + double_(i) ;
                    MLE2v(i) = double_(i) ;
                end
            end
            
            sum_MLE2 = sum_MLE2 + MLE2;
            model_choice.DoubleExpo.(name_embryo).logMLv = MLE2v;
            model_choice.DoubleExpo.(name_embryo).logML = MLE2;
            model_choice.DoubleExpo.(name_embryo).AIC = 2*MLE2 +2*fitting_results.DoubleExpo.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.DoubleExpo.(name_embryo).AICc = 2*MLE2 +2*fitting_results.DoubleExpo.nparam + 2*fitting_results.DoubleExpo.nparam*(fitting_results.DoubleExpo.nparam+1)/...
                (length(binranges) - fitting_results.DoubleExpo.nparam -1); %wagenmakers & farrell 2004
            model_choice.DoubleExpo.(name_embryo).BIC = 2*MLE2 +fitting_results.DoubleExpo.nparam * log(length(binranges) );
            AIC_double = AIC_double + model_choice.DoubleExpo.(name_embryo).AIC;
            BIC_double = BIC_double + model_choice.DoubleExpo.(name_embryo).BIC;
            AICc_double = AICc_double + model_choice.DoubleExpo.(name_embryo).AICc;
            
            clear double_
            clear binranges
            clear data
            clear cste
            clear Se
            clear MLE2v
        end
        
    end
    
    % commented lines below since I think this is not correct
    %     model_choice.DoubleExpo.AIC = 2*sum_MLE2 +2*fitting_results.DoubleExpo.nparam;
    %     model_choice.DoubleExpo.AICc = 2*sum_MLE2 +2*fitting_results.DoubleExpo.nparam + 2*fitting_results.DoubleExpo.nparam*(fitting_results.DoubleExpo.nparam+1)/...
    %         (mean(binranges_length)  - fitting_results.DoubleExpo.nparam -1); %wagenmakers & farrell 2004
    %     model_choice.DoubleExpo.BIC = 2*sum_MLE2 +fitting_results.DoubleExpo.nparam * log(mean(binranges_length) );
    model_choice.DoubleExpo.logML = sum_MLE2;
    model_choice.DoubleExpo.AIC = AIC_double;
    model_choice.DoubleExpo.AICc = AICc_double;
    model_choice.DoubleExpo.BIC = BIC_double;
    AIC = [ [AIC] model_choice.DoubleExpo.AIC];
    AICc = [ [AICc] model_choice.DoubleExpo.AICc];
    BIC = [ [BIC] model_choice.DoubleExpo.BIC];
    
    clear MLE2
    clear binranges_length
    clear AIC_double AICc_double BIC_double
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end

%-----------------
% model = double-expo with fixed T0

if ~isempty(find(ismember(models,'DoubleExpo_fixedT0'),1))
    
    binranges_length = [];
    AIC_double_ = 0;
    BIC_double_ = 0;
    AICc_double_ = 0;
    sum_MLE2_ = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            double__ = fitting_results.DoubleExpo_fixedT0.(name_embryo).double_;
            MLE2_ = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE2_ = MLE2_ + double__(i)-data(i,2)*log(double__(i)) + cste;
                    %MLE2_ = MLE2_ + double__(i)-data(i,2)*log(double__(i));
                    MLE2v_(i) = double__(i)-data(i,2)*log(double__(i));
                else
                    MLE2_ = MLE2_ + double__(i) ;
                    MLE2v_(i) = double__(i) ;
                end
            end
            
            sum_MLE2_ = sum_MLE2_ + MLE2_;
            model_choice.DoubleExpo_fixedT0.(name_embryo).logMLv = MLE2v_;
            model_choice.DoubleExpo_fixedT0.(name_embryo).logML = MLE2_;
            model_choice.DoubleExpo_fixedT0.(name_embryo).AIC = 2*MLE2_ +2*fitting_results.DoubleExpo_fixedT0.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.DoubleExpo_fixedT0.(name_embryo).AICc = 2*MLE2_ +2*fitting_results.DoubleExpo_fixedT0.nparam + ...
                2*fitting_results.DoubleExpo_fixedT0.nparam*(fitting_results.DoubleExpo_fixedT0.nparam+1)/...
                (length(binranges) - fitting_results.DoubleExpo_fixedT0.nparam -1); %wagenmakers & farrell 2004
            model_choice.DoubleExpo_fixedT0.(name_embryo).BIC = 2*MLE2_ +fitting_results.DoubleExpo_fixedT0.nparam * log(length(binranges) );
            AIC_double_ = AIC_double_ + model_choice.DoubleExpo_fixedT0.(name_embryo).AIC;
            BIC_double_ = BIC_double_ + model_choice.DoubleExpo_fixedT0.(name_embryo).BIC;
            AICc_double_ = AICc_double_ + model_choice.DoubleExpo_fixedT0.(name_embryo).AICc;
            
            clear double__
            clear binranges
            clear data
            clear cste
            clear Se
            clear MLE2v_
        end
        
    end
    
    model_choice.DoubleExpo_fixedT0.logML = sum_MLE2_;
    model_choice.DoubleExpo_fixedT0.AIC = AIC_double_;
    model_choice.DoubleExpo_fixedT0.AICc = AICc_double_;
    model_choice.DoubleExpo_fixedT0.BIC = BIC_double_;
    AIC = [ [AIC] model_choice.DoubleExpo_fixedT0.AIC];
    AICc = [ [AICc] model_choice.DoubleExpo_fixedT0.AICc];
    BIC = [ [BIC] model_choice.DoubleExpo_fixedT0.BIC];
    
    clear MLE2_
    clear binranges_length
    clear AIC_double_ AICc_double_ BIC_double_
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end


%-----------------
% model = double-expo with fixed T0P0

if ~isempty(find(ismember(models,'DoubleExpo_fixedT0P0'),1))
    
    binranges_length = [];
    AIC_double__ = 0;
    BIC_double__ = 0;
    AICc_double__ = 0;
    sum_MLE2__ = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            double___ = fitting_results.DoubleExpo_fixedT0P0.(name_embryo).double__;
            MLE2__ = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE2__ = MLE2__ + double___(i)-data(i,2)*log(double___(i)) + cste;
                    %MLE2__ = MLE2__ + double___(i)-data(i,2)*log(double___(i));
                    MLE2v__(i) = double___(i)-data(i,2)*log(double___(i));
                else
                    MLE2__ = MLE2__ + double___(i) ;
                    MLE2v__(i) = double___(i) ;
                end
            end
            
            sum_MLE2__ = sum_MLE2__ + MLE2__;
            model_choice.DoubleExpo_fixedT0P0.(name_embryo).logMLv = MLE2v__;
            model_choice.DoubleExpo_fixedT0P0.(name_embryo).logML = MLE2__;
            model_choice.DoubleExpo_fixedT0P0.(name_embryo).AIC = 2*MLE2__ +2*fitting_results.DoubleExpo_fixedT0P0.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.DoubleExpo_fixedT0P0.(name_embryo).AICc = 2*MLE2__ +2*fitting_results.DoubleExpo_fixedT0P0.nparam + ...
                2*fitting_results.DoubleExpo_fixedT0P0.nparam*(fitting_results.DoubleExpo_fixedT0P0.nparam+1)/...
                (length(binranges) - fitting_results.DoubleExpo_fixedT0P0.nparam -1); %wagenmakers & farrell 2004
            model_choice.DoubleExpo_fixedT0P0.(name_embryo).BIC = 2*MLE2__ +fitting_results.DoubleExpo_fixedT0P0.nparam * log(length(binranges) );
            AIC_double__ = AIC_double__ + model_choice.DoubleExpo_fixedT0P0.(name_embryo).AIC;
            BIC_double__ = BIC_double__ + model_choice.DoubleExpo_fixedT0P0.(name_embryo).BIC;
            AICc_double__ = AICc_double__ + model_choice.DoubleExpo_fixedT0P0.(name_embryo).AICc;
            
            clear double___
            clear binranges
            clear data
            clear cste
            clear Se
            clear MLE2v_
        end
        
    end
    
    model_choice.DoubleExpo_fixedT0P0.logML = sum_MLE2__;
    model_choice.DoubleExpo_fixedT0P0.AIC = AIC_double__;
    model_choice.DoubleExpo_fixedT0P0.AICc = AICc_double__;
    model_choice.DoubleExpo_fixedT0P0.BIC = BIC_double__;
    AIC = [ [AIC] model_choice.DoubleExpo_fixedT0P0.AIC];
    AICc = [ [AICc] model_choice.DoubleExpo_fixedT0P0.AICc];
    BIC = [ [BIC] model_choice.DoubleExpo_fixedT0P0.BIC];
    
    clear MLE2__
    clear binranges_length
    clear AIC_double__ AICc_double__ BIC_double__
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end

%-----------------
% model = triple-expo

if ~isempty(find(ismember(models,'TripleExpo'),1))
    
    binranges_length = [];
    AIC_triple = 0;
    BIC_triple = 0;
    AICc_triple = 0;
    sum_MLE3 = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            triple_ = fitting_results.TripleExpo.(name_embryo).triple;
            MLE3 = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE3 = MLE3 + triple_(i)-data(i,2)*log(triple_(i)) + cste;
                    %MLE3 = MLE3 + triple_(i)-data(i,2)*log(triple_(i));
                    MLE3v(i) = triple_(i)-data(i,2)*log(triple_(i));
                else
                    MLE3 = MLE3 + triple_(i);
                    MLE3v(i) = triple_(i);
                end
            end
            
            sum_MLE3 = sum_MLE3 + MLE3;
            model_choice.TripleExpo.(name_embryo).logMLv = MLE3v;
            model_choice.TripleExpo.(name_embryo).logML = MLE3;
            model_choice.TripleExpo.(name_embryo).AIC = 2*MLE3 +2*fitting_results.TripleExpo.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.TripleExpo.(name_embryo).AICc = 2*MLE3 +2*fitting_results.TripleExpo.nparam + 2*fitting_results.TripleExpo.nparam*(fitting_results.TripleExpo.nparam+1)/...
                (length(binranges) - fitting_results.TripleExpo.nparam -1); %wagenmakers & farrell 2004
            model_choice.TripleExpo.(name_embryo).BIC = 2*MLE3 +fitting_results.TripleExpo.nparam * log(length(binranges) );
            AIC_triple = AIC_triple + model_choice.TripleExpo.(name_embryo).AIC;
            BIC_triple = BIC_triple + model_choice.TripleExpo.(name_embryo).BIC;
            AICc_triple = AICc_triple + model_choice.TripleExpo.(name_embryo).AICc;
            
            clear binranges
            clear data
            clear triple_
            clear data_norm
            clear cste
            clear Se
            clear MLE3v
        end
        
    end
    
    % commented lines below since I think this is not correct
    %      model_choice.TripleExpo.AIC = 2*sum_MLE3 +2*fitting_results.TripleExpo.nparam;
    %      model_choice.TripleExpo.AICc = 2*sum_MLE3 +2*fitting_results.TripleExpo.nparam + 2*fitting_results.TripleExpo.nparam*(fitting_results.TripleExpo.nparam+1)/...
    %          (mean(binranges_length)  - fitting_results.TripleExpo.nparam -1); %wagenmakers & farrell 2004
    %      model_choice.TripleExpo.BIC = 2*sum_MLE3 +fitting_results.TripleExpo.nparam * log(mean(binranges_length) );
    model_choice.TripleExpo.logML = sum_MLE3;
    model_choice.TripleExpo.AIC = AIC_triple;
    model_choice.TripleExpo.AICc = AICc_triple;
    model_choice.TripleExpo.BIC = BIC_triple;
    AIC = [ [AIC] model_choice.TripleExpo.AIC ];
    AICc = [ [AICc] model_choice.TripleExpo.AICc ];
    BIC = [ [BIC] model_choice.TripleExpo.BIC ];
    
    clear MLE3
    clear binranges_length
    clear AIC_triple AICc_triple BIC_triple
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end

% ----------------------------
% triple expo model with fixed T0

if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1))
    
    binranges_length = [];
    AIC_triple_ = 0;
    BIC_triple_ = 0;
    AICc_triple_ = 0;
    sum_MLE3_ = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            triple__ = fitting_results.TripleExpo_fixedT0.(name_embryo).triple_;
            MLE3_ = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE3_ = MLE3_ + triple__(i)-data(i,2)*log(triple__(i)) + cste;
                    %MLE3 = MLE3 + triple_(i)-data(i,2)*log(triple_(i));
                    MLE3v_(i) = triple__(i)-data(i,2)*log(triple__(i));
                else
                    MLE3_ = MLE3_ + triple__(i);
                    MLE3v_(i) = triple__(i);
                end
            end
            
            sum_MLE3_ = sum_MLE3_ + MLE3_;
            model_choice.TripleExpo_fixedT0.(name_embryo).logMLv = MLE3v_;
            model_choice.TripleExpo_fixedT0.(name_embryo).logML = MLE3_;
            model_choice.TripleExpo_fixedT0.(name_embryo).AIC = 2*MLE3_ +2*fitting_results.TripleExpo_fixedT0.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.TripleExpo_fixedT0.(name_embryo).AICc = 2*MLE3_ +2*fitting_results.TripleExpo_fixedT0.nparam + ...
                2*fitting_results.TripleExpo_fixedT0.nparam*(fitting_results.TripleExpo_fixedT0.nparam+1)/...
                (length(binranges) - fitting_results.TripleExpo_fixedT0.nparam -1); %wagenmakers & farrell 2004
            model_choice.TripleExpo_fixedT0.(name_embryo).BIC = 2*MLE3_ +fitting_results.TripleExpo_fixedT0.nparam * log(length(binranges) );
            AIC_triple_ = AIC_triple_ + model_choice.TripleExpo_fixedT0.(name_embryo).AIC;
            BIC_triple_ = BIC_triple_ + model_choice.TripleExpo_fixedT0.(name_embryo).BIC;
            AICc_triple_ = AICc_triple_ + model_choice.TripleExpo_fixedT0.(name_embryo).AICc;
            
            clear binranges
            clear data
            clear triple__
            clear data_norm
            clear cste
            clear Se
            clear MLE3v_
        end
        
    end
    
    model_choice.TripleExpo_fixedT0.logML = sum_MLE3_;
    model_choice.TripleExpo_fixedT0.AIC = AIC_triple_;
    model_choice.TripleExpo_fixedT0.AICc = AICc_triple_;
    model_choice.TripleExpo_fixedT0.BIC = BIC_triple_;
    AIC = [ [AIC] model_choice.TripleExpo_fixedT0.AIC ];
    AICc = [ [AICc] model_choice.TripleExpo_fixedT0.AICc ];
    BIC = [ [BIC] model_choice.TripleExpo_fixedT0.BIC ];
    
    clear MLE3_
    clear binranges_length
    clear AIC_triple_ AICc_triple_ BIC_triple_
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end


% ----------------------------
% triple expo model with fixed T0 and fixed P0

if ~isempty(find(ismember(models,'TripleExpo_fixedT0P0'),1))
    
    binranges_length = [];
    AIC_triple__ = 0;
    BIC_triple__ = 0;
    AICc_triple__ = 0;
    sum_MLE3__ = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            triple___ = fitting_results.TripleExpo_fixedT0P0.(name_embryo).triple__;
            MLE3__ = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE3__ = MLE3__ + triple___(i)-data(i,2)*log(triple___(i)) + cste;
                    %MLE3 = MLE3 + triple_(i)-data(i,2)*log(triple_(i));
                    MLE3v__(i) = triple___(i)-data(i,2)*log(triple___(i));
                else
                    MLE3__ = MLE3__ + triple___(i);
                    MLE3v__(i) = triple___(i);
                end
            end
            
            sum_MLE3__ = sum_MLE3__ + MLE3__;
            model_choice.TripleExpo_fixedT0P0.(name_embryo).logMLv = MLE3v__;
            model_choice.TripleExpo_fixedT0P0.(name_embryo).logML = MLE3__;
            model_choice.TripleExpo_fixedT0P0.(name_embryo).AIC = 2*MLE3__ +2*fitting_results.TripleExpo_fixedT0P0.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.TripleExpo_fixedT0P0.(name_embryo).AICc = 2*MLE3__ +2*fitting_results.TripleExpo_fixedT0P0.nparam + ...
                2*fitting_results.TripleExpo_fixedT0P0.nparam*(fitting_results.TripleExpo_fixedT0P0.nparam+1)/...
                (length(binranges) - fitting_results.TripleExpo_fixedT0P0.nparam -1); %wagenmakers & farrell 2004
            model_choice.TripleExpo_fixedT0P0.(name_embryo).BIC = 2*MLE3__ +fitting_results.TripleExpo_fixedT0P0.nparam * log(length(binranges) );
            AIC_triple__ = AIC_triple__ + model_choice.TripleExpo_fixedT0P0.(name_embryo).AIC;
            BIC_triple__ = BIC_triple__ + model_choice.TripleExpo_fixedT0P0.(name_embryo).BIC;
            AICc_triple__ = AICc_triple__ + model_choice.TripleExpo_fixedT0P0.(name_embryo).AICc;
            
            clear binranges
            clear data
            clear triple___
            clear data_norm
            clear cste
            clear Se
            clear MLE3v__
        end
        
    end
    
    model_choice.TripleExpo_fixedT0P0.logML = sum_MLE3__;
    model_choice.TripleExpo_fixedT0P0.AIC = AIC_triple__;
    model_choice.TripleExpo_fixedT0P0.AICc = AICc_triple__;
    model_choice.TripleExpo_fixedT0P0.BIC = BIC_triple__;
    AIC = [ [AIC] model_choice.TripleExpo_fixedT0P0.AIC ];
    AICc = [ [AICc] model_choice.TripleExpo_fixedT0P0.AICc ];
    BIC = [ [BIC] model_choice.TripleExpo_fixedT0P0.BIC ];
    
    clear MLE3__
    clear binranges_length
    clear AIC_triple__ AICc_triple__ BIC_triple__
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end


%-----------------
% model = quadro-expo

if ~isempty(find(ismember(models,'QuadroExpo'),1))
    
    binranges_length = [];
    AIC_quadro = 0;
    BIC_quadro = 0;
    AICc_quadro = 0;
    sum_MLE4 = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        if ~isempty(data)
            binranges_length = [ [binranges_length ] length(binranges) ];
            quadro_ = fitting_results.QuadroExpo.(name_embryo).quadro;
            MLE4 = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE4 = MLE4 + quadro_(i)-data(i,2)*log(quadro_(i)) + cste;
                    %MLE4 = MLE4 + quadro_(i)-data(i,2)*log(quadro_(i));
                    MLE4v(i) = quadro_(i)-data(i,2)*log(quadro_(i));
                else
                    MLE4 = MLE4 + quadro_(i);
                    MLE4v(i) = quadro_(i);
                end
            end
            
            sum_MLE4 = sum_MLE4 + MLE4;
            model_choice.QuadroExpo.(name_embryo).logMLv = MLE4v;
            model_choice.QuadroExpo.(name_embryo).logML = MLE4;
            model_choice.QuadroExpo.(name_embryo).AIC = 2*MLE4 +2*fitting_results.QuadroExpo.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.QuadroExpo.(name_embryo).AICc = 2*MLE4 +2*fitting_results.QuadroExpo.nparam + 2*fitting_results.QuadroExpo.nparam*(fitting_results.QuadroExpo.nparam+1)/...
                (length(binranges) - fitting_results.QuadroExpo.nparam -1); %wagenmakers & farrell 2004
            model_choice.QuadroExpo.(name_embryo).BIC = 2*MLE4 +fitting_results.QuadroExpo.nparam * log(length(binranges) );
            AIC_quadro = AIC_quadro + model_choice.QuadroExpo.(name_embryo).AIC;
            BIC_quadro = BIC_quadro + model_choice.QuadroExpo.(name_embryo).BIC;
            AICc_quadro = AICc_quadro + model_choice.QuadroExpo.(name_embryo).AICc;
            
            clear binranges
            clear data
            clear quadro_
            clear data_norm
            clear cste
            clear Se
            clear MLE4v
        end
    end
    
    % commented lines below since I think this is not correct
    %      model_choice.QuadroExpo.AIC = 2*sum_MLE4 +2*fitting_results.QuadroExpo.nparam;
    %      model_choice.QuadroExpo.AICc = 2*sum_MLE4 +2*fitting_results.QuadroExpo.nparam + 2*fitting_results.QuadroExpo.nparam*(fitting_results.QuadroExpo.nparam+1)/...
    %          (mean(binranges_length)  - fitting_results.QuadroExpo.nparam -1); %wagenmakers & farrell 2004
    %      model_choice.QuadroExpo.BIC = 2*sum_MLE4 +fitting_results.QuadroExpo.nparam * log(mean(binranges_length) );
    model_choice.QuadroExpo.logML = sum_MLE4;
    model_choice.QuadroExpo.AIC = AIC_quadro;
    model_choice.QuadroExpo.AICc = AICc_quadro;
    model_choice.QuadroExpo.BIC = BIC_quadro;
    AIC = [ [AIC] model_choice.QuadroExpo.AIC ];
    AICc = [ [AICc] model_choice.QuadroExpo.AICc ];
    BIC = [ [BIC] model_choice.QuadroExpo.BIC ];
    
    clear MLE4
    clear binranges_length
    clear AIC_quadro AICc_quadro BIC_quadro
    
else
    AIC = [ [AIC] NaN ];
    AICc = [ [AICc] NaN ];
    BIC = [ [BIC] NaN ];
end

%--------------------------------------
% drift and diffusion model

if ~isempty(find(ismember(models,'Drift_diffusion'),1))
    
    binranges_length = [];
    AIC_dd = 0;
    BIC_dd = 0;
    AICc_dd = 0;
    sum_MLE = 0;
    
    for iEmbryo = 1 : nbEmbryo_givenCondition
        
        name_embryo = ['embryo' num2str(iEmbryo)];
        binranges = input.(name_embryo).data(:,1);
        data = input.(name_embryo).data;
        
        if ~isempty(data)
            d_d_ = fitting_results.Drift_diffusion.(name_embryo).d_d;
            binranges_length = [ [binranges_length ] length(binranges) ];
            MLE = 0;
            
            for i=1:length(binranges)
                if data(i,2) > 0
                    cste = 0;
                    for j = 1 : data(i,2)
                        cste = cste + log(j);
                    end
                    MLE = MLE + d_d_(i)-data(i,2)*log(d_d_(i)) + cste;
                    %MLE = MLE + d_d_(i)-data(i,2)*log(d_d_(i));
                    MLEv(i) = d_d_(i)-data(i,2)*log(d_d_(i));
                else
                    MLE = MLE + d_d_(i) ;
                    MLEv(i) = d_d_(i) ;
                end
            end
            
            sum_MLE = sum_MLE + MLE;
            model_choice.Drift_diffusion.(name_embryo).logMLv = MLEv;
            model_choice.Drift_diffusion.(name_embryo).logML = MLE;
            model_choice.Drift_diffusion.(name_embryo).AIC = 2*MLE +2*fitting_results.Drift_diffusion.nparam; % AIC = -2 log(MLE) + 2 * nb_param
            model_choice.Drift_diffusion.(name_embryo).AICc = 2*MLE +2*fitting_results.Drift_diffusion.nparam + 2*fitting_results.Drift_diffusion.nparam*(fitting_results.Drift_diffusion.nparam+1)/...
                (length(binranges) - fitting_results.Drift_diffusion.nparam -1); %wagenmakers & farrell 2004
            model_choice.Drift_diffusion.(name_embryo).BIC = 2*MLE +fitting_results.Drift_diffusion.nparam * log(length(binranges) );
            AIC_dd = AIC_dd + model_choice.Drift_diffusion.(name_embryo).AIC;
            BIC_dd = BIC_dd + model_choice.Drift_diffusion.(name_embryo).BIC;
            AICc_dd = AICc_dd + model_choice.Drift_diffusion.(name_embryo).AICc;
            
            clear d_d_
            clear data
            clear binranges
            clear cste
            clear Se
            clear MLEv
        end
        
    end
    
    % commented lines below since I think this is not correct
    %     model_choice.Drift_diffusion.AIC = 2*sum_MLE +2*fitting_results.Drift_diffusion.nparam; % AIC = -2 log(MLE) + 2 * nb_param
    %     model_choice.Drift_diffusion.AICc = 2*sum_MLE +2*fitting_results.Drift_diffusion.nparam + 2*fitting_results.Drift_diffusion.nparam*(fitting_results.Drift_diffusion.nparam+1)/...
    %        (mean(binranges_length) - fitting_results.Drift_diffusion.nparam -1); %wagenmakers & farrell 2004
    %     model_choice.Drift_diffusion.BIC = 2*sum_MLE +fitting_results.Drift_diffusion.nparam * log(mean(binranges_length) );
    model_choice.Drift_diffusion.logML = sum_MLE;
    model_choice.Drift_diffusion.AIC = AIC_dd;
    model_choice.Drift_diffusion.AICc = AICc_dd;
    model_choice.Drift_diffusion.BIC = BIC_dd;
    AIC = [ [AIC] model_choice.Drift_diffusion.AIC];
    AICc = [ [AICc] model_choice.Drift_diffusion.AICc];
    BIC = [ [BIC] model_choice.Drift_diffusion.BIC];
    
    clear MLE
    clear binranges_length
    clear AIC_dd AICc_dd BIC_dd
    
else
    AIC = [ [AIC] NaN];
    AICc = [ [AICc] NaN];
    BIC = [ [BIC]  NaN];
end



%% calculate PrM of indivudal embryo: may be useful to use it later to
% ponderate the influence of a given embryo for the averaging

weight_AIC_sum = 0; %normalization factor for model probabilities
weight_AICc_sum = 0;
weight_BIC_sum = 0;


[AIC_min,Index_min]=min(AIC);
[AICc_min,Index_AICc_min]=min(AICc);
[BIC_min,Index_BIC_min]=min(BIC);

if Index_min == 1
    model_choice.best_model_AIC = 'MonoExpo';
elseif Index_min ==2
    model_choice.best_model_AIC = 'MonoExpo_stretched';
elseif Index_min ==3
    model_choice.best_model_AIC = 'DoubleExpo';
elseif Index_min ==4
    model_choice.best_model_AIC = 'DoubleExpo_fixedT0';
elseif Index_min ==5
    model_choice.best_model_AIC = 'DoubleExpo_fixedT0P0';
elseif Index_min ==6
    model_choice.best_model_AIC = 'TripleExpo';
elseif Index_min ==7
    model_choice.best_model_AIC = 'TripleExpo_fixedT0';
elseif Index_min ==8
    model_choice.best_model_AIC = 'TripleExpo_fixedT0P0';
elseif Index_min ==9
    model_choice.best_model_AIC = 'QuadroExpo';
elseif Index_min ==10
    model_choice.best_model_AIC = 'Drift_diffusion';    
end

if Index_AICc_min == 1
    model_choice.best_model_AICc = 'MonoExpo';
elseif Index_AICc_min ==2
    model_choice.best_model_AICc = 'MonoExpo_stretched';
elseif Index_AICc_min ==3
    model_choice.best_model_AICc = 'DoubleExpo';
elseif Index_AICc_min ==4
    model_choice.best_model_AICc = 'DoubleExpo_fixedT0';
elseif Index_AICc_min ==5
    model_choice.best_model_AICc = 'DoubleExpo_fixedT0P0';
elseif Index_AICc_min ==6
    model_choice.best_model_AICc = 'TripleExpo';
elseif Index_AICc_min ==7
    model_choice.best_model_AICc = 'TripleExpo_fixedT0';
elseif Index_AICc_min ==8
    model_choice.best_model_AICc = 'TripleExpo_fixedT0P0';
elseif Index_AICc_min ==9
    model_choice.best_model_AICc = 'QuadroExpo';
elseif Index_AICc_min ==10
    model_choice.best_model_AICc = 'Drift_diffusion';    
end

if Index_BIC_min == 1
    model_choice.best_model_BIC = 'MonoExpo';
elseif Index_BIC_min ==2
    model_choice.best_model_BIC = 'MonoExpo_stretched';
elseif Index_BIC_min ==3
    model_choice.best_model_BIC = 'DoubleExpo';
elseif Index_BIC_min ==4
    model_choice.best_model_BIC = 'DoubleExpo_fixedT0';
elseif Index_BIC_min ==5
    model_choice.best_model_BIC = 'DoubleExpo_fixedT0P0';
elseif Index_BIC_min ==6
    model_choice.best_model_BIC = 'TripleExpo';
elseif Index_BIC_min ==7
    model_choice.best_model_BIC = 'TripleExpo_fixedT0';
elseif Index_BIC_min ==8
    model_choice.best_model_BIC = 'TripleExpo_fixedT0P0';
elseif Index_BIC_min ==9
    model_choice.best_model_BIC = 'QuadroExpo';
elseif Index_BIC_min ==10
    model_choice.best_model_BIC = 'Drift_diffusion';    
end

% to choose to keep BIC as the best criteria to get the best model

best_model = model_choice.best_model_BIC;


%% Calculate the normalized probablity of each model

if ~isempty(find(ismember(models,'MonoExpo'),1))
    weight_AIC_mono = exp(-0.5*(model_choice.MonoExpo.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_mono;
    weight_AICc_mono = exp(-0.5*(model_choice.MonoExpo.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_mono;
    weight_BIC_mono = exp(-0.5*(model_choice.MonoExpo.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_mono;
    clear AIC_mono
    clear AICc_mono
    clear BIC_mono
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    weight_AIC_double = exp(-0.5*(model_choice.DoubleExpo.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_double;
    weight_AICc_double = exp(-0.5*(model_choice.DoubleExpo.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_double;
    weight_BIC_double = exp(-0.5*(model_choice.DoubleExpo.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_double;
    clear AIC_double
    clear AICc_double
    clear BIC_double
end
if ~isempty(find(ismember(models,'DoubleExpo_fixedT0'),1))
    weight_AIC_double_ = exp(-0.5*(model_choice.DoubleExpo_fixedT0.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_double_;
    weight_AICc_double_ = exp(-0.5*(model_choice.DoubleExpo_fixedT0.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_double_;
    weight_BIC_double_ = exp(-0.5*(model_choice.DoubleExpo_fixedT0.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_double_;
    clear AIC_double_
    clear AICc_double_
    clear BIC_double_
end
if ~isempty(find(ismember(models,'DoubleExpo_fixedT0P0'),1))
    weight_AIC_double__ = exp(-0.5*(model_choice.DoubleExpo_fixedT0P0.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_double__;
    weight_AICc_double__ = exp(-0.5*(model_choice.DoubleExpo_fixedT0P0.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_double__;
    weight_BIC_double__ = exp(-0.5*(model_choice.DoubleExpo_fixedT0P0.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_double__;
    clear AIC_double__
    clear AICc_double__
    clear BIC_double__
end
if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    weight_AIC_monoS = exp(-0.5*(model_choice.MonoExpo_stretched.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_monoS;
    weight_AICc_monoS = exp(-0.5*(model_choice.MonoExpo_stretched.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_monoS;
    weight_BIC_monoS = exp(-0.5*(model_choice.MonoExpo_stretched.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_monoS;
    clear AIC_monoS
    clear AICc_monoS
    clear BIC_monoS
end
if ~isempty(find(ismember(models,'TripleExpo'),1))
    weight_AIC_triple = exp(-0.5*(model_choice.TripleExpo.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_triple;
    weight_AICc_triple = exp(-0.5*(model_choice.TripleExpo.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_triple;
    weight_BIC_triple = exp(-0.5*(model_choice.TripleExpo.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_triple;
    clear AIC_triple
    clear AICc_triple
    clear BIC_triple
end
if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1))
    weight_AIC_triple_ = exp(-0.5*(model_choice.TripleExpo_fixedT0.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_triple_;
    weight_AICc_triple_ = exp(-0.5*(model_choice.TripleExpo_fixedT0.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_triple_;
    weight_BIC_triple_ = exp(-0.5*(model_choice.TripleExpo_fixedT0.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_triple_;
    clear AIC_triple
    clear AICc_triple
    clear BIC_triple
end
if ~isempty(find(ismember(models,'TripleExpo_fixedT0P0'),1))
    weight_AIC_triple__ = exp(-0.5*(model_choice.TripleExpo_fixedT0P0.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_triple__;
    weight_AICc_triple__ = exp(-0.5*(model_choice.TripleExpo_fixedT0P0.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_triple__;
    weight_BIC_triple__ = exp(-0.5*(model_choice.TripleExpo_fixedT0P0.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_triple__;
    clear AIC_triple
    clear AICc_triple
    clear BIC_triple
end
if ~isempty(find(ismember(models,'QuadroExpo'),1))
    weight_AIC_quadro = exp(-0.5*(model_choice.QuadroExpo.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_quadro;
    weight_AICc_quadro = exp(-0.5*(model_choice.QuadroExpo.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_quadro;
    weight_BIC_quadro = exp(-0.5*(model_choice.QuadroExpo.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_quadro;
    clear AIC_quadro
    clear AICc_quadro
    clear BIC_quadro
end
if ~isempty(find(ismember(models,'Drift_diffusion'),1))
    weight_AIC_dd = exp(-0.5*(model_choice.Drift_diffusion.AIC-AIC_min));
    weight_AIC_sum = weight_AIC_sum + weight_AIC_dd;
    weight_AICc_dd = exp(-0.5*(model_choice.Drift_diffusion.AICc-AICc_min));
    weight_AICc_sum = weight_AICc_sum + weight_AICc_dd;
    weight_BIC_dd = exp(-0.5*(model_choice.Drift_diffusion.BIC-BIC_min));
    weight_BIC_sum = weight_BIC_sum + weight_BIC_dd;
    clear AIC_dd
    clear AICc_dd
    clear BIC_dd
end

%-------------------------------
% Calculate the normalized probablity of each model

if ~isempty(find(ismember(models,'MonoExpo'),1))
    model_choice.MonoExpo.PrM_AIC = weight_AIC_mono/weight_AIC_sum;
    model_choice.MonoExpo.PrM_AICc = weight_AICc_mono/weight_AICc_sum;
    model_choice.MonoExpo.PrM_BIC = weight_BIC_mono/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'DoubleExpo'),1))
    model_choice.DoubleExpo.PrM_AIC = weight_AIC_double/weight_AIC_sum;
    model_choice.DoubleExpo.PrM_AICc = weight_AICc_double/weight_AICc_sum;
    model_choice.DoubleExpo.PrM_BIC = weight_BIC_double/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'DoubleExpo_fixedT0'),1))
    model_choice.DoubleExpo_fixedT0.PrM_AIC = weight_AIC_double_/weight_AIC_sum;
    model_choice.DoubleExpo_fixedT0.PrM_AICc = weight_AICc_double_/weight_AICc_sum;
    model_choice.DoubleExpo_fixedT0.PrM_BIC = weight_BIC_double_/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'DoubleExpo_fixedT0P0'),1))
    model_choice.DoubleExpo_fixedT0P0.PrM_AIC = weight_AIC_double__/weight_AIC_sum;
    model_choice.DoubleExpo_fixedT0P0.PrM_AICc = weight_AICc_double__/weight_AICc_sum;
    model_choice.DoubleExpo_fixedT0P0.PrM_BIC = weight_BIC_double__/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    model_choice.MonoExpo_stretched.PrM_AIC = weight_AIC_monoS/weight_AIC_sum;
    model_choice.MonoExpo_stretched.PrM_AICc = weight_AICc_monoS/weight_AICc_sum;
    model_choice.MonoExpo_stretched.PrM_BIC = weight_BIC_monoS/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'TripleExpo'),1))
    model_choice.TripleExpo.PrM_AIC = weight_AIC_triple/weight_AIC_sum;
    model_choice.TripleExpo.PrM_AICc = weight_AICc_triple/weight_AICc_sum;
    model_choice.TripleExpo.PrM_BIC = weight_BIC_triple/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'TripleExpo_fixedT0'),1))
    model_choice.TripleExpo_fixedT0.PrM_AIC = weight_AIC_triple_/weight_AIC_sum;
    model_choice.TripleExpo_fixedT0.PrM_AICc = weight_AICc_triple_/weight_AICc_sum;
    model_choice.TripleExpo_fixedT0.PrM_BIC = weight_BIC_triple_/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'TripleExpo_fixedT0P0'),1))
    model_choice.TripleExpo_fixedT0P0.PrM_AIC = weight_AIC_triple__/weight_AIC_sum;
    model_choice.TripleExpo_fixedT0P0.PrM_AICc = weight_AICc_triple__/weight_AICc_sum;
    model_choice.TripleExpo_fixedT0P0.PrM_BIC = weight_BIC_triple__/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'QuadroExpo'),1))
    model_choice.QuadroExpo.PrM_AIC = weight_AIC_quadro/weight_AIC_sum;
    model_choice.QuadroExpo.PrM_AICc = weight_AICc_quadro/weight_AICc_sum;
    model_choice.QuadroExpo.PrM_BIC = weight_BIC_quadro/weight_BIC_sum;
end
if ~isempty(find(ismember(models,'Drift_diffusion'),1))
    model_choice.Drift_diffusion.PrM_AIC = weight_AIC_dd/weight_AIC_sum;
    model_choice.Drift_diffusion.PrM_AICc = weight_AICc_dd/weight_AICc_sum;
    model_choice.Drift_diffusion.PrM_BIC = weight_BIC_dd/weight_BIC_sum;
end


end

