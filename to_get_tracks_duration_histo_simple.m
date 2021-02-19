function [ tracks_duration_histo ] = to_get_tracks_duration_histo_simple( dataTracks_input,iBootstrap,nbEmbryo,tracks_duration_histo,timeDuration, area)

global param

maxLength= max(dataTracks_input);
minLength= min(dataTracks_input);
binranges = [minLength : 1 : maxLength];

if maxLength > 100 % 10 sec
    index_to_remove = [(100-minLength+1):(maxLength-minLength+1)];
    binranges(index_to_remove) = [];
end

[bincounts] = hist(dataTracks_input,binranges);

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
if length(index_to_remove2) <= length(bincounts) && length(index_to_remove2) <= length(binranges)
    binranges(index_to_remove2) = [];
    bincounts(index_to_remove2) = [];
end

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
if bincounts_above_max_one >= 1
    bincounts_above_max_one = bincounts_above_max_one / (max_max_index_one - max_index_one - 1);
end

if ~isnan(bincounts_above_max_one) && bincounts_above_max_one > 0
    index_to_remove3 = [(max_index_one):(max_max_index_one)];
    binranges(index_to_remove3) = [];
    bincounts(index_to_remove3) = [];
    binranges(max_index_one) = max_index_one + ( max_max_index_one - max_index_one)/2;
    bincounts(max_index_one) = bincounts_above_max_one;
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

tracks_duration_histo{iBootstrap,nbEmbryo}.bincounts = bincounts;
tracks_duration_histo{iBootstrap,nbEmbryo}.binranges = binranges./param.sp6; % in sec
tracks_duration_histo{iBootstrap,nbEmbryo}.number_tracks = sum(bincounts);
tracks_duration_histo{iBootstrap,nbEmbryo}.lengthTracks = dataTracks_input;
tracks_duration_histo{iBootstrap,nbEmbryo}.timeDuration_phase = timeDuration;
tracks_duration_histo{iBootstrap,nbEmbryo}.area = area;

clear maxLength
clear minimLength
clear bincounts
clear binranges
clear timeDuration
clear area


end

