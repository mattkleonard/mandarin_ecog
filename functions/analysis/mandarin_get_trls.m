function [trls,plot_trls] = mandarin_get_trls(dat,trls_org,ripple_flag,compare_blocks_flag,blocks,syls,spkrs,tones,acc)

trls = [];
if isempty(blocks)
    blocks = unique(dat.behav.block);
end
if compare_blocks_flag
    trls{1} = find((ismember(dat.behav.block,accBlock.good{tones}))' & ...
        (ismember(dat.behav.tone,tones)) & ...
        (ismember(dat.behav.accuracy,acc)));
    trls{2} = find((ismember(dat.behav.block,accBlock.bad{tones}))' & ...
        (ismember(dat.behav.tone,tones)) & ...
        (ismember(dat.behav.accuracy,acc)));
elseif ~ripple_flag
    if strcmpi(trls_org,'tones')
        for i = 1:length(tones)
            trls{i} = find((ismember(dat.behav.block,blocks)) & ...
                (ismember(dat.behav.tone,tones(i))) & ...
                (ismember(dat.behav.accuracy,acc)));
        end
    elseif strcmpi(trls_org,'syls')
        for i = 1:length(syls)
            trls{i} = find((ismember(dat.behav.block,blocks)) & ...
                (strcmpi(dat.behav.syllable,syls{i})) & ...
                (ismember(dat.behav.accuracy,acc)));
        end
    end
else
    for i = 1:length(tones)
        trls{i} = find((ismember(dat.behav.tone,tones(i))) & ...
            (ismember(dat.behav.accuracy,acc)));
    end
end

plot_trls = [];
for i = 1:length(trls)
    plot_trls = [plot_trls ; trls{i}];
end

