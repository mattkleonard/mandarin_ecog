function badTrials = Find_bad_trials(datarootdir,subj,current_block,stimOnsetTime,timeLim)
% written by Matt
% Adapted by Max on 9/10

bad = load([datarootdir '/' subj '/' subj '_' current_block '/Artifacts/badTimeSegments.mat']);
% load matrix with bad time segments raws are different time range with
% artifacts, first column is onset, second column is offset

badTrials = [];

for j = 1:size(bad.badTimeSegments,1)% number of bad time segments
    badRange = [bad.badTimeSegments(j,1) + timeLim(1) ...
        bad.badTimeSegments(j,2) + timeLim(2)];
    if any(stimOnsetTime >= badRange(1) & stimOnsetTime <= badRange(2)) % if the Stimulus presentation is within a bad time range
        badTrials = [badTrials find(stimOnsetTime >= badRange(1) & stimOnsetTime <= badRange(2))];
    end
end
end