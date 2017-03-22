# mandarin_ecog
Scripts/functions for preprocessing and analyzing ECoG and behavioral data from Mandarin tone learning tasks.

Step 1: Create dat.mat file (e.g., EC131_dat.mat). This file contains the following data:

subj: 'EC143'
hemi: 'lh'
behav:
    block: [216×1 double]
    trial: [216×1 double]
    syllable: {216×1 cell}
    speaker: {216×1 cell}
    tone: [216×1 double]
    choice: [216×1 double]
    rt: [216×1 double]
    accuracy: [216×1 double]
hg_stim: [chan x time x trials]
hg_button: [chan x time x trials]
hg_feedback: [chan x time x trials]
gridOrient: [16 x 16]
badChans: [nChans]

The dat.mat file is created using mandarin_preprocess_data_v2.m. The script finds events (using circular convolution between ANIN and stim.wav files), segments and z-scores data, and puts it into the dat.mat structure.

