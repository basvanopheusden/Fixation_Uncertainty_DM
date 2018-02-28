clear
%%
% question: for those subs with large LL diff, how many "hit" can you get
% per trial?
subj=1;
load(sprintf('Negstd2_fitreal_bads_subj_%d',subj),'thisFittedPara','ndT')


load('FixNumLNR100_fromzero')
load ProcessedData
D = ProcessedData;
AllSubjLabels = unique(D(:,13));
TrialLabels =find(D(:,13)==AllSubjLabels(subj));
SubFixNumLNR = FixNumLNR(TrialLabels);
SubLRating = D(TrialLabels,2);
SubRRating = D(TrialLabels,1);
SubRT = allRT(TrialLabels)-ndT;
SubChoice = D(TrialLabels,3);

ScalingFactor = ones(size(thisFittedPara));
hit = Negstd2_runonce(thisFittedPara,ScalingFactor,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT); 

