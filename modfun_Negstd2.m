function [ifDecided,DeltaUMat] = modfun_Negstd2(Para, ScalingFactor, thisFixNums, LRating, RRating)
% binRT: reduce all RT into K bins; ro: rich output
% adapted from IBS_richoutput.m.
% difference: 1) bin RT (as one more input variable) 2) the "hit"
% comparison line 3) do the real ibs

ObsVar = Para(1)/ScalingFactor(1);
A = Para(2)/ScalingFactor(2); % risk aversion factor
B = Para(3)/ScalingFactor(3);
PriorVar = Para(4)/ScalingFactor(4); %Integer that does not need to scale
Step =Para(5)/ScalingFactor(5);
k =Para(6)/ScalingFactor(6);
BoundaryFunc = @(t) B*exp(-(t/Step)^k);
BoundSeries = arrayfun(BoundaryFunc,(1:225));

FixNumsL = thisFixNums(1,:);
FixNumsR = thisFixNums(2,:);
LVarTerm =  1./(1/PriorVar + FixNumsL/ObsVar);
RVarTerm =  1./(1/PriorVar + FixNumsR/ObsVar);
FixSeries = FixNumsL - [0,FixNumsL(1:end-1)];
LMeanTerm = cumsum(FixSeries .* (LRating+sqrt(ObsVar)* randn(1, length(FixSeries))),2);
LMeanTerm = LMeanTerm./(1/PriorVar + FixNumsL/ObsVar) / ObsVar;
RMeanTerm = cumsum(~FixSeries .*(RRating+sqrt(ObsVar)* randn(1, length(FixSeries))),2);
RMeanTerm = RMeanTerm./(1/PriorVar + FixNumsR/ObsVar)/ ObsVar;
LUMat = LMeanTerm- A*sqrt(LVarTerm); % neg std model
RUMat = RMeanTerm- A*sqrt(RVarTerm);
DeltaUMat = LUMat - RUMat; % comparison exp!

MakeDecision = abs(DeltaUMat)-BoundSeries(1:length(FixSeries)); % >0 means decision has been made

ifDecided = MakeDecision>0;
end
