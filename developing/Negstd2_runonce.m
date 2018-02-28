function hit = Negstd2_runonce(Para4,ScalingFactor, FixNumLNR, LRating, RRating, Choice,ReactionTime)

ObsVar = Para4(1)/ScalingFactor(1);
A = Para4(2)/ScalingFactor(2);
B = Para4(3)/ScalingFactor(3);
PriorVar = Para4(4)/ScalingFactor(4);
Step =Para4(5)/ScalingFactor(5);
k =Para4(6)/ScalingFactor(6);
BoundaryFunc = @(t) B*exp(-(t/Step)^k);

feature accel on
hit = NaN(size(LRating));
nrepeat = 1;
%profile on
for trial = 1:length(LRating)
    AllRT = NaN(1,nrepeat);
    AllChoice = NaN(1,nrepeat);
    thisFixNums = FixNumLNR{trial};
    if isempty(thisFixNums)
        continue
    else
        FixNumsL = thisFixNums(1,:);
        FixNumsR = thisFixNums(2,:);
        LVarTerm =  1./(1/PriorVar + FixNumsL/ObsVar);
        RVarTerm =  1./(1/PriorVar + FixNumsR/ObsVar);
        FixSeries = FixNumsL - [0,FixNumsL(1:end-1)];
        LMeanTerm = cumsum(bsxfun(@times, FixSeries,bsxfun(@plus, LRating(trial),sqrt(ObsVar)* randn(1, length(FixSeries)))),2);
        LMeanTerm = bsxfun(@times, 1./(1/PriorVar + FixNumsL/ObsVar), LMeanTerm) / ObsVar;
        RMeanTerm = cumsum(bsxfun(@times, ~FixSeries,bsxfun(@plus, RRating(trial), sqrt(ObsVar)* randn(1, length(FixSeries)))),2);
        RMeanTerm = bsxfun(@times, 1./(1/PriorVar + FixNumsR/ObsVar), RMeanTerm) / ObsVar;
        LUMat = bsxfun(@minus, LMeanTerm, A*sqrt(LVarTerm) ); % neg std model
        RUMat = bsxfun(@minus, RMeanTerm, A*sqrt(RVarTerm) );
        DeltaUMat = LUMat - RUMat; % comparison exp!
        BoundSeries = arrayfun(BoundaryFunc,(1:length(FixSeries)));
        MakeDecision = bsxfun(@minus,abs(DeltaUMat),BoundSeries); % >0 means decision has been made
        
        ifDecided = MakeDecision>0;
        ifDecided = cumsum(ifDecided,2)>0; % make sure once decision is made, all following moments are 1.
        FirstDecision =  ifDecided-[zeros(nrepeat,1),ifDecided(:,1:end-1)];
        [HaveDecided,RT] = find(FirstDecision);
        DecisionMoments = sub2ind(size(MakeDecision), HaveDecided, RT);
        AllRT(HaveDecided) = RT;        
        allDecision = DeltaUMat(DecisionMoments)>0; % >0 = choose left
        AllChoice(HaveDecided) = allDecision;
        hit(trial) = (AllRT == ReactionTime(trial)) && (AllChoice == Choice(trial));
                    
    end
end



end