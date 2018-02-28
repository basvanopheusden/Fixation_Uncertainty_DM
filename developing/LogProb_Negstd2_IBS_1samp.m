function [SumLogProb_ibs] = LogProb_Negstd2_IBS_1samp(Para4,ScalingFactor, FixNumLNR, LRating, RRating, Choice,ReactionTime)

ObsVar = Para4(1)/ScalingFactor(1);
A = Para4(2)/ScalingFactor(2);
B = Para4(3)/ScalingFactor(3);
PriorVar = Para4(4)/ScalingFactor(4);
Step =Para4(5)/ScalingFactor(5);
k =Para4(6)/ScalingFactor(6);
BoundaryFunc = @(t) B*exp(-(t/Step)^k);

LogProb_ub = NaN(1,length(LRating));


L=[0 cumsum(1./(1:1e5))]; % I'm still adding a cut-off here
feature accel on

nhit = NaN(length(LRating),1);
nrerun = zeros(size(LRating)); % record to see how often do you need rerun
nrepeat = 1;
%profile on
for trial = 1:length(LRating)
    AllRT = -ones(1,nrepeat);
    AllChoice = -ones(1,nrepeat);
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
        temp = 1:length(AllRT);
        hitTime = temp(logical((abs(AllRT- ReactionTime(trial)) < 1).*(AllChoice == Choice(trial))));
        
        if numel(hitTime)==0 % I'm still adding a cut-off here.
            LogProb_ub(trial) = L(end);
            
        elseif numel(hitTime)>1
            hitTime = [hitTime(1), hitTime(2:end)-hitTime(1:end-1)];            
            LogProb_ub(trial) = mean(L(hitTime));                       

        elseif numel(hitTime)==1
            LogProb_ub(trial)=L(hitTime);
        end
                   
        
        
    end
end

SumLogProb_ibs = sum(LogProb_ub);


end