function [ output_args ] = ibs_eachtrial_LPwrapper(Para,Func,ScalingFactor, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin,LLcutoff,Nrep,savefile)
% to be called by optimization function, to record every step of calling
lambda = Para(end)/ScalingFactor(end); % used by ibs.
modelfun = str2func(Func);
persistent funtable timeRec
if isempty(Para)
    output_args = funtable;
    funtable = []; % so next time calling you'll get a new table
else
    allLs = zeros(Nrep,length(LRating));
    for krep = 1:Nrep
        tic
        L = zeros(1,length(LRating));
        n = ones(1,length(LRating));
        feature accel on
        
        Tarr = 1:225;
        cutoff = LLcutoff;%4.6; % -log(1/100) = 4.605, 100 for all possible responses
        
        trials_completed = false(1,length(LRating));
        while ~all(trials_completed) && mean(L)<cutoff
            for trial = find(~trials_completed)
                thisFixNums = FixNumLNR{trial};
                thisLR = LRating(trial);
                thisRR = RRating(trial);
                if rand()<lambda
                    ishit = rand()<1/(2*length(RTbin));% Maybe replace with smarter lapse
                else
                    [ifDecided,DeltaUMat] = modelfun(Para, ScalingFactor, thisFixNums, thisLR, thisRR);
                    if sum(ifDecided)==0
                        ishit = false;
                    else
                        RT = Tarr(ifDecided);
                        RT = RT(1);
                        Choi = DeltaUMat(RT)>0;
                        
                        RTinb = find(RT<RTbin(2:end));
                        RTinb = RTinb(1);
                        trRTinb = find(ReactionTime(trial)<RTbin(2:end));
                        trRTinb = trRTinb(1);
                        ishit = (RTinb==trRTinb && Choi == Choice(trial));
                    end
                end
                if ishit
                    trials_completed(trial)=true;
                else
                    L(trial)=L(trial)+1/n(trial);
                    n(trial)=n(trial)+1;
                end
            end
            
        end
        allLs(krep,:)=L;

    end
    output_args = allLs;
            
    funtable = [funtable;Para,-sum(mean(allLs,1)),mean(std(allLs,1))];
    Tcalc = toc;
    timeRec = [timeRec;Tcalc,mean(n),max(n)];
    if ~isempty(savefile)
        save(savefile,'funtable','timeRec')
    end
end

end

