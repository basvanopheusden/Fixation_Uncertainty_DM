
clear
FMod='Negstd2';
nrun=1;

Budget = 4; %on average, 4 repitition per trial
nrep={};
for whichSubj = 1:39
    
    for ndT=0:4
        try
            savefile=['~/Google Drive/Zhiwei-project-mats/' FMod '/trialLL/' FMod '_estimLL_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'_run',num2str(nrun)];
            load(savefile)
        catch
            disp(savefile)
        end
        stdcpr=[mean(std(mean(Probs,1),0,3)),std(Probs(:))]; % compare std over diff parameters, or over all trials. std between trial should be larger than within
        trProbs=mean(mean(exp(-Probs),1),3);
        
        reps=sqrt(trProbs.*li2(trProbs));
        nrep{whichSubj,ndT+1}=ceil(reps/sum(reps)*Budget*length(reps));
    end
end
%%
% quickly check the consistency
% pick 2 random pars and compare with 2 random subs
within=[];
between=[];
for whichSubj=1:38
    for ndT=1:4
        try
        incpr=sqrt(mean((nrep{whichSubj,ndT}-nrep{whichSubj,ndT+1}).^2));
        catch
            continue
        end
        if ~isnan(incpr)
            within=[within,incpr];
        end
        
        l1=length(nrep{whichSubj,ndT});
        l2=length(nrep{whichSubj+1,ndT});
        if l1<l2
            betcpr=sqrt(mean((nrep{whichSubj,ndT}-nrep{whichSubj+1,ndT}(1:l1)).^2));
        else
            betcpr=sqrt(mean((nrep{whichSubj,ndT}(1:l2)-nrep{whichSubj+1,ndT}).^2));
        end
        
        if ~isnan(betcpr)
            between=[between,betcpr];
        end
    end
end

