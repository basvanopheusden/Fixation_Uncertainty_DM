clear
neval=NaN(39,1);
finalstd=NaN(39,1);
tictoc=NaN(39,1);
FMod='Negstd2';
for whichSubj=1:39
    for nrun=1:5
        for ndT=0:4
            clear LogProb100ms ProbFluct funtable timeRec
            %load(['~/Google Drive/Zhiwei-project-mats/' FMod '/ibs_constrep_BADS/' FMod '_repibs_BADS_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'_run',num2str(nrun)])
            try
                load(['~/Google Drive/Zhiwei-project-mats/' FMod '/ibs_constrep_BADS/' FMod '_repibs_BADS_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'_run',num2str(nrun)])
            catch
                disp(whichSubj)
            end
            if exist('timeRec','var')
                tictoc(whichSubj)=mean(timeRec(:,1));
                if exist('LogProb100ms','var')
                    neval(whichSubj)=length(funtable);
                    finalstd(whichSubj)=std(ProbFluct);
                end
            else
                continue
                
            end
        end
    end
end
%%
[median(tictoc(~isnan(tictoc))),max(tictoc(~isnan(tictoc))),min(tictoc(~isnan(tictoc)))]
[median(finalstd(~isnan(neval))),min(finalstd(~isnan(neval))),max(finalstd(~isnan(neval)))]
[median(neval(~isnan(neval))),min(neval(~isnan(neval))),max(neval(~isnan(neval)))]