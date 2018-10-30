function par_cluster_save2(fname,Rejection,Data,pars,optionsModel,optionsReject)
% dirty and inflexible function for saving clustering output within a parfor loop

save(fname,'Rejection','pars','optionsModel','optionsReject','Data','-v7.3')