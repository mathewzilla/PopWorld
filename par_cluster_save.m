function par_cluster_save(fname,Full,Connected,clusterpars)
% dirty and inflexible function for saving clustering output within a parfor loop

save(fname,'Full','Connected','clusterpars','-v7.3')