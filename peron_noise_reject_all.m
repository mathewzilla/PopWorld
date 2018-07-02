% peron_noise_reject
%
% Script to apply null model based Noise Rejection to data from Peron et al 2015

%% Big loop to perform Data_Noise_Rejection on all subvolumes in Peron 2015 (Calcium to start. Events now too)


clear all;

% analysis parameters from Noise_Rejection example repo
pars.N = 100;           % repeats of permutation
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval
pars.Model = 'Poiss';   % Poiss or 'WCM' . % which null model
pars.C = 100;             % conversion factor for real-valued weights (set=1 for integers, 'all' to use full data range)
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble?
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default
optionsReject.Interval = 'CI';

% topdir = ['/Volumes/05/Peron_2015/Peron_ssc-2/ssc-2'];
% topdir = ['/media/mathew/Data_1/Peron_ssc-2/ssc-2'];
topdir = ['/media/mathew/Data_1/Peron_ssc-2/with_events'];
% topdir = ['/mnt/isilon/Hlab/Mat_Evans/Peron_ssc_events'];
animals = {'an171923';'an194181';'an194672';'an197522';'an198503';'an229716';'an229717';'an229719'};
%%
for i = 1:numel(animals)
    %     cd([topdir,'/',animals{i}]);
    files = dir([topdir,'/',animals{i}]);
    for j = 3:numel(files)
        
        load([topdir,'/',animals{i},'/',files(j).name]);
        % Special case for animal 4 
        if i == 4
            s = dat;
        end
        
        % Event traces are always the last (N-1)/2 entries
        N = numel(s.timeSeriesArrayHash.value);
        ev_files = ((N-1)/2)+2 : N;
        
        %         for k = 2:numel(s.timeSeriesArrayHash.value)
        for k = 1:numel(ev_files)
            
            %             fname = [files(j).name(1:end-9),'_sv_',num2str(k-1)]
            % Special case for animal 4
            if i == 4
                fname = [files(j).name(1:end-4),'_events_sv_',num2str(k)]
            else
                fname = [files(j).name(1:end-9),'_events_sv_',num2str(k)]
            end
            
            %             data = s.timeSeriesArrayHash.value{k}.valueMatrix;
            data = s.timeSeriesArrayHash.value{ev_files(k)}.valueMatrix;
            % clean up
            data(find(isnan(data))) = 0;
            
            % Create CXY
            A = corrcoef(data');
            
            % clean up nans again
            A(find(isnan(A))) = 0;
            
            % make undirected if necessary
            A = (A + A') / 2; % make undirected
            
            % Restrict to positive (for now)
            A(find(A<0)) = 0;
            
            % clean-up A, get largest component, and store as basis for all further analysis
            % all indices are with reference to Data.A
            Data = {};
            [Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);
            
            % get expected distribution of eigenvalues under null model
            [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = RndPoissonConfigModel(Data.A,pars.N,pars.C,optionsModel);
            
            %% decompose nodes into signal and noise
            B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
            
            % find low-dimensional projection
            [Data.Dspace,Data.ixpos,Data.Dn,Data.EigEst,Data.Nspace,Data.ixneg,Data.Dneg,Data.NEigEst] = LowDSpace(B,Data.Emodel,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors
            
            % compute dimensions based on just positive eigenvalues
            egs = eig(B);  % eigenspectra of data modularity matrix
            egs = sort(egs,'descend'); % sort eigenvalues into descending order
            Data.PosDn = sum(egs > pars.eg_min);
            
            % node rejection within low-dimensional projection
            Rejection = NodeRejection(B,Data.Emodel,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections
            
            % new signal matrix
            Data.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);
            
            % connected signal matrix: find largest component, and use that - store
            % others
            [Data.Asignal_comp,ixRetain,Data.SignalComps,Data.SignalComp_sizes] = prep_A(Data.Asignal);
            Data.ixSignal_comp = Rejection.ixSignal(ixRetain);  % original node indices
            
            % and then strip out leaves - nodes with single links
            K = sum(Data.Asignal_comp);
            ixLeaves = find(K==1); ixKeep = find(K > 1);
            
            Data.ixSignal_Final = Data.ixSignal_comp(ixKeep);
            Data.ixSignal_Leaves = Data.ixSignal_comp(ixLeaves);
            Data.Asignal_final = Data.Asignal_comp(ixKeep,ixKeep);
            
            %% compare to standard configuration model
            [Control.Emodel,diagnostics,Vmodel] = RndPoissonConfigModel(Data.A,pars.N,pars.C);
            Control.P = expectedA(Data.A);
            
            B = Data.A - Control.P;
            
            % compute groups based on just positive eigenvalues
            egs = eig(B);  % eigenspectra of data modularity matrix
            egs = sort(egs,'descend'); % sort eigenvalues into descending order
            Control.PosDn = sum(egs > pars.eg_min);
            
            % compute groups based on estimated bounds
            [Control.Dspace,~,Control.Dn,Control.EigEst] = LowDSpace(B,Control.Emodel,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors
            
            %% save
            save(['Results_batch1/Rejected_' fname],'Rejection','Data','Control','pars','optionsModel','optionsReject')
        end
    end
end

%% Summarise rejection results into a table
fnames = dir('Results_batch1/');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    
    if any(strfind(fnames(iF).name,'Rejected'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(10:end-4); % strip out 'Rejected' and .mat
        Data = {};
        load(['Results_batch1/' fnames(iF).name]);
        result(netCtr).Network_Size = numel(Data.ixRetain);
        result(netCtr).Signal_Size_WCM = numel(Data.ixSignal_Final);
        result(netCtr).WCM_Dn = Data.PosDn;
        result(netCtr).WCM_RejectionDn = Data.Dn;
        result(netCtr).Config_Dn = Control.PosDn;
        result(netCtr).Config_RejectionDn = Control.Dn;
        result(netCtr).Signal_Components = numel(Data.SignalComp_sizes);
    end
end

Network_Rejection_Table = struct2table(result);
save('Results_batch1/Network_Rejection_Table2','Network_Rejection_Table');

%% Cluster noise rejection results
clear all
blnLabels = 0;      % write node labels? Omit for large networks
fontsize = 6;

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;

files = dir('Results_batch1/Rejected_*');
% files = dir('Results/Rejected_*');

parfor i = 383:566 %numel(files)
    fname = files(i).name(10:end-4)
    
    % load data
    temp_data = load(['Results_batch1/Rejected_', fname,'.mat']);
    Data = temp_data.Data;
    
%     load(['Results/Rejected_', fname,'.mat'])
    
    %% cluster - with noise rejection
    % construct new null model
    P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
    

    
    % then cluster
    Connected = {};
    if Data.Dn > 0
        [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
            ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps);
    else
        Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
    end
    
    % Louvain algorithm
    if Data.Dn > 0
    [Connected.LouvCluster,Connected.LouvQ,~,~] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    else
        Connected.LouvCluster = []; Connected.LouvQ = 0;
    end
    %% cluster - without noise rejection
    Full = {};
    if Data.Dn > 0
        [QmaxCluster,Qmax,ConsCluster, ConsQ,~] = ...
            ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn);
        Full.QmaxCluster = QmaxCluster; Full.Qmax = Qmax; Full.ConsCluster = ConsCluster; Full.ConsQ = ConsQ; 

    else
        Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
    end
    
    [LouvCluster,LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    Full.LouvCluster = LouvCluster; Full.LouvQ = LouvQ;
    %% Save
    par_cluster_save(['Results_batch1/Clustered_',fname,'.mat'],Full,Connected,clusterpars)
%     save(['Results_batch1/Clustered_' fname],'Full','Connected','clusterpars')
end

%% Cluster using output of Control noise rejection
clear all
blnLabels = 0;      % write node labels? Omit for large networks
fontsize = 6;

clusterpars.nreps = 100;
clusterpars.nLouvain = 5;

files = dir('Results/Rejected_*'); %dir('Results_batch1/Rejected_*');

for i = [1,3]; %1:numel(files);
    fname = files(i).name(10:end-4);
    
    % load data
%     load(['Results_batch1/Rejected_', fname,'.mat']
    temp_data = load(['Results/Rejected_', fname,'.mat'])
    Data = temp_data.Data;
    
    %% cluster - with noise rejection
    % construct new null model
%     P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
    
    % Null model from Control rejection
    P = Control.P;
    
    % then cluster
    Connected = {};
    if Data.Dn > 0
        [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
            ConsensusCommunityDetect(Data.A,P,1+Control.Dn,1+Control.Dn,clusterpars.nreps);
    else
        Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
    end
    % Louvain algorithm
    [Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    
    %% cluster - without noise rejection
    Full = {};
    if Data.Dn > 0
        [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
            ConsensusCommunityDetect(Data.A,Data.ExpA,1+Control.Dn,1+Control.Dn);
    else
        Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
    end
    
    [Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    
    %% Save
    save(['Results/Clustered_Config_' fname],'Full','Connected','clusterpars')
end

%% Load Rejection table and plot various features
load('Results_batch1/Network_Rejection_Table2')

%% Identify ca vs events data

%% Re-generate NRstats equivalent and match up SVs
% ALSO computing eigenvalues
topdir = ['/media/mathew/Data_1/Peron_ssc-2/with_events'];
% topdir = ['/mnt/isilon/Hlab/Mat_Evans/Peron_ssc_events'];
animals = {'an171923';'an194181';'an194672';'an197522';'an198503';'an229716';'an229717';'an229719'};

for i = 1:numel(animals)
    %     cd([topdir,'/',animals{i}]);
    files = dir([topdir,'/',animals{i}]);
    for j = 3:numel(files)
        
        load([topdir,'/',animals{i},'/',files(j).name]);
        % Special case for animal 4
        if i == 4
            s = dat;
        end
        
        % Event traces are always the last (N-1)/2 entries
        % Data traces are entry 2 : end of first half 
        N = numel(s.timeSeriesArrayHash.value);
        ev_files = ((N-1)/2)+2 : N;
        dat_files = 2 : ((N-1)/2) + 1;
        
        %         for k = 2:numel(s.timeSeriesArrayHash.value)
        for k = 1:numel(dat_files) %ev_files)
            
            %             fname = [files(j).name(1:end-9),'_sv_',num2str(k-1)]
            % Special case for animal 4
            if i == 4
%                 fname = [files(j).name(1:end-4),'_events_sv_',num2str(k)]
                fname = [files(j).name(1:end-4),'_data_s_sv_',num2str(k)]
            else
%                 fname = [files(j).name(1:end-9),'_events_sv_',num2str(k)]
                fname = [files(j).name(1:end-9),'_data_s_sv_',num2str(k)]
            end
            
            clear statsTable
            statsTable.fname = fname;                                      % filename
            statsTable.method = 'calcium'; %'Peron';                                   % Method
            
            statsTable.Animal = i;                                         % Animal
            statsTable.Session = j-2;                                      % Session
            statsTable.Subvolume = k;                                      % Subvolume
            
            data = s.timeSeriesArrayHash.value{dat_files(k)}.valueMatrix;
            
            [N,T] = size(data);
            
            statsTable.N = N;                                           % Original N cells
            statsTable.T = T;                                           % Duration of recording
            
            % Performance
            trials = s.timeSeriesArrayHash.descrHash{dat_files(k)}.value{1,2};
            trial_data = s.trialTypeMat(:,find(ismember(s.trialIds,trials)));
            
            % P (correct)
            statsTable.Pcorrect = (sum(trial_data(1,:)) + sum(trial_data(2,:))) / length(trial_data);
            
            % P (no lick) - measure of task engagement
            statsTable.Pnolick = (sum(trial_data(5,:)) + sum(trial_data(6,:))) / length(trial_data);
            
            % eigs
            data(find(isnan(data))) = 0;
            [V,D] = eig(cov(data'));
            egs = sort(diag(D),1,'descend');
            explained = 100*cumsum(egs)/sum(egs);
            
            statsTable.egs = egs;
            statsTable.explained = explained;
            
            
            
            % Save these basic stats in a format that is easy to load
            % alongside noise rejection/clustering data
            save(['Results_batch1/StatsTable_',fname,'.mat'],'statsTable');
            
        end
    end
end

%% Summarise stats and rejection results into table (but not eigs - just D to 90% variance) 
fnames = dir('Results_batch1/Rejected_*');
stat_fnames = dir('Results_batch1/StatsTable_*');
nF = numel(fnames);
nS = numel(stat_fnames);
if nF ~= nS
    display('Uh Oh: nF and nS do not match')
end

clear result
netCtr = 0;
for iF = 1:nF

    iS = 0;
    % Check if there's a matching stat file (they should all be present
    for j = 1:nS
        if strcmp(stat_fnames(j).name(12:end-4),fnames(iF).name(10:end-4))
            iS = j;
        end
    end
    
    if iS == 0
        display(['stat file ',num2str(iF),' not found']);
    else
        Data = {};
        load(['Results_batch1/',fnames(iF).name])
        
        result(iF).NetworkName = fnames(iF).name(10:end-4); % strip out 'Rejected' and .mat
        
        result(iF).Network_Size = numel(Data.ixRetain);
        result(iF).Signal_Size_WCM = numel(Data.ixSignal_Final);
        result(iF).WCM_Dn = Data.PosDn;
        result(iF).WCM_RejectionDn = Data.Dn;
        result(iF).Config_Dn = Control.PosDn;
        result(iF).Config_RejectionDn = Control.Dn;
        result(iF).Signal_Components = numel(Data.SignalComp_sizes);
        
        % stats stuff
        load(['Results_batch1/',stat_fnames(iS).name])
        
        result(iF).Animal = statsTable.Animal;
        result(iF).Session = statsTable.Session;
        result(iF).Subvolume = statsTable.Subvolume;
        result(iF).N = statsTable.N;
        result(iF).T = statsTable.T;
        result(iF).Pcorrect = statsTable.Pcorrect;
        result(iF).Pnolick = statsTable.Pnolick;
        result(iF).method = statsTable.method;
        
        % First eigenvalue above 90%
        
        result(iF).eig90 = find(statsTable.explained > 90,1,'first');
    end
    
    
    
end

Network_Rejection_Table = struct2table(result);
save('Results_batch1/Network_Rejection_Table_wStats','Network_Rejection_Table');


%% NRstats FYI
% load('/Users/mathew/work/Peron_crcns/noiserejection/NRstats.mat')
% load('~/work/Peron_crcns/noiserejection/NRstats.mat')
% NRstats fields are:
% 1. Animal, 2. Session, 3. Subvolume, 4. Nretain, 5. Npos, 6. Pretain, 7. normalized Npos
% 8. Original N cells, 9. Duration of recording, 10. P (correct), 11. P (no lick) - measure of task engagement

animals = {'an171923';'an194181';'an194672';'an197522';'an198503';'an229716';'an229717';'an229719'};


%% NO LONGER NEEDED:
%% First, add columns to Network_Rejection_Table corresponding to Animal, session and subvolume
this_a = 0;
this_sess = 0;
for i = 1:height(Network_Rejection_Table)
    
    this_row = Network_Rejection_Table.NetworkName{i};
    an = this_row(1:8);
    a = find(ismember(animals,an)); % Animal
    
    sess = this_row(10:20);
    % Need to iterate over sessions to get a session number
    
    % If this animal and session was seen on the last trial, keep same s.
    % If this animal but a different session was seen, iterate s.
    % If this animal was not seen on the last trial, s = 1;
    if a == this_a
        if sess == this_sess
        else
            s = s+1;
        end
    else
        s = 1;
    end
    
    % Catch for sessions with b in date
    if strcmp(sess(end),'b')
        sv = str2num(this_row(32:end));
    else
        sv = str2num(this_row(31:end));
    end
    
    % Split logic depending on data vs events (N.B. should have specified the
    % same length identifier, but didn't)
    data_YN = this_row(21:24);
    if strcmp(data_YN,'data')
        meth = 'calcium';
    else
        meth = 'Peron';
    end
    
    this_a = a;
    this_sess = sess;
    
    try
        Network_Rejection_Table.Animal(i) = a;
        Network_Rejection_Table.Session(i) = s;
        Network_Rejection_Table.Subvolume(i) = sv;
        Network_Rejection_Table.QC(i) = 1;
        Network_Rejection_Table.method{i} = meth;
        
        display([this_row,', a = ',num2str(a),', s = ',num2str(s),', sv = ',num2str(sv)])
    catch
        Network_Rejection_Table.QC(i) = 0;
    end
end

