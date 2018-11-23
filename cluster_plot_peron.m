%% Code to load and plot peron crcns clustering results

%% Re-generate CCons and save

f = dir('/Volumes/Extras/197522_rejection/Rejected*');
cc = dir('/Volumes/Extras/197522_rejection/CCons_*');
% f(2) = [];
for i = 12: numel(f) %linspace(2,numel(f),numel(f)/2)
    
    % LOAD EVENTS
    enw = f(i).name
    
    % check if this dataset has been clustered already
    oldnews = 0;
    for j = 1:9
        cc_enw = cc(j).name;
        if strcmp(enw(10:end-4),cc_enw(7:end-8))
            oldnews = 1;
        end
    end
    
    if ~oldnews
        % Load rejection data (if not already done)
        load(['/Volumes/Extras/197522_rejection/',enw]);%'/Users/mathew/work/PopWorld/197522_rejection/Rejected_an197522_2013_02_18_events_sv_1.mat')
        
        if Data.Dn > 0
            % Compute consensus clustering matrix
            clusterpars.nreps = 50;
            
            clusterpars.explore = 'explore';  % allow consensus to use more groups than specified by spectral rejection
            
            W = Data.Asignal_final;
            P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model
            B = W - P;
            m = sum(sum(W))/2;
            [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr,~,C] = ...
                ConsensusCommunityDetect(W,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps,[],clusterpars.explore);
            
            for iQ = 1:size(C,2)
                Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering
            end
            
            Allowed = (Q > 0);       % only take consensus using Q-positive clusterings...
            CCons = makeConsensusMatrix(C(:,Allowed));
            
            % Save C and Consensus matrix
            save(['/Volumes/Extras/197522_rejection/CCons_',enw(10:end),'.mat'],'CCons','C')
        end
        
    end
end



%% Look at cluster ID for data vs events versions
% f = dir('example_clustering/*.mat');
f = dir('/Volumes/Extras/197522_rejection/Rejected*'); %dir('Clustering_Results/Clustered_an197522*');
f(1:2) = [];
for i = 1:numel(f)
    try
        if strcmp(f(i).name(30:35),'data_s') % 'events'
            % LOAD EVENTS
            enw = f(i).name;
            %     load(['example_clustering/',f(i-1).name]);
            load(['Clustering_Results/Clustered_',enw(10:end)]);
            
            clf
            subplot(1,2,1)
            yyaxis left
            %             plot(Connected.ConsCluster,'linewidth',2);
            %             plot(conv(Connected.ConsCluster,gausswin(15),'same'),'linewidth',2);
            plot(conv(Connected.QmaxCluster,gausswin(15),'same'),'linewidth',2);
            %             plot(conv(Connected.LouvCluster{1}{1},gausswin(15),'same'),'linewidth',2);
            ylabel('Cluster ID')
            title('Calcium')
            
            % Load data
            load (['/Volumes/05/Peron_2015/Peron_ssc_events/an197522/',enw(10:end-16),'.mat'])
            yyaxis right
            plot(dat.timeSeriesArrayHash.value{1,2}.ids,'linewidth',2);
            ylabel('ROI Id')
            
            %     load(['example_clustering/',f(i).name]);
            % LOAD CALCIUM
            enw_2 = enw;
            enw_2(30:35) = 'events';
            load(['Clustering_Results/Clustered_',enw_2(10:end)]);
            subplot(1,2,2)
            yyaxis left
            %             plot(Connected.ConsCluster,'linewidth',2);
            %             plot(conv(Connected.ConsCluster,gausswin(15),'same'),'linewidth',2);
            plot(conv(Connected.QmaxCluster,gausswin(15),'same'),'linewidth',2);
            %             plot(conv(Connected.LouvCluster{1}{1},gausswin(15),'same'),'linewidth',2);
            
            ylabel('Cluster ID')
            yyaxis right
            plot(dat.timeSeriesArrayHash.value{1,2}.ids,'linewidth',2);
            ylabel('ROI ID')
            
            title('Peron events')
            
            suptitle(['Qmax ',enw(10:end-4)])
            
            drawnow;
            %             pause
            print(['Figures/clustering/image_plane_check/Qmax_',enw(10:end-4)],'-dpdf','-bestfit')
        end
        
    catch
        display(['Something wrong with ',enw])
    end
    
end

%% Image plot version (use C from consensus)


%% Look at CCons while we're here
% f = dir('example_clustering/*.mat');
f = dir('/Volumes/Extras/197522_rejection/CCons*'); %dir('Clustering_Results/Clustered_an197522*');
for i = 1:numel(f)
    try
        if strcmp(f(i).name(27:32),'data_s') % 'events'
            
            % LOAD CCons for calcium data
            enw = f(i).name;
            load(['/Volumes/Extras/197522_rejection/',enw]);
            
            clf
            subplot(1,2,1)
            imagesc(CCons);
            title('Calcium')
            
            %             % Load data
            %             load (['/Volumes/05/Peron_2015/Peron_ssc_events/an197522/',enw(7:end-20),'.mat'])
            %             ids = dat.timeSeriesArrayHash.value{1,2}.ids;
            
            %             ylabel('ROI Id')
            
            %     load(['example_clustering/',f(i).name]);
            % LOAD CALCIUM
            enw_2 = enw;
            enw_2(27:32) = 'events';
            load(['/Volumes/Extras/197522_rejection/',enw_2]);
            subplot(1,2,2)
            imagesc(CCons);
            title('Peron events')
            
            suptitle(['CCons 1 ',enw(7:end-4)])
            
            drawnow;
            %             pause
            print(['Figures/clustering/image_plane_check/CCons_compare_',enw(7:end-8)],'-dpdf','-bestfit')
        end
        
    catch
        %         display(['Something wrong with ',enw])
    end
    
end

%% Distance dependent PCC distributions. All then separate in planes + pairs of planes
% enw = 'CCons_an197522_2013_03_07_data_s_sv_1.mat.mat';
enw = 'CCons_an197522_2013_03_07_events_sv_1.mat.mat';
% Load data
load(['/Volumes/05/Peron_2015/Peron_ssc_events/an197522/',enw(7:end-20),'.mat'])
% Load rejection results for Data.Asignal_final etc
load(['/Volumes/Extras/197522_rejection/Rejected_',enw(7:end-4)]);

% IDs of interest - from noise rejection
ix = Data.ixRetain;

% Get roi IDs for each plane
clear plane_ids
pid_v = [];
roi_order = [];
for i = 1:3
    plane_ids{i} = dat.timeSeriesArrayHash.value{2}.imagingPlane{i}.ids;
    pid_v = [pid_v; i*ones(numel(plane_ids{i}),1)];
    roi_order = [roi_order, 1:numel(plane_ids{i})];
end

retained_ids = dat.timeSeriesArrayHash.value{2}.ids(ix);

% work out which plane each cell is in
cell_z = pid_v(ix);

for i = 1:numel(ix)
    cell_x(i) = mean(dat.timeSeriesArrayHash.descrHash{1,2}.value{1,1}(cell_z(i)).rois(1,roi_order(ix(i))).cornersXY(1,:));
    cell_y(i) = mean(dat.timeSeriesArrayHash.descrHash{1,2}.value{1,1}(cell_z(i)).rois(1,roi_order(ix(i))).cornersXY(2,:));
end

%% Plot cell center locations
% NB. I have forgotten the logic behind 12.75 pixel spacing between planes (8.8.18)
clf

for p = 1:3
    this_p = find(cell_z==p);
    plot3(cell_x(this_p),cell_y(this_p),12.75*cell_z(this_p),'.','markersize',20);
    hold all
end
axis equal
xlabel('X position (pixels)')
ylabel('Y position (pixels)')
zlabel('Z position (scaled)')

print(['Figures/clustering/distance_dependence/Overview_',enw(10:end-4)],'-dpdf','-bestfit')
% set(get(gca,'ZLabel'),'Rotation',0)

%% Compute pairwise distance + angle between each cell
clear D A
A = zeros(numel(ix));
D = zeros(numel(ix));
for i = 1:numel(ix)
    for j = 1:numel(ix)
        if j > i
            a = [cell_x(i), cell_y(i)];
            b = [cell_x(j), cell_y(j)];
            D(i,j) = norm(b-a);
            A(i,j) = atan2(b(2)-a(2),b(1)-a(1));
        end
    end
end


%% Extract relevant PCC/distance numbers into long vectors
idx = find(triu(ones(size(D)),1)); % upper triangular above diagonal;
allDs = D(idx);
allAs = A(idx);
allCs = Data.Asignal_final(idx);

%% Distributions
figure(5);
clf
subplot(1,3,1);
[h,x] = hist(allDs,50);
bar(x,h,'facecolor','k')
xlabel('Distance (pixels)')
ylabel('Count')
axis square

subplot(1,3,2);
[h,x] = hist(allAs,50);
bar(x,h,'facecolor','k')
xlabel('Angle (radians)')
ylabel('Count')
axis square

subplot(1,3,3);
[h,x] = hist(allCs,50);
bar(x,h,'facecolor','k')
xlabel('Pearson correlation coefficient')
ylabel('Count')
axis square

print(['Figures/clustering/distance_dependence/DistDists_',enw(7:end-8)],'-dpdf','-bestfit')

%% Distance dependent cluster ID stats
figure(6)
plot(allDs,allCs,'k.');
xlabel('Distance (pixels)')
ylabel('PCC')

print(['Figures/clustering/distance_dependence/Scatter_',enw(7:end-8)],'-dpdf','-bestfit')
%% Polar plot but with correlations as radial magnitude
rand_ix = randperm(numel(allCs));
sub_ix = rand_ix(1:1e3);

% figure(7)
% polarplot(allAs(sub_ix),allCs(sub_ix),'.')

%% Polar plot but with distance as radial magnitude and colour as normalized PCC
[sorted_pcc, pcc_order] = sort(allCs(sub_ix));
cmap = jet(numel(sub_ix));

% clf
for i = 1:numel(sub_ix)
    polarplot(allAs(sub_ix(pcc_order(i))),allDs(sub_ix(pcc_order(i))),'.','color',cmap(i,:),'markersize',10);
    hold all
end

% print(['Figures/clustering/distance_dependence/Polar_1Ka_',enw(7:end-8)],'-dpdf','-bestfit')

%% Histogram version
% Sort distance into 50 bins of equal amount
[sortD,Dorder] = sort(allDs);
D_bins = linspace(sortD(1),sortD(end),50);
D_edges(1) = 1;%D_bins(1);
D_edges(50) = numel(sortD);%D_bins(50);

for i = 2:49
    [~,mi] = min(abs(sortD - D_bins(i)));
    D_edges(i) = mi;
end

%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%

%% Load typical example data + clustering output


% load('/Users/mathew/work/PopWorld/Clustering_Results/Clustered_an197522_2013_03_07_data_s_sv_1.mat')
% load('/Users/mathew/work/PopWorld/Results/Clustered_Peron_example.mat')
% load('/Users/mathew/work/PopWorld/Results/Rejected_Peron_example.mat')

% load('/Users/mathew/work/PopWorld/197522_rejection/Rejected_an197522_2013_02_18_events_sv_1.mat')
% load('Clustering_Results/Clustered_an197522_2013_02_18_events_sv_1.mat')
% load('/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_02_18.mat'); % 2013_03_07


clear all


FNAMES = {'an171923_2012_06_04_data_s_sv_4';
    'an194181_2013_01_31_data_s_sv_2';
    'an194672_2013_02_24_data_s_sv_1';
    'an197522_2013_03_08_data_s_sv_2';
    'an198503_2013_03_15_data_s_sv_5';
    'an229716_2013_12_06_data_s_sv_6';
    'an229717_2013_12_01_data_s_sv_3';
    'an229719_2013_12_01_data_s_sv_6'};


% load('/Users/mathew/work/PopWorld/Results_reject_preround/Network_Rejection_Table_wStats_preround3.mat')
load('Results_reject_preround/Network_Rejection_Table_wEvents.mat');

% Results_path = '/Volumes/Extras/197522_rejection/Rejection_tests/';
% fname = 'Pre_conv_1K_20_02'; % 'Pre_conv_1K_21_02';

% load([Results_path,'Scaling_Results_Table'])
% load([Results_path,'/Rejected_',fname],'Data','Rejection')
% load([Results_path,'/Clustered_',fname],'Full','Connected')


% fname = 'an171923_2012_06_04_data_s_sv_4';
% fname = 'an194181_2013_01_31_data_s_sv_2';
% fname = 'an194672_2013_02_24_data_s_sv_1';
% fname = 'an197522_2013_03_08_data_s_sv_2';
% fname = 'an198503_2013_03_15_data_s_sv_5';
% fname = 'an229716_2013_12_06_data_s_sv_6';
% fname = 'an229717_2013_12_01_data_s_sv_3';
% fname = 'an229719_2013_12_01_data_s_sv_6';
for i = 41:47; figure(i); clf; end
% POS_mice = [1,2,4,6,8]; NEG_mice = [3,5,7]; % Ca
POS_mice = [4,5,6,8]; NEG_mice = [1,2,3,7]; % Ca

for A = 1:8
% A = 4;
    fname = FNAMES{A};
    
    % Clustering results if available
    % load(['/Users/mathew/work/PopWorld/Clustering_Results_preround/Clustered_',fname,'.mat']); % _data_s_sv_1
    
    %     load(['/Volumes/05/Peron_2015/Peron_ssc_events/an197522/',fname(1:19),'.mat'])
    %     load(['/Volumes/05/Peron_2015/Peron_ssc-2/ssc-2/',fname(1:8),'/',fname(1:19),'_data_struct.mat'])
    if A~=4
        load(['/media/mathew/Data_1/Peron_ssc-2/with_events/',fname(1:8),'/',fname(1:19),'_sess.mat']);
        dat = s; clear s;
    else
        load(['/media/mathew/Data_1/Peron_ssc-2/with_events/',fname(1:8),'/',fname(1:19),'.mat']);
    end
    %     dat = s; clear s;
    sv = str2num(fname(end));
    
    N = numel(dat.timeSeriesArrayHash.value);
    ca_files = 2 : ((N+1)/2);
    ev_files = ((N+1)/2)+1 : N;
    %% Re-run noise rejection
    
%     if ~(exist(['Results_reject_preround/Rejected_',fname,'.mat']))
    if ~(exist(['Results_reject_preround/Rejected_',fname,'_events.mat']))
        data = dat.timeSeriesArrayHash.value{ev_files(sv)}.valueMatrix;
        data(find(isnan(data))) = 0;
        [~,~,Data,Rejection,~] = local_NR(data);
        
%         save(['Results_reject_preround/Rejected_',fname,'.mat'],'Rejection','Data','-v7.3')
        save(['Results_reject_preround/Rejected_',fname,'_events.mat'],'Rejection','Data','-v7.3')

    else
        data = dat.timeSeriesArrayHash.value{ev_files(sv)}.valueMatrix;
        data(find(isnan(data))) = 0;
%         load(['Results_reject_preround/Rejected_',fname,'.mat'])
        load(['Results_reject_preround/Rejected_',fname,'_events.mat'])
    end
    
    %% Compute eig decomposition of B
    B = Data.A - Data.ExpA;
    [V,egs] = eig(B);  % eigenspectra of data modularity matrix
    egs = diag(egs); % extract vector from diagonal
    [egs,ix] = sort(egs,'descend'); % sort eigenvalues into descending order
    V = V(:,ix);  % sort eigenvectors accordingly
    
    %% For first 5 egs, find high and low projecting cells
    clear negx posx
    for i = 1:5
        %    negx{i} = find(V(:,i)<= mean(V(:,i))-2*std(V(:,i)));
        %    posx{i} = find(V(:,i)>= mean(V(:,i))+2*std(V(:,i)));
        
        negx{i} = find(V(:,i)<= prctile(V(:,i),5));
        posx{i} = find(V(:,i)>= prctile(V(:,i),95));
    end
    
    %% Trial bounaries and touches
    T = dat.timeSeriesArrayHash.value{ev_files(sv)}.time;
    x = dat.timeSeriesArrayHash.value{1,ev_files(sv)}.trial;
    trialId = find(diff(x));
    Touch_G = dat.eventSeriesArrayHash.value{2}.eventTimes{1};
    Touch_NG = dat.eventSeriesArrayHash.value{2}.eventTimes{2};
    
    
    %% Reconstruct activity based on V
    % data(find(isnan(data))) = 0;
    cmap = lines(5);
    cmap2 = redblue(5);
%     for i = 1:5
%         NEG = -V(negx{i},i)' * data(Data.ixRetain(negx{i}),:);
%         % NEG(find(isnan(NEG))) = 0;
%         
%         POS = V(posx{i},i)' * data(Data.ixRetain(posx{i}),:);
%         % POS(find(isnan(POS))) = 0;
%         
%         figure(i); clf;
%         plot(T,zscore(POS),'color',cmap(1,:));
%         hold all
%         plot(T,conv(zscore(POS),ones(200,1),'same')/200,'color',cmap2(2,:),'linewidth',2)
%         plot(T,zscore(NEG),'color',cmap(2,:));
%         plot(T,conv(zscore(NEG),ones(200,1),'same')/200,'color',cmap2(4,:),'linewidth',2)
%         title(['Eig ',num2str(i),'. ',num2str(numel(posx{i})),' pos, ',num2str(numel(negx{i})),' neg cells'])
%         
%         % Trial stuff
%         plot([T(trialId);T(trialId)],[-2.5*ones(1,numel(trialId));-3.5*ones(1,numel(trialId))],'k')
%         
%         plot(Touch_G,-2*ones(1,numel(Touch_G)),'r*')
%         plot(Touch_NG,-2.1*ones(1,numel(Touch_NG)),'b*')
%         legend('positive projecting','Smoothed Pos','negative projecting','Smoothed Neg')
%         xlim([T(1),T(end)])
%         %     pause
% %         print(['Figures/V1_LowDprojection/LowD_prctile_eig_',num2str(i),'_reconstruction_',fname],'-dpdf','-bestfit')
%         
%     end
%     
    %% PSTHs of NEG and POS groups
    x = dat.timeSeriesArrayHash.value{1,ev_files(sv)}.trial;
    trials = unique(x);
    
    % Trial type - Left or right
    L = unique([find(dat.trialTypeMat(1,trials)), find(dat.trialTypeMat(3,trials))]);
    R = unique([find(dat.trialTypeMat(2,trials)), find(dat.trialTypeMat(4,trials))]);
    
    % Outcome - Correct, incorrect
    C = unique([find(dat.trialTypeMat(1,trials)), find(dat.trialTypeMat(2,trials))]);
    IC = unique([find(dat.trialTypeMat(3,trials)), find(dat.trialTypeMat(4,trials))]);
    
    % Hit, miss, correct rejection, false alarm
    hit = L(ismember(L,C));
    CR = R(ismember(R,C));
    miss = L(ismember(L,IC));
    FA = R(ismember(R,IC));
    
    clear psth_ca psth_ca_L  psth_ca_R  psth_ca_C  psth_ca_IC
    clear psth_ca_hit psth_ca_CR psth_ca_miss psth_ca_FA
    
    for p = 1:5
        neg_ca = [];
        pos_ca = [];
        
        NEG = -V(negx{p},p)' * data(Data.ixRetain(negx{p}),:);
        
        POS = V(posx{p},p)' * data(Data.ixRetain(posx{p}),:);
        
        for i = 1:numel(trials)
            neg_ca(i,:) = [NEG(find(x == trials(i))),nan(1,100-numel(find(x == trials(i))))];
            pos_ca(i,:) = [POS(find(x == trials(i))),nan(1,100-numel(find(x == trials(i))))];
        end
        
        psth_ca(p,:) = nanmean(pos_ca(:,1:50),1); psth_ca(p+6,:) = nanmean(neg_ca(:,1:50),1);
        
        psth_ca_L(p,:) = nanmean(pos_ca(L,1:50),1); psth_ca_L(p+6,:) = nanmean(neg_ca(L,1:50),1);
        psth_ca_R(p,:) = nanmean(pos_ca(R,1:50),1); psth_ca_R(p+6,:) = nanmean(neg_ca(R,1:50),1);
        psth_ca_C(p,:) = nanmean(pos_ca(C,1:50),1); psth_ca_C(p+6,:) = nanmean(neg_ca(C,1:50),1);
        psth_ca_IC(p,:) = nanmean(pos_ca(IC,1:50),1); psth_ca_IC(p+6,:) = nanmean(neg_ca(IC,1:50),1);
        
        psth_ca_hit(p,:) = nanmean(pos_ca(hit,1:50),1); psth_ca_hit(p+6,:) = nanmean(neg_ca(hit,1:50),1);
        psth_ca_CR(p,:) = nanmean(pos_ca(CR,1:50),1); psth_ca_CR(p+6,:) = nanmean(neg_ca(CR,1:50),1);
        psth_ca_miss(p,:) = nanmean(pos_ca(miss,1:50),1); psth_ca_miss(p+6,:) = nanmean(neg_ca(miss,1:50),1);
        psth_ca_FA(p,:) = nanmean(pos_ca(FA,1:50),1); psth_ca_FA(p+6,:) = nanmean(neg_ca(FA,1:50),1);
        
    end
    
    
    
    %%%%%%%%%%%%%%%%% NEWER CoSyNe abstract code %%%%%%
    %% Calcium times
    x = dat.timeSeriesArrayHash.value{1,ev_files(sv)}.trial;
    trials = unique(x);
    
    trialStarts = [];
    for i = 1:numel(trials)
        trialStarts(i) = dat.trialStartTimes(find(dat.trialIds == trials(i)));
        ca_t(i,:) = [dat.timeSeriesArrayHash.value{ev_files(sv)}.time(x==trials(i))-trialStarts(i),nan(1,100-numel(find(x == trials(i))))];
    end
    
    ca_ts = nanmean(ca_t(:,1:50));
    
    
    %% Event array
    events = []; touches = []; trialStarts = []; dat.trialStartTimes(trials);
    for i = 1:numel(trials)
        trialStarts(i) = dat.trialStartTimes(find(dat.trialIds == trials(i)));
        % Pole up
        x = find(dat.eventSeriesArrayHash.value{1}.eventTrials == trials(i),1,'first');
        events(1,i) = max([0,dat.eventSeriesArrayHash.value{1}.eventTimes(x) - trialStarts(i)]);
        % Pole down
        x = find(dat.eventSeriesArrayHash.value{1}.eventTrials == trials(i),1,'last');
        events(2,i) = max([0,dat.eventSeriesArrayHash.value{1}.eventTimes(x) - trialStarts(i)]);
        % Lick left
        x = find(dat.eventSeriesArrayHash.value{3}.eventTrials == trials(i),1,'first');
        events(3,i) = max([0,dat.eventSeriesArrayHash.value{3}.eventTimes(x) - trialStarts(i)]);
        % Lick right
        x = find(dat.eventSeriesArrayHash.value{4}.eventTrials == trials(i),1,'first');
        events(4,i) = max([0,dat.eventSeriesArrayHash.value{4}.eventTimes(x) - trialStarts(i)]);
        % Reward left
        x = find(dat.eventSeriesArrayHash.value{5}.eventTrials == trials(i),1,'first');
        events(5,i) = max([0,dat.eventSeriesArrayHash.value{5}.eventTimes(x) - trialStarts(i)]);
        % Reward right
        x = find(dat.eventSeriesArrayHash.value{6}.eventTrials == trials(i),1,'first');
        events(6,i) = max([0,dat.eventSeriesArrayHash.value{6}.eventTimes(x) - trialStarts(i)]);
        % reward cue
        x = find(dat.eventSeriesArrayHash.value{7}.eventTrials == trials(i),1,'first');
        events(7,i) = max([0,dat.eventSeriesArrayHash.value{7}.eventTimes(x) - trialStarts(i)]);
        x = find(dat.eventSeriesArrayHash.value{2}.eventTrials{1} == trials(i),1,'first');
        events(8,i) = max([0,dat.eventSeriesArrayHash.value{2}.eventTimes{1}(x) - trialStarts(i)]);
        x = find(dat.eventSeriesArrayHash.value{2}.eventTrials{2} == trials(i),1,'first');
        events(9,i) = max([0,dat.eventSeriesArrayHash.value{2}.eventTimes{2}(x) - trialStarts(i)]);
        
    end
    % Medians
    for i = 1:size(events,1)
        med_events(i) = median(events(i,find(events(i,:))));
    end
    
    %% Plot psths of all groups
    
    % loop to plot event lines on top (pole up,down,reward port)
    
%     for f = [20,21,22]
%         figure(f);clf;
%         if f ~= 20
%             for j = 1:4
%                 subplot(2,2,j)
%                 event_subset = [1,2,7];
%                 for i = 1:numel(event_subset)
%                     plot(med_events(event_subset(i))*ones(2,1),[0,1.5],'linewidth',2,'color','k'); %[0,0.1]
%                     hold all
%                 end
%             end
%         else
%             event_subset = [1,2,7];
%             for i = 1:numel(event_subset)
%                 plot(med_events(event_subset(i))*ones(2,1),[0,1.5],'linewidth',2,'color','k');
%                 hold all
%             end
%         end
%         
%     end
    
    %%
%     cmap3 = brewermap(10,'Paired');
%   
%     for i = 1:5
%         figure(20);
%         plot(ca_ts,psth_ca(i,:)','color',cmap3((2*i)-1,:),'linewidth',2)
%         plot(ca_ts,psth_ca(i+6,:)','color',cmap3(2*i,:),'linewidth',2)
%         title(['Population PSTH (',num2str(numel(trials)),' trials)'])
%         xlim([0,7000])
%         legend('Pole up','Pole down','Reward cue')
%         xlabel('Time (ms)')
%         ylabel('Cluster mean Ca2+')
%         % ylim([0,0.1])
%         % pause
%     end
%     print(['Figures/V1_LowDprojection/PSTH_all_prctile_',fname,'_Events'],'-dpdf','-bestfit')
    %%
    % Two 2x2 plots of the trial type/ performance data
%     for i = 1:5
%         figure(21);
%         ax(1) = subplot(2,2,1);
%         plot(ca_ts,psth_ca_L(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Left (',num2str(numel(L)),' trials)'])
%         plot(ca_ts,psth_ca_L(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         ax(2) = subplot(2,2,2);
%         plot(ca_ts,psth_ca_R(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Right (',num2str(numel(R)),' trials)'])
%         plot(ca_ts,psth_ca_R(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         ax(3) = subplot(2,2,3);
%         plot(ca_ts,psth_ca_C(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Correct (',num2str(numel(C)),' trials)'])
%         plot(ca_ts,psth_ca_C(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         xlabel('Time (ms)')
%         ylabel('Cluster mean Ca2+')
%         ax(4) = subplot(2,2,4);
%         plot(ca_ts,psth_ca_IC(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Incorrect (',num2str(numel(IC)),' trials)'])
%         plot(ca_ts,psth_ca_IC(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         linkaxes(ax)
%         xlim([0,7000])
%         % ylim([0,0.1])
%         
%         figure(22);
%         bx(1) = subplot(2,2,1);
%         plot(ca_ts,psth_ca_hit(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Correct Left (',num2str(numel(hit)),' trials)'])
%         plot(ca_ts,psth_ca_hit(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         bx(2) = subplot(2,2,2);
%         plot(ca_ts,psth_ca_CR(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Correct Right (',num2str(numel(CR)),' trials)'])
%         plot(ca_ts,psth_ca_CR(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         bx(3) = subplot(2,2,3);
%         plot(ca_ts,psth_ca_miss(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Incorrect Left (',num2str(numel(miss)),' trials)'])
%         plot(ca_ts,psth_ca_miss(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         xlabel('Time (ms)')
%         ylabel('Cluster mean Ca2+')
%         bx(4) = subplot(2,2,4);
%         plot(ca_ts,psth_ca_FA(i,:)','color',cmap3((2*i)-1,:),'linewidth',2); title(['Incorrect Right (',num2str(numel(FA)),' trials)'])
%         plot(ca_ts,psth_ca_FA(i+6,:)','color',cmap3(2*i,:),'linewidth',2);
%         linkaxes(bx)
%         xlim([0,7000])
%         % ylim([0,0.1])
%         % pause
%     end
%     
%     figure(21);
%     print(['Figures/V1_LowDprojection/PSTH_TrialType_prctile_',fname,'_Events'],'-dpdf','-bestfit')
%     
%     figure(22);
%     print(['Figures/V1_LowDprojection/PSTH_Outcome_prctile_',fname,'_Events'],'-dpdf','-bestfit')
    
    %% One eigenvector per animal

        figure(41);
        cmap3 = brewermap(10,'Paired');
        ax(A) = subplot(4,2,A)
        event_subset = [1,2,7];
        for i = 1:numel(event_subset)
            plot(med_events(event_subset(i))*ones(2,1),[0,1.5],'linewidth',2,'color','k'); %[0,0.1]
            hold all
        end
        
        plot(ca_ts,psth_ca(1,:)','color',cmap3(1,:),'linewidth',2)
        plot(ca_ts,psth_ca(7,:)','color',cmap3(2,:),'linewidth',2)
        
        title([fname(1:8),' (',num2str(numel(trials)),' trials)'])
        xlim([0,7000])
%         legend('Pole up','Pole down','Reward cue','Positive','Negative');
        xlabel('Time (ms)')
%         ylabel('Cluster mean Ca2+')
        ylabel('Ca2+')
        % ylim([0,0.1])
        % pause
    
    %% One eigenvector per animal, 'positive' pop only, overlaid
    % Also store population pole up times, ca_ts and psth_ca for plotting later
        figure(45);
        cmap4 = brewermap(8,'Dark2');

        event_subset = [1,2,7];
        for i = 1:numel(event_subset)
            plot(med_events(event_subset(i))*ones(2,1),[0,1.5],'linewidth',2,'color',cmap4(A,:)); %[0,0.1]
            pop_med_events(A,i) = med_events(event_subset(i));
            hold all
        end
        
        if any(ismember(POS_mice,A))
            plot(ca_ts,psth_ca(1,:)' - psth_ca(1,1),'color',cmap4(A,:),'linewidth',2)
            pop_ca_ts(A,:) = ca_ts;
            pop_psth_ca(A,:) = psth_ca(1,:)' - psth_ca(1,1);
        elseif any(ismember(NEG_mice,A))
            plot(ca_ts,psth_ca(7,:)' - psth_ca(7,1),'color',cmap4(A,:),'linewidth',2)
            pop_ca_ts(A,:) = ca_ts;
            pop_psth_ca(A,:) = psth_ca(7,:)' - psth_ca(7,1);
        end
        
        xlim([0,7000])
        xlabel('Time (ms)')
%         ylabel('Cluster mean Ca2+')
        ylabel('Ca2+')
        % ylim([0,0.1])
        % pause
    
        
    %% Continuous singals - touch, curvature and angle
    %% Load whisker angle and curvature
    Ang = dat.timeSeriesArrayHash.value{1}.valueMatrix(1,:);
    Ang(isnan(Ang)) = 0;
    K = dat.timeSeriesArrayHash.value{1}.valueMatrix(2,:);
    K(isnan(K)) = 0;
    
    % whisker variable times
    W_t = dat.timeSeriesArrayHash.value{1}.time;
    
    %% Compute downsampled version of A and K such that the value at each
    % time point is the mean of the higher sampling rate data
    % How about we just convolve the original time series with a 142ms boxcar,
    % emulating averaging?
    
    % Hilbert amplitude first
    t = Ang;
    t = t';
    % fix nans
    t(find(isnan(t))) = 0;
    % whisk amplitude
    clear ts
    ts = timeseries(t,(1:length(t))./1000);
    bandpass = [6 30];
    theta_filt = idealfilter(ts,bandpass,'pass');
    % Hilbert transform to find amplitude and phase
    H = hilbert(theta_filt.data);
    h = squeeze(abs(H));
    
    
    A_conv = conv(h,ones(142,1),'same')/142;
    K_conv = conv(K,ones(142,1),'same')/142;
    clear A_ds K_ds
    for i = 1:numel(T)
        [~,this_t] = min(abs(W_t-T(i)));
        A_ds(i) = A_conv(this_t);
        K_ds(i) = K_conv(this_t);
    end
    
    
    %% Create tuning curves of for Ang and K
    for p = 1%:5
        NEG = -V(negx{p},p)' * data(Data.ixRetain(negx{p}),:);
        POS = V(posx{p},p)' * data(Data.ixRetain(posx{p}),:);
%         figure(7+p);clf
        % Run either NEG or POS depending on mouse
        % POS_mice = [1,2,4,6,8]; NEG_mice = [3,5,7];
%         for j = 1:2
%             if j == 1
%                 ca = POS;
% %                 figure(42);
%             elseif j == 2
%                 ca = NEG;
% %                 figure(43);
%             end
            if any(ismember(POS_mice,A))
                ca = POS;
            elseif any(ismember(NEG_mice,A))
                ca = NEG;
            end
            figure(46);
            
            numbins = 15;
            bin_edges = round(linspace(1,numel(K_ds),numbins+1));
            [~,sort_K] = sort((K_ds));
            clear sum_ca_K  mean_ca_K std_ca_K mean_kappa
            for i = 1:numel(bin_edges)-1
                sum_ca_K(i) = sum(ca(sort_K(bin_edges(i):bin_edges(i+1))));
                mean_ca_K(i) = nanmean(ca(sort_K(bin_edges(i):bin_edges(i+1))));
                std_ca_K(i) = nanstd(ca(sort_K(bin_edges(i):bin_edges(i+1))));
                mean_kappa(i) = nanmean((K_ds(sort_K(bin_edges(i):bin_edges(i+1)))));
            end
            
%                         ax(2) = subplot(2,2,(j*2)-1);
%                         subplot(8,2,(A*2)-1); cla;
            subplot(2,1,1); hold all;
%                         myeb(mean_kappa,mean_ca_K,std_ca_K/sqrt(bin_edges(2)-1),[0,0,1,0])
%                         myeb(mean_kappa,mean_ca_K,std_ca_K/sqrt(bin_edges(2)-1),cmap4(A,:),cmap4(A,:))
            plot(mean_kappa,mean_ca_K- mean_ca_K(8),'linewidth',2,'color',cmap4(A,:))
            title([fname(1:8)])
            ylabel('Mean Normalized Ca (bin)');
            xlabel('Whisker curvature');
            
            [~,sort_A] = sort((A_ds));
            clear sum_ca_A  mean_ca_A std_ca_A mean_amp
            for i = 1:numel(bin_edges)-1
                sum_ca_A(i) = sum(ca(sort_A(bin_edges(i):bin_edges(i+1))));
                mean_ca_A(i) = nanmean(ca(sort_A(bin_edges(i):bin_edges(i+1))));
                std_ca_A(i) = nanstd(ca(sort_A(bin_edges(i):bin_edges(i+1))));
                mean_amp(i) = nanmean((A_ds(sort_A(bin_edges(i):bin_edges(i+1)))));
            end
            
            pop_mean_ca_A(A,:) = mean_ca_A;
            pop_mean_amp(A,:) = mean_amp;
            N_amp(A) = numel(K_ds);
%                         ax(1) = subplot(2,2,j*2);
%                         subplot(8,2,(A*2)); cla
            subplot(2,1,2); hold all
%                         myeb(mean_amp,mean_ca_A,std_ca_A/sqrt(bin_edges(2)-1),[1,0.1,0,0])
%                         myeb(mean_amp,mean_ca_A,std_ca_A/sqrt(bin_edges(2)-1),cmap4(A,:),cmap4(A,:))
            plot(mean_amp,mean_ca_A-mean_ca_A(1),'linewidth',2,'color',cmap4(A,:))
            ylabel('Mean Normalized Ca (bin)');
            xlabel('Whisker amplitude');
        end
%         suptitle(['V',num2str(p),'. ',fname(1:8)])
%         print(['Figures/V1_LowDprojection/Kappa_tuning_V',num2str(p),'_',fname,'_Events'],'-dpdf','-bestfit')
%     end
    
    %% Touch triggered average
    % First find touches on relevant trials
    x = dat.timeSeriesArrayHash.value{1,ev_files(sv)}.trial;
    trials = unique(x);
    
    % Find the first touch on every trial
    first_touch = [];
    for i = 1:numel(trials)
        tr = trials(i);
        tt = [];
        x1 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{1} ==tr,1,'first');
        tt1 = dat.eventSeriesArrayHash.value{2}.eventTimes{1}(x1);
        x2 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{2} ==tr,1,'first');
        tt2 = dat.eventSeriesArrayHash.value{2}.eventTimes{2}(x2);
        tt = [tt1,tt2];
        if find(tt)
            first_touch = [first_touch, min(tt(find(tt)))];
        end
        
    end
    
    % Closest frame to touch time
    ftf = [];
    for i = 1:numel(first_touch)
        [mn,im] = min(abs(dat.timeSeriesArrayHash.value{ev_files(sv)}.time - first_touch(i)));
        ftf(i) = im; % dat.timeSeriesArrayHash.value{2}.time(im);
    end
    
%     figure(31); clf;
    figure(47);
    for p = 1%:5
        NEG = -V(negx{p},p)' * data(Data.ixRetain(negx{p}),:);
        POS = V(posx{p},p)' * data(Data.ixRetain(posx{p}),:);
        
        clear c_POS c_NEG
        for i = 1:numel(ftf)
            c_POS(i,:) = POS(ftf(i)-7 : ftf(i)+14); %12); % Formerly -2, +12
            c_NEG(i,:) = NEG(ftf(i)-7 : ftf(i)+14); %12);
        end
        
        POS_tta = nanmean(c_POS,1);
        NEG_tta = nanmean(c_NEG,1);
        
%         subplot(3,2,p);
%         subplot(4,2,A); cla
%         myeb(1:22,mean(zscore(c_POS')'),std(zscore(c_POS')')/sqrt(length(c_POS)),[0,1,0,0])
%         myeb(1:22,mean(zscore(c_NEG')'),std(zscore(c_NEG')')/sqrt(length(c_NEG)),[0,0.1,1,0])

        if any(ismember(POS_mice,A))
%             myeb(1:22,mean(zscore(c_POS')'),std(zscore(c_POS')')/sqrt(length(c_POS)),cmap4(A,:),cmap4(A,:))
            tth = mean(zscore(c_POS')');
            plot(1:22,tth-tth(1),'linewidth',2,'color',cmap4(A,:))
            pop_mean_tth(A,:) = tth-tth(1);
            N_tth(A) = numel(ftf);
        elseif any(ismember(NEG_mice,A))
%             myeb(1:22,mean(zscore(c_NEG')'),std(zscore(c_NEG')')/sqrt(length(c_NEG)),cmap4(A,:),cmap4(A,:))
            tth = mean(zscore(c_NEG')');
            plot(1:22,tth-tth(1),'linewidth',2,'color',cmap4(A,:))
            pop_mean_tth(A,:) = tth-tth(1);
            N_tth(A) = numel(ftf);
        end
        hold all
%         plot([8,8],[-1,1],'k--')
        
        title(['Touch triggered average. V',num2str(p)])
%         title([fname(1:8)])
        %     pause
    end
%     print(['Figures/V1_LowDprojection/Touch_tuning_',fname,'_Events'],'-dpdf','-bestfit')
end

%%
figure(48); clf;
leg_str = []; for i = 1:8; plot([1,1],[1,2],'linewidth',2,'color',cmap4(i,:)); hold all; leg_str = [leg_str;'Mouse ',num2str(i)];end
legend(leg_str)

%% All mice Amplitude and touch coding plots 
figure(48); clf;
subplot(1,2,1);
plot(pop_mean_amp',pop_mean_ca_A','linewidth',2,'color',[.7,.7,.7]); hold all
plot(mean(pop_mean_amp),mean(pop_mean_ca_A),'linewidth',2,'color',[1,0.1,0])
ylabel('Population Ca2+')
xlabel('Whisker amplitude (degrees)')
title('Whisking amplitude')
axis square
%            N_amp
subplot(1,2,2);
plot(pop_mean_tth','linewidth',2,'color',[.7,.7,.7]); hold all
plot(mean(pop_mean_tth),'linewidth',2,'color',[0.1,0,1])
set(gca,'Xtick',linspace(3,23,5),'XTicklabel',linspace(-5,15,5));
ylabel('Population Ca2+')
xlabel('Time to touch (frames)')
title('Touch')
plot([8,8],[-0.5,3],'k--','linewidth',2)
axis square

%% Pop psth 
% pop_med_events
% pop_ca_ts
% pop_psth_ca
figure(49); clf;
mpme = mean(pop_med_events);
fill([mpme(1);mpme(2)*ones(2,1);mpme(1)],[-0.5,-0.5,3,3],[1,0.8,0.8]); hold all
fill([mpme(2);mpme(3)*ones(2,1);mpme(2)],[-0.5,-0.5,3,3],[0.8,0.8,1]); hold all
for i = 1:3
    plot(mpme(i)*ones(2,1),[-0.5,3],'k--','linewidth',2); 
end

plot(pop_ca_ts'-pop_med_events(:,2)'+mean(pop_med_events(:,2)),pop_psth_ca','linewidth',2,'color','k')
ylabel('Population Ca2+')
xlabel('Time (ms)')
xlim([0,7000])
% plot(mean([pop_ca_ts'-pop_med_events(:,2)']')+mean(pop_med_events(:,2)),mean(pop_psth_ca),'linewidth',2,'color',[.7,0,.7])
%% OLDER CLUSTER PLOT CODE %%%%%%%%%%%%%%%%%%%%%%
%% High skewness neurons
ca = dat.timeSeriesArrayHash.value{2}.valueMatrix;
ev = dat.timeSeriesArrayHash.value{3}.valueMatrix;
S = skewness(ca');
S(isnan(S)) = [];
%% Image data by order of skewness
[~,skew_order] = sort(S);
imagesc(ca(skew_order,:))

%% Plot 10 most and least skewed cells
figure
plot(ca(skew_order(numel(S)-10:numel(S)),:)','r')
hold all
plot(ca(skew_order(1:10),:)','k')
%% Visulise clusters of most skewy neurons
% skewy = find(S>2);
skewy = 1:numel(S);
% temp_C = Full.QmaxCluster(skewy);
temp_C = Full.ConsCluster(skewy);
% temp_C = Full.LouvCluster{1}{1}(skewy);
unique_C = unique(temp_C);
skewy_C = zeros(numel(skewy),1);
for i = 1:numel(unique_C)
    this_C = find(temp_C ==unique_C(i));
    skewy_C(this_C) = i;
    
end

% Identify groups in similarity order
figure
[H,h,Ix] = plotClusterMap(Data.A(skewy,skewy),skewy_C,[],[],'S');

% Ordered clusters
n = 1;
newG = zeros(size(skewy_C));
for i = 1:numel(unique(skewy_C))
    this_c = skewy_C(Ix(n));
    these_c = find(skewy_C(Ix) == this_c);
    newG(these_c) = i;
    n = these_c(end) + 1;
end

%% Image first 50 or so
% figure
imagesc(ca(skewy(Ix(1:500)),:)); % 1:50
% imagesc(ML_s(skewy(Ix(1:500)),:)); % 1:50
%% Zooming in
figure(19);
imagesc(ca(skewy(Ix(20:29)),:)); % 1:50 %20:29 %13:21
figure(20);
plot(ca(skewy(Ix(20:29)),:)')

%% Plot each group in turn
% figure(21); clf
T = dat.timeSeriesArrayHash.value{2}.time;
for i = 1: numel(unique(newG))
    figure
    this_c = find(newG == i);
    %     plot(T,ca(skewy(Ix(this_c)),:)');
    s_ca = ca(skewy(Ix(this_c)),:);
    s_ca(isnan(s_ca)) = 0;
    %     imagesc(T,1:numel(this_c),x')
    plot(T,100*nanmean(s_ca))
    title(['Cluster ',num2str(i)])
    %     pause
    %% Get trial boundaries and plot on top
    hold all
    x = dat.timeSeriesArrayHash.value{1,2}.trial;
    trialId = find(diff(x));
    %     plot([T(trialId);T(trialId)],[-0.5*ones(1,numel(trialId));-1.5*ones(1,numel(trialId))],'k')
    plot([T(trialId);T(trialId)],[-2.5*ones(1,numel(trialId));-15*ones(1,numel(trialId))],'k')
    
    %% Plot touches
    Touch_G = dat.eventSeriesArrayHash.value{2}.eventTimes{1};
    Touch_NG = dat.eventSeriesArrayHash.value{2}.eventTimes{2};
    
    plot(Touch_G,-2*ones(1,numel(Touch_G)),'r*')
    plot(Touch_NG,-2.1*ones(1,numel(Touch_NG)),'b*')
    
    %     ylim([-15,numel(this_c)])
    ylim([-15,100*max(nanmean(s_ca)) + 1])
    %     pause
    
end

%% Load whisker angle and curvature
A = dat.timeSeriesArrayHash.value{1}.valueMatrix(1,:);

K = dat.timeSeriesArrayHash.value{1}.valueMatrix(2,:);

% whisker variable times
W_t = dat.timeSeriesArrayHash.value{1}.time;

%% Dense touch time array
T_t = zeros(2,length(W_t));
for i = 1:numel(Touch_G)
    this_t = Touch_G(i);
    [~,ix_t] = min(abs(W_t-this_t)) ;
    T_t(1,ix_t) = 1;
end
for i = 1:numel(Touch_NG)
    this_t = Touch_NG(i);
    [~,ix_t] = min(abs(W_t-this_t)) ;
    T_t(2,ix_t) = 1;
end




%% Load sparse event version of data + generate event array
% load('/Volumes/05/Peron_2015/Deconvolution_test/ML_peron.mat')
% ncells = numel(ML_peron);
% nt = numel(ML_peron{1}.fit);
% ML_s = zeros(ncells,nt);
% for t = 1:ncells
%     try
%         spikes = 7*ML_peron{t}.spikest1; % ML spike output is in dt space
%
%         % Generate dense array of spikes
%         spk_dense = zeros(nt,1);
%         spk_dense(round(spikes)) = 1;
%         ML_s(t,:) = spk_dense;
%         %     events = conv(spk_dense,kernel)/sum(kernel);
%         %     ML_e(t,:) = events(1:nt);
%     end
% end

%% PSTHs of groups
% First get trial type info

%% PSTHs
x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

psth_ca = []; psth_ca_L = []; psth_ca_R = []; psth_ca_C = []; psth_ca_IC = [];

% Trial type - Left or right
L = unique([find(dat.trialTypeMat(1,trials)), find(dat.trialTypeMat(3,trials))]);
R = unique([find(dat.trialTypeMat(2,trials)), find(dat.trialTypeMat(4,trials))]);

% Outcome - Correct, incorrect
C = unique([find(dat.trialTypeMat(1,trials)), find(dat.trialTypeMat(2,trials))]);
IC = unique([find(dat.trialTypeMat(3,trials)), find(dat.trialTypeMat(4,trials))]);

% Hit, miss, correct rejection, false alarm
hit = L(ismember(L,C));
CR = R(ismember(R,C));
miss = L(ismember(L,IC));
FA = R(ismember(R,IC));

clear psth_ca psth_ca_L  psth_ca_R  psth_ca_C  psth_ca_IC
clear psth_ca_hit psth_ca_CR psth_ca_miss psth_ca_FA

for p = 1:numel(unique(newG))
    c_ca = [];
    s_ca = [];
    
    this_c = find(newG == p);
    
    s_ca = ca(this_c,:);
    s_ca(isnan(s_ca)) = 0;
    
    pop_mean = nanmean(s_ca);
    
    for i = 1:numel(trials)
        c_ca(i,:) = [pop_mean(find(x == trials(i))),nan(1,100-numel(find(x == trials(i))))];
        
    end
    
    psth_ca(p,:) = nanmean(c_ca(:,1:50),1);
    
    
    
    psth_ca_L(p,:) = nanmean(c_ca(L,1:50),1);
    psth_ca_R(p,:) = nanmean(c_ca(R,1:50),1);
    
    psth_ca_C(p,:) = nanmean(c_ca(C,1:50),1);
    psth_ca_IC(p,:) = nanmean(c_ca(IC,1:50),1);
    
    
    psth_ca_hit(p,:) = nanmean(c_ca(hit,1:50),1);
    psth_ca_CR(p,:) = nanmean(c_ca(CR,1:50),1);
    
    psth_ca_miss(p,:) = nanmean(c_ca(miss,1:50),1);
    psth_ca_FA(p,:) = nanmean(c_ca(FA,1:50),1);
    
end



%% Create long dense vectors of all event data

% whisker variable times
W_t = dat.timeSeriesArrayHash.value{1}.time;
Touch_G = dat.eventSeriesArrayHash.value{2}.eventTimes{1};
Touch_NG = dat.eventSeriesArrayHash.value{2}.eventTimes{2};


T_t = zeros(2,length(W_t));
for i = 1:numel(Touch_G)
    this_t = Touch_G(i);
    [~,ix_t] = min(abs(W_t-this_t)) ;
    T_t(1,ix_t) = 1;
    
end
for i = 1:numel(Touch_NG)
    this_t = Touch_NG(i);
    [~,ix_t] = min(abs(W_t-this_t)) ;
    T_t(2,ix_t) = 1;
    
end

%% Event array (repeated code)
% Events are as follows:
%     'times when the pole was accessible to the whiskers.'
%     'times when the pole was touched by whiskers.'
%     'left lickport contact'
%     'right lickport contact'
%     'left water reward delivery'
%     'right water reward delivery'
%     'auditory cue signaling animal to collect reward'
%     'First touch left'
%     'First touch right'
events = [];
touches = [];

trialStarts = []; dat.trialStartTimes(trials);

for i = 1:numel(trials)
    trialStarts(i) = dat.trialStartTimes(find(dat.trialIds == trials(i)));
    
    
    % Pole up
    x = find(dat.eventSeriesArrayHash.value{1}.eventTrials == trials(i),1,'first');
    events(1,i) = max([0,dat.eventSeriesArrayHash.value{1}.eventTimes(x) - trialStarts(i)]);
    
    % Pole down
    x = find(dat.eventSeriesArrayHash.value{1}.eventTrials == trials(i),1,'last');
    events(2,i) = max([0,dat.eventSeriesArrayHash.value{1}.eventTimes(x) - trialStarts(i)]);
    
    % Lick left
    x = find(dat.eventSeriesArrayHash.value{3}.eventTrials == trials(i),1,'first');
    events(3,i) = max([0,dat.eventSeriesArrayHash.value{3}.eventTimes(x) - trialStarts(i)]);
    
    % Lick right
    x = find(dat.eventSeriesArrayHash.value{4}.eventTrials == trials(i),1,'first');
    events(4,i) = max([0,dat.eventSeriesArrayHash.value{4}.eventTimes(x) - trialStarts(i)]);
    
    % Reward left
    x = find(dat.eventSeriesArrayHash.value{5}.eventTrials == trials(i),1,'first');
    events(5,i) = max([0,dat.eventSeriesArrayHash.value{5}.eventTimes(x) - trialStarts(i)]);
    
    % Reward right
    x = find(dat.eventSeriesArrayHash.value{6}.eventTrials == trials(i),1,'first');
    events(6,i) = max([0,dat.eventSeriesArrayHash.value{6}.eventTimes(x) - trialStarts(i)]);
    
    % reward cue
    x = find(dat.eventSeriesArrayHash.value{7}.eventTrials == trials(i),1,'first');
    events(7,i) = max([0,dat.eventSeriesArrayHash.value{7}.eventTimes(x) - trialStarts(i)]);
    
    x = find(dat.eventSeriesArrayHash.value{2}.eventTrials{1} == trials(i),1,'first');
    events(8,i) = max([0,dat.eventSeriesArrayHash.value{2}.eventTimes{1}(x) - trialStarts(i)]);
    
    x = find(dat.eventSeriesArrayHash.value{2}.eventTrials{2} == trials(i),1,'first');
    events(9,i) = max([0,dat.eventSeriesArrayHash.value{2}.eventTimes{2}(x) - trialStarts(i)]);
    
    
end

% events(:,1) = zeros(9,1);
clf
plot(events');
legend('Pole up','Pole down','Lick left','Lick right','Reward left','Reward right','Reward cue','Touch left','Touch right')

% Medians
for i = 1:size(events,1)
    med_events(i) = median(events(i,find(events(i,:))));
end

%% Distribution plot of each

figure(24); clf
% subplot(2,1,1);
% plot(ca_ts,psth_ca','linewidth',2)
% xlim([0,7000])

bmap = brewermap(9,'Paired');
edges = linspace(0,5000,101);
% subplot(2,1,2);
for i = 1: size(events,1)
    h = histc(events(i,find(events(i,:))),edges);
    stairs(edges,h/sum(h),'linewidth',2,'color',bmap(i,:)); hold all
    %      [F,XI] = ksdensity(events(i,find(events(i,:)))); hold all
    %      plot(XI,F/sum(F),'linewidth',2,'color',bmap(i,:))
    
end

legend('Pole up','Pole down','Lick left','Lick right','Reward left','Reward right','Reward cue','Touch left','Touch right')
xlim([0,7000])
xlabel('Time (ms)')
ylabel('Frequency (normalized trials)')

%% Plot psth with line for event times
for i = 1:size(events,1)
    med_events(i) = median(events(i,find(events(i,:))));
end

%%
figure(25);
clf
bmap2 = brewermap(size(psth_ca,1),'Spectral');
ca_leg = {};
for i = 1:size(psth_ca,1)
    plot(ca_ts,psth_ca(i,:),'linewidth',2,'color',bmap2(i,:))
    ca_leg = {ca_leg{:},['Grp ' num2str(i)]};
    xlim([0,7000])
    hold all
end

for i = 1:numel(med_events)
    plot(med_events(i)*ones(2,1),[0,0.1],'linewidth',2,'color',bmap(i,:))
end

legend(ca_leg{:},'Pole up','Pole down','Lick left','Lick right','Reward left','Reward right','Reward cue','Touch left','Touch right')

%% Times that correspond to each trial


%% Create imaging-rate versions of event + whisker data
A = dat.timeSeriesArrayHash.value{1}.valueMatrix(1:2,:);

K = dat.timeSeriesArrayHash.value{1}.valueMatrix(3:4,:);




%% Create a psth of event + whisker data

x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);
trialStarts = [];
for i = 1:numel(trials)
    trialStarts(i) = dat.trialStartTimes(find(dat.trialIds == i));
end


for i = 1:numel(trials)
    c_ca(i,:) = [pop_mean(find(x == trials(i))),nan(1,100-numel(find(x == trials(i))))];
end

psth_ca(p,:) = nanmean(c_ca(:,1:50),1);



psth_ca_L(p,:) = nanmean(c_ca(L,1:50),1);
psth_ca_R(p,:) = nanmean(c_ca(R,1:50),1);
psth_ca_C(p,:) = nanmean(c_ca(C,1:50),1);
psth_ca_IC(p,:) = nanmean(c_ca(IC,1:50),1);
psth_ca_hit(p,:) = nanmean(c_ca(hit,1:50),1);
psth_ca_CR(p,:) = nanmean(c_ca(CR,1:50),1);
psth_ca_miss(p,:) = nanmean(c_ca(miss,1:50),1);
psth_ca_FA(p,:) = nanmean(c_ca(FA,1:50),1);

%% Correlation timescales of each group?

%% Aligned raster of some kind
% imagesc(Data.A(skewy(Ix)),Data.A(skewy(Ix)))



%% LOCAL NOISE REJECTION FUNCTION
function [A,B,Data,Rejection,egs] = local_NR(data)
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
% Set links below max(A(:))/C to zero
upper_A = triu(A,1);
link_thresh = max(upper_A(:))/pars.C;
A(find(A<link_thresh)) = 0;

% clean-up A, get largest component, and store as basis for all further analysis
% all indices are with reference to Data.A
Data = {};
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);
% get expected distribution of eigenvalues under null model
% [Data.E,Data.D,Vmodel,Data.ExpA] = RndPoissonConfigModel(Data.A,pars.N,pars.C,optionsModel);
[Data.E,Data.D,Vmodel,Data.ExpA] = poissonSparseWCMReal(Data.A,pars.N,pars.C,optionsModel);

% decompose nodes into signal and noise
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model
% find low-dimensional projection
[Data.Dspace,Data.ixpos,Data.Dn,Data.EigEst,Data.Nspace,Data.ixneg,Data.Dneg,Data.NEigEst] = LowDSpace(B,Data.E,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors
% compute dimensions based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order
Data.PosDn = sum(egs > pars.eg_min);
% node rejection within low-dimensional projection
Rejection = NodeRejection(B,Data.E,pars.I,Vmodel,optionsReject); % N.B. also calls LowDSpace function to find projections
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
end