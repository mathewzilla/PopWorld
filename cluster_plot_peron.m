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


%% Load typical example data + clustering output


% load('/Users/mathew/work/PopWorld/Clustering_Results/Clustered_an197522_2013_03_07_data_s_sv_1.mat')
% load('/Users/mathew/work/PopWorld/Results/Clustered_Peron_example.mat')
% load('/Users/mathew/work/PopWorld/Results/Rejected_Peron_example.mat')

% load('/Users/mathew/work/PopWorld/197522_rejection/Rejected_an197522_2013_02_18_events_sv_1.mat')
% load('Clustering_Results/Clustered_an197522_2013_02_18_events_sv_1.mat')
% load('/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_02_18.mat'); % 2013_03_07

clear all
Results_path = '/Volumes/Extras/197522_rejection/Rejection_tests/';
fname = 'Pre_conv_1K_20_02'; % 'Pre_conv_1K_21_02';

load([Results_path,'Scaling_Results_Table'])
load([Results_path,'/Rejected_',fname],'Data','Rejection')
load([Results_path,'/Clustered_',fname],'Full','Connected')
load('/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_02_20.mat')

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
A = dat.timeSeriesArrayHash.value{1}.valueMatrix(1:2,:);

K = dat.timeSeriesArrayHash.value{1}.valueMatrix(3:4,:);

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

%% Calcium times
x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

trialStarts = [];
for i = 1:numel(trials)
    trialStarts(i) = dat.trialStartTimes(find(dat.trialIds == trials(i)));
    ca_t(i,:) = [dat.timeSeriesArrayHash.value{2}.time(x==trials(i))-trialStarts(i),nan(1,100-numel(find(x == trials(i))))];
end

ca_ts = nanmean(ca_t(:,1:50));
%% Plot psths of all groups

% loop to plot event lines on top (pole up,down,reward port)

for f = [20,21,22]
    figure(f);clf;
    if f ~= 20
        for j = 1:4
            subplot(2,2,j)
            event_subset = [1,2,7];
            for i = 1:numel(event_subset)
                plot(med_events(event_subset(i))*ones(2,1),[0,0.1],'linewidth',2,'color','k');
                hold all
            end
        end
    else
        event_subset = [1,2,7];
        for i = 1:numel(event_subset)
            plot(med_events(event_subset(i))*ones(2,1),[0,0.1],'linewidth',2,'color','k');
            hold all
        end
    end
    
end


figure(20); 
plot(ca_ts,psth_ca','linewidth',2)
title('Population PSTH')
 xlim([0,7000])
legend('Pole up','Pole down','Reward cue')
xlabel('Time (ms)')
ylabel('Cluster mean Ca2+')
ylim([0,0.1])

% Two 2x2 plots of the trial type/ performance data
figure(21); 
ax(1) = subplot(2,2,1);
plot(ca_ts,psth_ca_L','linewidth',2); title('Left')
ax(2) = subplot(2,2,2);
plot(ca_ts,psth_ca_R','linewidth',2); title('Right')
ax(3) = subplot(2,2,3);
plot(ca_ts,psth_ca_C','linewidth',2); title('Correct')
xlabel('Time (ms)')
ylabel('Cluster mean Ca2+')
ax(4) = subplot(2,2,4);
plot(ca_ts,psth_ca_IC','linewidth',2); title('Incorrect')
linkaxes(ax)
xlim([0,7000])
ylim([0,0.1])

figure(22); 
bx(1) = subplot(2,2,1);
plot(ca_ts,psth_ca_hit','linewidth',2); title('Correct Left')
bx(2) = subplot(2,2,2);
plot(ca_ts,psth_ca_CR','linewidth',2); title('Correct Right')
bx(3) = subplot(2,2,3);
plot(ca_ts,psth_ca_miss','linewidth',2); title('Incorrect Left')
xlabel('Time (ms)')
ylabel('Cluster mean Ca2+')
bx(4) = subplot(2,2,4);
plot(ca_ts,psth_ca_FA','linewidth',2); title('Incorrect Right')
linkaxes(bx)
xlim([0,7000])
ylim([0,0.1])

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

%% Event array
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