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

load('/Users/mathew/work/PopWorld/197522_rejection/Rejected_an197522_2013_02_18_events_sv_1.mat')
load('Clustering_Results/Clustered_an197522_2013_02_18_events_sv_1.mat')
load('/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_02_18.mat'); % 2013_03_07

%% High skewness neurons
ca = dat.timeSeriesArrayHash.value{2}.valueMatrix;
ev = dat.timeSeriesArrayHash.value{3}.valueMatrix;
S = skewness(ca');
S(isnan(S)) = [];
%% Image data by order of skewness
[~,skew_order] = sort(S);
imagesc(ca(skew_order,:))

%% Plot 10 most and least skewed cells
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
% imagesc(ca(skewy(Ix(1:50)),:)); % 1:50
imagesc(ML_s(skewy(Ix(1:500)),:)); % 1:50
%% Zooming in
figure(10);
imagesc(ca(skewy(Ix(20:29)),:)); % 1:50 %20:29 %13:21
figure(11);
plot(ca(skewy(Ix(20:29)),:)')

%% Plot each group in turn
figure(16); clf
T = dat.timeSeriesArrayHash.value{2}.time;
for i = 1: numel(unique(newG))
    clf
    this_c = find(newG == i)
    plot(T,ca(skewy(Ix(this_c)),:)');
    title(['Cluster ',num2str(i)])
    pause
    
end


%% Get trial boundaries and plot on top
hold all
x = dat.timeSeriesArrayHash.value{1,2}.trial;
trialId = find(diff(x));
plot([T(trialId);T(trialId)],[-0.5*ones(1,numel(trialId));-1.5*ones(1,numel(trialId))],'k')

%% Plot touches
Touch_G = dat.eventSeriesArrayHash.value{2}.eventTimes{1};
Touch_NG = dat.eventSeriesArrayHash.value{2}.eventTimes{2};

plot(Touch_G,-2*ones(1,numel(Touch_G)),'r*')
plot(Touch_NG,-2.1*ones(1,numel(Touch_NG)),'b*')

%% Load sparse event version of data + generate event array
load('/Volumes/05/Peron_2015/Deconvolution_test/ML_peron.mat')
ncells = numel(ML_peron);
nt = numel(ML_peron{1}.fit);
ML_s = zeros(ncells,nt);
for t = 1:ncells
    try
        spikes = 7*ML_peron{t}.spikest1; % ML spike output is in dt space
        
        % Generate dense array of spikes
        spk_dense = zeros(nt,1);
        spk_dense(round(spikes)) = 1;
        ML_s(t,:) = spk_dense;
        %     events = conv(spk_dense,kernel)/sum(kernel);
        %     ML_e(t,:) = events(1:nt);
    end
end

%% PSTHs of groups

%% Aligned raster of some kind
% imagesc(Data.A(skewy(Ix)),Data.A(skewy(Ix)))