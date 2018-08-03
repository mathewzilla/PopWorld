%% Code to load and plot peron crcns clustering results

%% Re-generate CCons and save

f = dir('197522_rejection/Rejected*');
% f(2) = [];
for i = 1: numel(f) %linspace(2,numel(f),numel(f)/2)
        
            % LOAD EVENTS
            enw = f(i).name;
            % Load rejection data
            load(['197522_rejection/',enw]);%'/Users/mathew/work/PopWorld/197522_rejection/Rejected_an197522_2013_02_18_events_sv_1.mat')
            
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
            save(['~/work/PopWorld/197522_rejection/CCons_',enw(10:end),'.mat'],'CCons','C')
            end

end


%% Look at cluster ID for data vs events versions
% f = dir('example_clustering/*.mat');
f = f = dir('197522_rejection/Rejected*'); %dir('Clustering_Results/Clustered_an197522*');
f(end) = [];
for i = 1:25 %linspace(2,numel(f),numel(f)/2)
    try
        if strcmp(f(i).name(31:36),'events')
            % LOAD EVENTS
            enw = f(i).name;
            %     load(['example_clustering/',f(i-1).name]);
            load(['Clustering_Results/',enw]);
            subplot(2,2,1)
            plot(conv(Full.LouvCluster{1}{1},gausswin(15),'same'));
            title('Peron events')
            enw_2 = [enw(1:30),'data_s',enw(37:45)];
            
            
            
            
            % Load data
            load /Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_02_18.mat
            
                        %     load(['example_clustering/',f(i).name]);
            % LOAD CALCIUM
            load(['Clustering_Results/',enw_2]);
            subplot(2,2,2)
            plot(conv(Full.LouvCluster{1}{1},gausswin(15),'same'));
            title('Calcium')
            suptitle(f(i).name)
            
        end
    end
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