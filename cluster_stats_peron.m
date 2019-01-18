% cluster_stats_peron.m
% 
% Simple script to load clustering data into a table, then plot
% distributions etc

%% gather all results into a Table...
clear all; close all;

fnames = dir('Clustering_Results_preround/Clustered*');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    try
    if any(strfind(fnames(iF).name,'Clustered'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(11:end-4); % strip out 'Clustered' and .mat
        load(['Clustering_Results_preround/' fnames(iF).name]);
        % keyboard
        result(netCtr).Signal_Consensus_Grps = numel(unique(Connected.ConsCluster)); 
        result(netCtr).Signal_Qmax_Grps = numel(unique(Connected.QmaxCluster)); 
        n = cellfun(@(x) numel(unique(x{1})),Connected.LouvCluster);    % number of groups in each Louvain clustering
        result(netCtr).Signal_Louvain_MeanGrps = mean(n); 
        result(netCtr).Signal_Louvain_RangeGrps = range(n); 
        if isfield(Connected,'VI_Louvain')
            ix = find(~tril(ones(size(Connected.VI_Louvain))));
            result(netCtr).Signal_Louvain_MeanVI = mean(Connected.VI_Louvain(ix));
        end
        result(netCtr).Full_Consensus_Grps = numel(unique(Full.ConsCluster)); 
        result(netCtr).Full_Qmax_Grps = numel(unique(Full.QmaxCluster)); 
        n = cellfun(@(x) numel(unique(x{1})),Full.LouvCluster);    % number of groups in each Louvain clustering
        result(netCtr).Full_Louvain_MeanGrps = mean(n); 
        result(netCtr).Full_Louvain_RangeGrps = range(n); 
        if isfield(Full,'VI_Louvain')
            ix = find(~tril(ones(size(Full.VI_Louvain))));
            result(netCtr).Full_Louvain_MeanVI = mean(Full.VI_Louvain(ix));
        end
    end
    catch
        display(['Problem with ',fnames(iF).name])
    end
end

Network_Clustering_Table = struct2table(result);
save('Clustering_Results_preround/Network_Clustering_Table','Network_Clustering_Table');

%% Load Network_Clustering_Table and Rejection table
% load('Results_reject_preround/Network_Rejection_Table_wStats_preround3.mat')
load('Results_reject_preround/Network_Rejection_Table_wEvents.mat')
load('Clustering_Results_preround/Network_Clustering_Table.mat')

%% Combine Clustering and Rejection results into a single table
clear result
n = 0;
for i = 1:height(Network_Clustering_Table)
    if Network_Clustering_Table.Signal_Qmax_Grps(i)>0
        n = n+1;
        result(n).NetworkName = Network_Clustering_Table.NetworkName{i};
        result(n).Signal_Consensus_Grps = Network_Clustering_Table.Signal_Consensus_Grps(i);
        result(n).Signal_Qmax_Grps = Network_Clustering_Table.Signal_Qmax_Grps(i);
        result(n).Signal_Louvain_MeanGrps = Network_Clustering_Table.Signal_Louvain_MeanGrps(i);
        result(n).Full_Consensus_Grps = Network_Clustering_Table.Full_Consensus_Grps(i);
        result(n).Full_Qmax_Grps = Network_Clustering_Table.Full_Qmax_Grps(i);
        result(n).Full_Louvain_MeanGrps = Network_Clustering_Table.Full_Louvain_MeanGrps(i);
        
        for j = 1:height(Network_Rejection_Table)
            if strcmp(Network_Rejection_Table.NetworkName{j},Network_Clustering_Table.NetworkName{i})
                result(n).Network_Size = Network_Rejection_Table.Network_Size(j);
                result(n).Signal_Size_WCM = Network_Rejection_Table.Signal_Size_WCM(j);
                result(n).WCM_Dn = Network_Rejection_Table.WCM_Dn(j);
                result(n).WCM_RejectionDn = Network_Rejection_Table.WCM_RejectionDn(j);
                result(n).Animal = Network_Rejection_Table.Animal(j);
                result(n).Session = Network_Rejection_Table.Session(j);
                result(n).Subvolume = Network_Rejection_Table.Subvolume(j);
                result(n).N = Network_Rejection_Table.N(j);
                result(n).T = Network_Rejection_Table.T(j);
                result(n).Pcorrect = Network_Rejection_Table.Pcorrect(j);
                result(n).method = Network_Rejection_Table.method(j);
                result(n).eig90 = Network_Rejection_Table.eig90(j);
                result(n).SVID_N = Network_Rejection_Table.SVID_N(j);
                result(n).sv_unique = Network_Rejection_Table.sv_unique(j);
                result(n).sv_order = Network_Rejection_Table.sv_order(j);
                result(n).Learning = Network_Rejection_Table.Learning(j);
            end
        end
    end
end

Network_Combined_Table = struct2table(result);

%% Simple animal based colours, with no colour indicating learning session
% Dots: black/grey based on learning/not learning
% Edges - 8 animals

colours = varycolor(8);

% Dot colour based on method
dotcolour = zeros(height(Network_Combined_Table),3);
edgecolour = zeros(height(Network_Combined_Table),3);
colour_ID = ones(height(Network_Combined_Table),1);

for n = 1:height(Network_Combined_Table)
    dotcolour(n,:) = [.5,.5,.5];
    edgecolour(n,:) = colours(Network_Combined_Table.Animal(n),:);
    colour_ID(n) = Network_Combined_Table.Animal(n);
end

%% Plot Eig90 vs Retained D vs Consensus
var1 = Network_Combined_Table.WCM_RejectionDn;
var2 = Network_Combined_Table.eig90;
var3 = Network_Combined_Table.Signal_Consensus_Grps;

labels = {'D_{rejection}';'Eig_{90}';'Consensus'};

% plot_rejection_pairs(va1,va2,labels);
figure(1); clf;
subplot(1,3,1)
labels = {'D_{rejection}';'Eig_{90}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,3,2)
labels = {'D_{rejection}';'Consensus'};
h(1) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

subplot(1,3,3)
labels = {'Eig_{90}';'Consensus'};
h(1) = plot_rejection_pairs(var2,var3,labels,colour_ID,edgecolour,dotcolour)






















%% Local plotting functions (required Matlab 2016a or later)

%% Stripped down single plot figure
function h = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)
cla; hold all
C = unique(colour_ID);
for c = 1:numel(C)
    these_d = find(colour_ID == C(c));
    h = plot(var1(these_d),var2(these_d),'o','markeredgecolor',edgecolour(these_d(1),:),'markerfacecolor',dotcolour(these_d(1),:),'markersize',5);
end
axis square
xlabel(labels{1});
ylabel(labels{2});
end
