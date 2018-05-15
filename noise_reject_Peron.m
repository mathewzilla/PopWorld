% Matlab script to load Peron et al 2015 crcns data and run noise rejection with different methods, saving the outputs
clear all; close all;
% Add Noise_Rejection to path if necessary
if ~exist('NetworkNoiseRejection')
    addpath(genpath('~/work/NetworkNoiseRejection'))
end

fnames = {'Peron_example';'Peron_example_events';'Peron_example_WCM';'Peron_example_events_WCM'};
models = {'Poiss';'Poiss';'WCM';'WCM'};
datasets = [2,3,2,3];

for i = 1:4
fname = fnames{i}
% analysis parameters from Noise_Rejection example repo
pars.N = 100;           % repeats of permutation
% pars.alpha = 0; %0.95; % 0.95; % 0;         % confidence interval on estimate of maxiumum eigenvalue for null model; set to 0 for mean
pars.I = 0;      % interval 
pars.Model = models{i};   % Poiss or 'WCM' . % which null model
pars.C = 100;             % conversion factor for real-valued weights (set=1 for integers, 'all' to use full data range)
pars.eg_min = 1e-2;      % given machine error, what is acceptable as "zero" eigenvalue

% null model options
optionsModel.Expected = 1;    % compute the expectation over the null model graph ensemble? 
optionsModel.NoLoops = 1;     % prevent self-loops in the null model?

% NodeRejection options
optionsReject.Weight = 'linear'; % 'linear' is default
optionsReject.Norm = 'L2';       % L2 is default
optionsReject.Interval = 'CI';

% load data-file
% load '/media/mathew/Data_1/Peron_ssc-2/with_events/an197522/an197522_2013_03_07.mat'
load '/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_03_07.mat'

% get relevant data matrix
data = dat.timeSeriesArrayHash.value{1,datasets(i)}.valueMatrix;
% ev = dat.timeSeriesArrayHash.value{1,3}.valueMatrix;

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
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(A);

% get expected distribution of eigenvalues under null model
switch pars.Model
    case 'Poiss'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = RndPoissonConfigModel(Data.A,pars.N,pars.C,optionsModel);
    case 'WCM'
        [Data.Emodel,diagnostics,Vmodel,Data.ExpA] = WeightedConfigModel(Data.A,pars.N,pars.C,optionsModel);
    otherwise
        error('Unrecognised null model specified')
end

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
save(['Results/Rejected_' fname],'Rejection','Data','Control','pars','optionsModel','optionsReject')

end

%% Simple visualizations of Matrix, expected matrix, retained neurons etc

B = Data.A - Data.ExpA;

figure
ax(1) = subplot(1,3,1); imagesc(Data.A); axis square; title('A')
ax(2) = subplot(1,3,2); imagesc(Data.ExpA); axis square; title('ExpA')
ax(3) = subplot(1,3,3); imagesc(B); axis square; title('B (A - ExpA)')
linkaxes(ax)

%% Eigenvectors
[V,egs] = eig(B);  % eigenspectra of data modularity matrix
egs = diag(egs); % extract vector from diagonal
[egs,ix] = sort(egs,'descend'); % sort eigenvalues into descending order 
V = V(:,ix);  % sort eigenvectors accordingly

%% Plot matrices again, but ordered by projection on to 3 largest eigenvectors
figure
for i = 1:3
[eig_sort,eig_order] = sort(V(:,i),'descend');
ax((i*3)-3+1) = subplot(3,3,(i*3)-3+1); imagesc(Data.A(eig_order,eig_order)); axis square; 
ax((i*3)-3+2) = subplot(3,3,(i*3)-3+2); imagesc(Data.ExpA(eig_order,eig_order)); axis square;
ax((i*3)-3+3) = subplot(3,3,(i*3)-3+3); imagesc(B(eig_order,eig_order)); axis square; 
end

%% B ordered by tsne
tsne_B = tsne(B,'NumDimensions',1);
[tsne_sort,tsne_order] = sort(tsne_B,'descend');
clf
imagesc(B(tsne_order,tsne_order))

%% Image raw data but in tsne order
load '/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_03_07.mat'
% get relevant data matrix
data = dat.timeSeriesArrayHash.value{1,datasets(i)}.valueMatrix;
% clean up
data(find(isnan(data))) = 0;
% only retained cells
data = data(Data.ixRetain,:);
imagesc(data(tsne_order,:))
%% Cluster with Louvain (for speed)
fnames = {'Peron_example';'Peron_example_events';'Peron_example_WCM';'Peron_example_events_WCM'};

for k = 2
    fname = fnames{k}
    blnLabels = 0;      % write node labels? Omit for large networks
    fontsize = 6;
    
    clusterpars.nreps = 100;
    clusterpars.nLouvain = 5;
    
    % load data
    load(['Results/Rejected_' fname])
    
    %% Cluster with louvain method (for speed)
    numConnected = length(Data.ixSignal_Final);
    % With noise rejection
    [Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    
    % Full dataset
    [Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only
    
    % Save results
    
    %% Plot
    
    % Signal subset
    for i=1:numel(Connected.LouvCluster)
        CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
        [H,h,Ix] = plotClusterMap(Data.Asignal_final,CLou,[],[],'S');
        plotorder = Data.ixSignal_Final(Ix);
        title(['Louvain ' num2str(i)]);
        if blnLabels
            % Add node labels
            set(gca,'Ytick',1:numConnected);
            set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        end
        % set(gca,'XTickLabelRotation',90);
        for j = i+1:numel(Connected.LouvCluster)
            CLou2 = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
            Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numConnected);
        end
    end
    
    % Full dataset
    for i=1:numel(Full.LouvCluster)
        CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
        [HL,h,Ix] = plotClusterMap(Data.A,CLou,[],[],'S');
        title(['Full Louvain ' num2str(i)]);
        plotorder = Ix;
        if blnLabels
            % Add node labels
            set(gca,'Ytick',1:numel(Data.ixRetain));
            set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        end
        for j = i+1:numel(Full.LouvCluster)
            CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
            Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numel(Data.ixRetain));
        end
    end
    
    save(['Results/Clustered_' fname],'Full','Connected','clusterpars')
    
    
% end




%% cluster - with noise rejection
% consensus modularity
% [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,Ngrps,~] = allevsplitConTransitive(Data.Aconnected);

% construct new null model
P = Data.ExpA(Data.ixSignal_Final,Data.ixSignal_Final); % extract relevant part of null model

% or make one
% [Signal.Emodel,~,Vmodel,Signal.ExpA] = RndPoissonConfigModel(Data.Aconnected,pars.N,pars.C,optionsModel);
% P = Data.Aconnected - Signal.ExpA;  % modularity matrix using chosen null model
% % find low-dimensional projection
% [~,~,Signal.Dn,~] = LowDSpace(P,Signal.Emodel,pars.alpha); % to just obtain low-dimensional projection

% then cluster
if Data.Dn > 0
    [Connected.QmaxCluster,Connected.Qmax,Connected.ConsCluster,Connected.ConsQ,ctr] = ...
                                        ConsensusCommunityDetect(Data.Asignal_final,P,1+Data.Dn,1+Data.Dn,clusterpars.nreps);
else
    Connected.QmaxCluster = []; Connected.Qmax = 0; Connected.ConsCluster = []; Connected.ConsQ = 0;
end
% Louvain algorithm
[Connected.LouvCluster,Connected.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.Asignal_final,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% cluster - without noise rejection
if Data.Dn > 0
    [Full.QmaxCluster,Full.Qmax,Full.ConsCluster,Full.ConsQ,~] = ...
                                                ConsensusCommunityDetect(Data.A,Data.ExpA,1+Data.Dn,1+Data.Dn);
else
    Full.QmaxCluster = []; Full.Qmax = 0; Full.ConsCluster = []; Full.ConsQ = 0;
end

[Full.LouvCluster,Full.LouvQ,allCn,allIters] = LouvainCommunityUDnondeterm(Data.A,clusterpars.nLouvain,1);  % run 5 times; return 1st level of hierarchy only

%% plot sorted into group order
numConnected = length(Data.ixSignal_Final);

if Data.Dn > 0
    [H,h,Ix] = plotClusterMap(Data.Asignal_final,Connected.ConsCluster,[],[],'S');
    title('Consensus clustering')
    plotorder = Data.ixSignal_Final(Ix);
    
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
    % compare to the Qmax solution at the requested number of groups
    [H,h,Ix] = plotClusterMap(Data.Asignal_final,Connected.QmaxCluster,[],[],'S');
    title('Qmax clustering')
    if blnLabels
        % Add node labelss
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
        % set(gca,'XTickLabelRotation',90);
    end
   
end

for i=1:numel(Connected.LouvCluster)
    CLou = Connected.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [H,h,Ix] = plotClusterMap(Data.Asignal_final,CLou,[],[],'S');
    plotorder = Data.ixSignal_Final(Ix);
    title(['Louvain ' num2str(i)]);
    if blnLabels
        % Add node labels
        set(gca,'Ytick',1:numConnected);
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
    % set(gca,'XTickLabelRotation',90);
    for j = i+1:numel(Connected.LouvCluster)
        CLou2 = Connected.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Connected.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numConnected);
    end
end

%% without noise rejection
if Data.Dn > 0
    [H,h,Ix] = plotClusterMap(Data.A,Full.ConsCluster,[],[],'S');
    title('Consensus clustering of all')
    plotorder = Ix;
    
    if blnLabels
        % Add node labels
        [srt,I] = sort(Full.ConsCluster,'ascend');
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
end

% Louvain algorithm
for i=1:numel(Full.LouvCluster)
    CLou = Full.LouvCluster{i}{1};  % Repeat#, Level of Hierarchy
    [HL,h,Ix] = plotClusterMap(Data.A,CLou,[],[],'S');
    title(['Full Louvain ' num2str(i)]);
    plotorder = Ix;
    if blnLabels
        % Add node labels
        set(gca,'Ytick',1:numel(Data.ixRetain));
        set(gca,'Yticklabel',Data.nodelabels(plotorder,:),'Fontsize',fontsize);
    end
    for j = i+1:numel(Full.LouvCluster)
        CLou2 = Full.LouvCluster{j}{1};  % Repeat#, Level of Hierarchy
        Full.VI_Louvain(i,j) = VIpartitions(CLou,CLou2) ./ log(numel(Data.ixRetain));
    end
end

save(['Results/Clustered_' fname],'Full','Connected','clusterpars')

end