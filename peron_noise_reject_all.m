% peron_noise_reject
%
% Script to apply null model based Noise Rejection to data from Peron et al 2015

%% Big loop to perform Data_Noise_Rejection on all subvolumes in Peron 2015 (Calcium only for now)


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
topdir = ['/media/mathew/Data_1/Peron_ssc-2/ssc-2'];
% topdir = ['/mnt/isilon/Hlab/Mat_Evans/Peron_ssc_events'];
animals = {'an171923';'an194181';'an194672';'an197522';'an198503';'an229716';'an229717';'an229719'};
Rejection_Results = {};
for i = 1:numel(animals)
    %     cd([topdir,'/',animals{i}]);
    files = dir([topdir,'/',animals{i}]);
    for j = 3:numel(files)
        load([topdir,'/',animals{i},'/',files(j).name]);
        for k = 2:numel(s.timeSeriesArrayHash.value)
            
            fname = [files(j).name(1:end-9),'_sv_',num2str(k-1)]
            
            data = s.timeSeriesArrayHash.value{k}.valueMatrix;
            
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

% Summarise into a table
fnames = dir('Results_batch1/');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    
    if any(strfind(fnames(iF).name,'Rejected'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(10:end-4); % strip out 'Rejected' and .mat
        load(['Results/' fnames(iF).name]);
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
save('Results_batch1/Network_Rejection_Table','Network_Rejection_Table');

%% ksdensity and scatter of various variables
% TO Do plot median with label
figure(1); clf;
subplot(2,2,1);
[y,x] = ksdensity(NRstats(:,4));
plot(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(400,8e-3,['Peak = ',num2str(x(mx),'%.3g')])
text(400,7.4e-3,['Median = ',num2str(median(NRstats(:,4)),'%.3g')])
xlabel('Number of retained neurons (Nretain)')
ylabel('Probability density')
% title('Nretain')

subplot(2,2,2);
[y,x] = ksdensity(NRstats(:,5));
plot(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(60,6.5e-2,['Peak = ',num2str(x(mx),'%.3g')])
text(60,6e-2,['Median = ',num2str(median(NRstats(:,5)),'%.3g')])
xlabel('Dimensionality of retained neurons (Npos)')
ylabel('Probability density')
% title('Npos')

subplot(2,2,3);
[y,x] = ksdensity(NRstats(:,6));
plot(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(0.4,4,['Peak = ',num2str(x(mx),'%.3g')])
text(0.4,3.7,['Median = ',num2str(median(NRstats(:,6)),'%.3g')])
xlabel('Proportion of retained neurons (Pretain)')
ylabel('Probability density')
% title('Pretain')

subplot(2,2,4);
[y,x] = ksdensity(NRstats(:,7));
plot(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(0.125,8,['Peak = ',num2str(x(mx),'%.3g')])
text(0.125,7.5,['Median = ',num2str(median(NRstats(:,7)),'%.3g')])
xlabel('Normalized dimensionality (Npos/Nretain)')
ylabel('Probability density')
% title('Normalized dimensionality')

print(gcf,'-depsc','-painters','noiserejection/figs/NRksdensity.eps')
%% Does dataset size or duration affect noise rejection or dimensionality?
figure(2);clf;

colours = varycolor(8);

for n = 1:length(NRstats)
    % N
    subplot(4,2,1);
    plot(NRstats(n,8),NRstats(n,4),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Number of neurons')
    ylabel('Nretain')
    
    subplot(4,2,3);
    plot(NRstats(n,8),NRstats(n,6),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Number of neurons')
    ylabel('Pretain')
    
    subplot(4,2,5);
    plot(NRstats(n,8),NRstats(n,5),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Number of neurons')
    ylabel('Npos')
    
    subplot(4,2,7);
    plot(NRstats(n,8),NRstats(n,7),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Number of neurons')
    ylabel('Npos/Nretain')
    
    % T
    subplot(4,2,2);
    plot(NRstats(n,9),NRstats(n,4),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Duration of recording (frames)')
    ylabel('Nretain')
    
    subplot(4,2,4);
    plot(NRstats(n,9),NRstats(n,6),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Duration of recording (frames)')
    ylabel('Pretain')
    
    subplot(4,2,6);
    plot(NRstats(n,9),NRstats(n,5),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Duration of recording (frames)')
    ylabel('Npos')
    
    subplot(4,2,8);
    plot(NRstats(n,9),NRstats(n,7),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('Duration of recording (frames)')
    ylabel('Npos/Nretain')
end

suptitle('Noise rejection values Vs dataset size')
% print(gcf,'-depsc','-painters','noiserejection/figs/NRvsNT.eps')

%% Correlates with performance?
figure(3); clf;

colours = varycolor(8);

for n = 1:length(NRstats)
    % N
    subplot(4,2,1);
    plot(NRstats(n,10),NRstats(n,4),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(correct)')
    ylabel('Nretain')
    
    subplot(4,2,3);
    plot(NRstats(n,10),NRstats(n,6),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(correct)')
    ylabel('Pretain')
    
    subplot(4,2,5);
    plot(NRstats(n,10),NRstats(n,5),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(correct)')
    ylabel('Npos')
    
    subplot(4,2,7);
    plot(NRstats(n,10),NRstats(n,7),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(correct)')
    ylabel('Npos/Nretain')
    
    % T
    subplot(4,2,2);
    plot(1-NRstats(n,11),NRstats(n,4),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(licking)')
    ylabel('Nretain')
    
    subplot(4,2,4);
    plot(1-NRstats(n,11),NRstats(n,6),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(licking)')
    ylabel('Pretain')
    
    subplot(4,2,6);
    plot(1-NRstats(n,11),NRstats(n,5),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(licking)')
    ylabel('Npos')
    
    subplot(4,2,8);
    plot(1-NRstats(n,11),NRstats(n,7),'.','color',colours(NRstats(n,1),:))
    hold all
    xlabel('P(licking)')
    ylabel('Npos/Nretain')
end

suptitle('Noise rejection values Vs performance')
% print(gcf,'-depsc','-painters','noiserejection/figs/NRvsPerf.eps')

%% Restrict time series analysis to learning sub volumes (one per animal for 4 animals)

a = 2:5;
SV  = 1;
Sess = {[2:10]; [2:14]; [2:13,15:17]; [2:11,13]}; % {[1:10]; [1:14]; [2:13,15:17]; [1:11,13]};
L_sess = [];
for i = 1:4;
    sessions = Sess{i};
    this_m = find(NRstats(:,1) == a(i));
    these_sv = find(NRstats(this_m,3) == 1);
    these_sessions = find(ismember(NRstats(this_m(these_sv),2),sessions));
    L_sess_index = this_m(these_sv(these_sessions));
    L_sess = [L_sess; NRstats(L_sess_index,:)];
end

%% Performance over time
figure(4);clf;
hold all;
for i = 2:5;
    a = find(L_sess(:,1) == i);
    plot(L_sess(a,2),L_sess(a,10),'color',colours(i,:),'linewidth',2);
end
legend(animals{2:5})
xlabel('Training session')
ylabel('P(correct)')

% print(gcf,'-depsc','-painters','noiserejection/figs/Lsess_Perf.eps')

%% Plot per animal Npos, Pretain and Nretain over time
figure(5);
cmap = varycolor(8);
clf
for i = 2:5;
    a = find(L_sess(:,1) == i);
    subplot(3,2,1);
    plot(L_sess(a,2),L_sess(a,5),'color',cmap(i,:),'linewidth',2);
    ylabel('N pos')
    hold all
    axis square;
    subplot(3,2,2);
    plot(L_sess(a,2),L_sess(a,4),'color',cmap(i,:),'linewidth',2);
    ylabel('N retain')
    hold all
    axis square;
    subplot(3,2,3);
    plot(L_sess(a,2),L_sess(a,6),'color',cmap(i,:),'linewidth',2);
    ylabel('P(retain)')
    hold all
    axis square;
    subplot(3,2,4);
    plot(L_sess(a,2),L_sess(a,7),'color',cmap(i,:),'linewidth',2);
    ylabel('Npos/Nretain')
    hold all
    axis square;
    subplot(3,2,5);
    plot(L_sess(a,2),L_sess(a,10),'color',cmap(i,:),'linewidth',2);
    ylabel('Performance')
    xlabel('Training session')
    hold all
    axis square;
    subplot(3,2,6);
    plot(L_sess(a,2),L_sess(a,9),'color',cmap(i,:),'linewidth',2);
    ylabel('Duration of recording')
    xlabel('Training session')
    hold all
    axis square;
end

print(gcf,'-depsc','-painters','noiserejection/figs/Lsess_NRvstime2.eps')
%% Correlates with performance
figure(6);
clf
for i = 2:5;
    a = find(L_sess(:,1) == i);
    % N
    subplot(2,2,1); hold all;
    %     plot(L_sess(a,4),L_sess(a,10),'color',colours(i,:)); hold all
    %     plot(L_sess(a,4),L_sess(a,10),'.','color',colours(i,:))
    %     xlabel('Nretain'); ylabel('P(correct)')
    %     plot(L_sess(a,10),L_sess(a,4),'color',colours(i,:));
    plot(L_sess(a,10),L_sess(a,4),'.','color',colours(i,:),'markersize',20)
    xlabel('P(correct)'); ylabel('Nretain')
    
    subplot(2,2,2); hold all;
    %     plot(L_sess(a,6),L_sess(a,10),'color',colours(i,:)); hold all
    %     plot(L_sess(a,6),L_sess(a,10),'.','color',colours(i,:))
    %     xlabel('Pretain'); ylabel('P(correct)')
    %     plot(L_sess(a,10),L_sess(a,6),'color',colours(i,:));
    plot(L_sess(a,10),L_sess(a,6),'.','color',colours(i,:),'markersize',20)
    xlabel('P(correct)'); ylabel('Pretain')
    
    subplot(2,2,3); hold all;
    %     plot(L_sess(a,5),L_sess(a,10),'color',colours(i,:)); hold all
    %     plot(L_sess(a,5),L_sess(a,10),'.','color',colours(i,:))
    %     xlabel('Npos'); ylabel('P(correct)')
    %     plot(L_sess(a,10),L_sess(a,5),'color',colours(i,:));
    plot(L_sess(a,10),L_sess(a,5),'.','color',colours(i,:),'markersize',20)
    xlabel('P(correct)'); ylabel('Npos')
    
    subplot(2,2,4); hold all;
    %     plot(L_sess(a,7),L_sess(a,10),'color',colours(i,:)); hold all
    %     plot(L_sess(a,7),L_sess(a,10),'.','color',colours(i,:))
    %     xlabel('Npos/Nretain'); ylabel('P(correct)')
    %     plot(L_sess(a,10),L_sess(a,7),'color',colours(i,:));
    plot(L_sess(a,10),L_sess(a,7),'.','color',colours(i,:),'markersize',20)
    xlabel('P(correct)'); ylabel('Npos/Nretain')
    
end

print(gcf,'-depsc','-painters','noiserejection/figs/Lsess_NRvperf2.eps')

%% NR correlates of basic stats for learning sessions

figure(7);clf;

colours = varycolor(8);

for i = 2:5;
    a = find(L_sess(:,1) == i);
    % N
    subplot(4,2,1);
    plot(L_sess(a,8),L_sess(a,4),'.','color',colours(i,:),'markersize',20)
    hold all
    xlabel('Number of neurons')
    ylabel('Nretain')
    
    subplot(4,2,3);
    plot(L_sess(a,8),L_sess(a,6),'.','color',colours(i,:),'markersize',20)
    hold all
    xlabel('Number of neurons')
    ylabel('Pretain')
    
    subplot(4,2,5);
    plot(L_sess(a,8),L_sess(a,5),'.','color',colours(i,:),'markersize',20)
    hold all
    xlabel('Number of neurons')
    ylabel('Npos')
    
    subplot(4,2,7);
    plot(L_sess(a,8),L_sess(a,7),'.','color',colours(i,:),'markersize',20)
    hold all
    xlabel('Number of neurons')
    ylabel('Npos/Nretain')
    
    % T
    subplot(4,2,2);
    %     plot(L_sess(a,9),L_sess(a,4),'color',colours(i,:))
    hold all
    plot(L_sess(a,9),L_sess(a,4),'.','color',colours(i,:),'markersize',20)
    xlabel('Duration of recording (frames)')
    ylabel('Nretain')
    
    subplot(4,2,4);
    %     plot(L_sess(a,9),L_sess(a,6),'color',colours(i,:))
    hold all
    plot(L_sess(a,9),L_sess(a,6),'.','color',colours(i,:),'markersize',20)
    
    xlabel('Duration of recording (frames)')
    ylabel('Pretain')
    
    subplot(4,2,6);
    %     plot(L_sess(a,9),L_sess(a,5),'color',colours(i,:))
    hold all
    plot(L_sess(a,9),L_sess(a,5),'.','color',colours(i,:),'markersize',20)
    
    xlabel('Duration of recording (frames)')
    ylabel('Npos')
    
    subplot(4,2,8);
    %     plot(L_sess(a,9),L_sess(a,7),'color',colours(i,:))
    hold all
    plot(L_sess(a,9),L_sess(a,7),'.','color',colours(i,:),'markersize',20)
    xlabel('Duration of recording (frames)')
    ylabel('Npos/Nretain')
end

suptitle('Noise rejection values Vs dataset size')
% print(gcf,'-depsc','-painters','noiserejection/figs/NRvsNT.eps')
print(gcf,'-depsc','-painters','noiserejection/figs/Lsess_NRvNT2.eps')

%% 3 space and time
figure(8);clf;

colours = varycolor(8);

for i = 2:5;
    a = find(L_sess(:,1) == i);
    % N + T
    subplot(2,2,1);
    plot3(L_sess(a,9),L_sess(a,8),L_sess(a,4),'.','color',colours(i,:),'markersize',20)
    hold all
    plot3(L_sess(a,9),L_sess(a,8),zeros(numel(a),1),'.','color',[.5,.5,.5])
    xlabel('Duration of recording (frames)')
    ylabel('Number of neurons')
    zlabel('Nretain')
    axis square; grid on; box off
    
    subplot(2,2,2);
    plot3(L_sess(a,9),L_sess(a,8),L_sess(a,6),'.','color',colours(i,:),'markersize',20)
    hold all
    plot3(L_sess(a,9),L_sess(a,8),zeros(numel(a),1),'.','color',[.5,.5,.5])
    xlabel('Duration of recording (frames)')
    ylabel('Number of neurons')
    zlabel('Pretain')
    axis square; grid on; box off
    
    subplot(2,2,3);
    plot3(L_sess(a,9),L_sess(a,8),L_sess(a,5),'.','color',colours(i,:),'markersize',20)
    hold all
    plot3(L_sess(a,9),L_sess(a,8),zeros(numel(a),1),'.','color',[.5,.5,.5])
    xlabel('Duration of recording (frames)')
    ylabel('Number of neurons')
    zlabel('Npos')
    axis square; grid on; box off
    
    subplot(2,2,4);
    plot3(L_sess(a,9),L_sess(a,8),L_sess(a,7),'.','color',colours(i,:),'markersize',20)
    hold all
    plot3(L_sess(a,9),L_sess(a,8),zeros(numel(a),1),'.','color',[.5,.5,.5])
    xlabel('Duration of recording (frames)')
    ylabel('Number of neurons')
    zlabel('Npos/Nretain')
    axis square; grid on; box off
end

print(gcf,'-depsc','-painters','noiserejection/figs/Lsess_NRvNT3D.eps')
