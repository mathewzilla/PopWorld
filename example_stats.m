3% example_stats.m
% Script to compute stats from example data, and plot them
% e.g. Event rate, N silent cells, N dimensions 

cd ~/work/Peron_crcns/
% load('/mnt/isilon/Hlab/Mat_Evans/Peron_ssc_events/an197522/an197522_2013_03_07.mat')
load('/Volumes/05/Peron_2015/Peron_ssc_events/an197522/an197522_2013_03_07.mat');
% load('/media/mathew/Data_1/Peron_ssc-2/with_events/an197522/an197522_2013_03_07.mat')

local_save_path = '/media/mathew/Data_1/Peron_shuffles/compile_test_gather/';
data_path = '/Volumes/05/Peron_2015/Deconvolution_test/example_data/';
% data_path = '~/work/Peron_crcns/example_data/';

data_ID = {'ca';'ev';'LZ_e';'ML_e';'S2P_e';'Y_e';'LZ_t';'S2P_t10';'ML_t';'ML_p';'LZ_k';'LZ_t2';'S2P_k10';'S2P_k6';'ML_t_hand';'ML_e2';'S2P_t6'}; %{'LZ_t';'S2P_t10';'S2P_t6';'ML_t';'ML_p'}; %{'ca';'ev';'LZ_e';'ML_e';'S2P_e';'Y_e'};
% 'LZ_k';'LZ_t2';'S2P_k10';'S2P_k6';'ML_t_hand';'ML_e2';'S2P_t6'
methods = {'ca';'ev';'LZ';'ML';'S2P';'Y';'LZ_t';'S2P_t10';'ML_t';'ML_p';'LZ_k';'LZ_t2';'S2P_k10';'S2P_k6';'ML_t_hand';'ML_e2';'S2P_t6'}; % {'LZ_t';'S2P_t10';'S2P_t6';'ML_t';'ML_p'}; %{'ca';'ev';'LZ';'ML';'S2P';'Y'};
meth_names = {'Calcium','Peron','LZero_{kernel_{old}}','MLSpike_{kernel_{hand}}','Suite2P_{kernel_{hand}}','Yaksi','LZero_{events_{hand}}','Suite2P_{events_{PCC}}','MLSpike_{events_{ER}}','MLSpike_{pspike}','LZero_{kernel_{ER}}','LZero_{events_{ER}}','Suite2P_{kernel_{PCC}}','Suite2P_{kernel_{ER}}','MLSpike_{events_{hand}}','MLSpike_{kernel_{ER}}','Suite2P_{events_{ER}}'};

nmeth = numel(methods);
x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

[ncells,nt] = size(dat.timeSeriesArrayHash.value{1,2}.valueMatrix);

L = find(dat.trialTypeMat(1,trials));
R = find(dat.trialTypeMat(2,trials));

% Load event rates for best methods
gmeths = [1,2,6,13,4,11,8,15,12];
meths = gmeths; nmeths = numel(meths);

load example_data/event_rates_best.mat
%% Eigenvalues and variance explained

clear egs explained

for j = 1:nmeths % 1:6;
    meth = methods{meths(j)}
    % Load data
    load([data_path,data_ID{meths(j)},'.mat'])
%     load(['~/work/Peron_crcns/example_data/',data_ID{j},'.mat'])
    data = eval(data_ID{meths(j)});
    data(find(isnan(data))) = 0;

    [V,D] = eig(cov(data'));
    
    
    egs(j,:) = sort(diag(D),1,'descend');
    explained(j,:) = 100*cumsum(egs(j,:))/sum(egs(j,:));
end

% save example_eigs/example_eigs egs explained
%% Plot
figure(1);
clf
colours = distinguishable_colors(nmeth);
for i = 1:nmeths
    plot(explained(i,:),'color',colours(meths(i),:),'linewidth',2);
    hold all
end
ylabel('Variance explained (%)')
xlabel('N eigenvectors')
legend(meth_names(meths))

% Add line for 80% variance
plot([0,ncells],[80,80],'k--')

% For each method, work out nDims for >80 explained variance
for i = 1:nmeths
    th(i) = find(explained(i,:) >= 80,1,'first');
    plot(th(i)*ones(1,2),[0,80],'color',colours(meths(i),:))
end

%% Correlation matrices + distributions
clf
CXY = {};
for j = 1:nmeths % 1:6;
    meth = methods{meths(j)}
    % Load data
    load([data_path,data_ID{meths(j)},'.mat'])
%     load(['~/work/Peron_crcns/example_data/',data_ID{j},'.mat'])
    data = eval(data_ID{meths(j)});
    
    pcc_data = corrcoef(data');
    pcc_data(find(eye(ncells))) = 0;
    pcc_data(find(isnan(pcc_data))) = 0;
    CXY{j} = pcc_data;
    ax(j) = subplot(3,3,j);
    imagesc(pcc_data)
    title(meth_names{meths(j)})
    axis square
    axis off
end

xlim([405,425])
ylim([405,425])

linkaxes(ax)

%% Ksdensity plot - first 6 methods
figure(12); clf;
n = 6; for n = 1:6; plot(-0.5,0,'linewidth',2,'color',cmap(meths(n),:)); hold all; end
plot([0,0],[0,80],'k--')
n = 0;
for n = 1:6; % 1:nmeth
    meth_names{meths(n)}
%     n = n+1;

    T = triu(CXY{n},1);
    pcc_1 = T(find(T));

    % subplot(2,1,2)
    figure(12);
    hold all
    [h,xi] = ksdensity(pcc_1);
    plot(xi,h,'linewidth',2,'color',cmap(meths(n),:));
    plot([prctile(pcc_1,[5,50,95]);prctile(pcc_1,[5,50,95])],[60+((5-n)*4)-4,60+((5-n)*4);60+((5-n)*4)-4,60+((5-n)*4);60+((5-n)*4)-4,60+((5-n)*4)]','linewidth',2,'color',cmap(meths(n),:))
    plot(prctile(pcc_1,[5,95]),[60+((5-n)*4)-2,60+((5-n)*4)-2],'linewidth',2,'color',cmap(meths(n),:))
    
    
    
end

figure(12);
xlim([-0.25,0.25])
legend(meth_names(meths(1:6)),'location','best')

%% Ksdensity of CXY - spike inference
figure(13); clf;
for n = 7:9; plot(-0.5,0,'linewidth',2,'color',cmap(meths(n),:)); hold all; end
plot([0,0],[0,85],'k--')
n = 0;
for n = 7:9; % 1:nmeth
    meth_names{meths(n)}
%     n = n+1;

    T = triu(CXY{n},1);
    pcc_1 = T(find(T));

    % subplot(2,1,2)
    figure(13);
    hold all
    [h,xi] = ksdensity(pcc_1);
    plot(xi,h,'linewidth',2,'color',cmap(meths(n),:));
    plot([prctile(pcc_1,[5,50,95]);prctile(pcc_1,[5,50,95])],[90+((5-n)*4)-4,90+((5-n)*4);90+((5-n)*4)-4,90+((5-n)*4);90+((5-n)*4)-4,90+((5-n)*4)]','linewidth',2,'color',cmap(meths(n),:))
    plot(prctile(pcc_1,[5,95]),[90+((5-n)*4)-2,90+((5-n)*4)-2],'linewidth',2,'color',cmap(meths(n),:))
    
    
    
end

figure(13);
xlim([-0.25,0.25])
legend(meth_names(meths(7:9)))

%% CXY data in a single array
X = ones(ncells);
Y = triu(X);
Ti = find(Y);

CXY_all = [];


for n = 1:9; % 1:nmeth
    meth_names{meths(n)}
%     n = n+1;
    T = CXY{n}(Ti);
    CXY_all(n,:) = T;
end

%% Scatter - all 9 methods
siz = 5;
alph = 0.1;
figure(14); clf;
n = 9; for n = 1:9; plot(-0.5,0,'linewidth',2,'color',cmap(meths(n),:)); hold all; end
plot([0,0],[0,80],'k--')
n = 0;


for n = 1:9;
    figure(14);
    hold all
    
    npairs = numel(T);
    % plot(0.1*randn(npairs,1)+ones(npairs,1),pcc_1,'.','markersize',5,'color',cmap(1,:))
    scatter(0.1*randn(npairs,1)+n+1*ones(npairs,1),CXY_all(n,:),siz,cmap(meths(n),:),'filled');%,'markerfacealpha',alph,'markeredgealpha',alph)
    plot(n+[0.75,1.25;0.75,1.25;0.75,1.25]',[prctile(CXY_all(n,:),[5,50,95]);prctile(CXY_all(n,:),[5,50,95])],'k','linewidth',2)
    
end
figure(14);
ylim([-0.4,1]);
xlim([1.25,10.75])
plot([0,11],[0,0],'k--')
set(gca,'XTick',2:n+1,'XTickLabel',meth_names(meths),'XTickLabelrotation',45)

legend(meth_names(meths),'location','bestoutside')

%% Boxplot - all 9 methods
clf
h = boxplot(CXY_all','boxstyle','outline','colors',cmap(meths,:),'orientation','horizontal','symbol','o','notch','on')
set(gca, 'YTick', [1:numel(gmeths)]);
set(gca,'Ydir','reverse')
set(h,{'linew'},{2})
% set(gca, 'XTickLabels',meth_names(gmeths),'XTickLabelrotation',45);% {'Calcium','Peron','LZero','MLSpike','Suite2P','Yaksi','ML_evnt','S2P_evnt','LZ_evnt'})
xlabel('Pairwise correlation coefficient')
hold all
plot([0,0],[0,10],'k--')

% Proper tick labels
set(gca, 'yticklabel', []) %Remove tick labels
% Get tick mark positions
yTicks = get(gca, 'ytick');
xTicks = get(gca, 'xtick');
ax = axis; %Get left most x-position
HorizontalOffset = 0.1;
% Reset the ytick labels in desired font
for i = 1:length(yTicks)
%Create text box and set appropriate properties
     text(ax(1) - HorizontalOffset,yTicks(i),meth_names(gmeths(i)),...
         'HorizontalAlignment','Right');   
end


%% CXY of CXYs, all 9 then separate for spike inference and first 6 methods
cxycxy = zeros(9,9);
for j = 1:9;
    for k = 1:9;
        if j == k
            
        else
            pcc = corrcoef(CXY_all(j,:),CXY_all(k,:));
            cxycxy(j,k) = pcc(2,1);
        end
    end
end
clf
imagesc(cxycxy)
axis square
set(gca,'XTick',1:9,'XTickLabel',meth_names(meths),'XTickLabelrotation',45)
set(gca,'YTick',1:9,'YTickLabel',meth_names(meths))

figure(2);
subplot(1,2,1);
imagesc(cxycxy(1:6,1:6))
axis square
set(gca,'XTick',1:6,'XTickLabel',meth_names(meths(1:6)),'XTickLabelrotation',45)
set(gca,'YTick',1:6,'YTickLabel',meth_names(meths(1:6)))

subplot(1,2,2);
imagesc(cxycxy(7:9,7:9))
axis square
set(gca,'XTick',1:3,'XTickLabel',meth_names(meths(7:9)),'XTickLabelrotation',45)
set(gca,'YTick',1:3,'YTickLabel',meth_names(meths(7:9)))
%% Firing rate
% All usable methods, ordered by approach (Orig, kernel_best, events_best, kernel_strawman, events_strawman)
% Orig = 1,2,6 (Calcium, Peron, Yaksi)
% Kernel = 3,4,5,11,13,14,16 (LZ_old, ML_hand, S2P_hand, LZ_ER, S2P_PCC, S2P_ER, ML_ER)
% Events = 7,8,9,12,15,17 (LZ_hand, S2P_PCC, ML_ER, LZ_ER, ML_hand, S2P_ER)

%% Event rate of kernel/orig version of data
kmeths = [1,2,6,3,4,5,11,13,14,16];
emeths = [7,8,9,12,15,17];

% Good methods to use in the paper
gmeths = [1,2,6,13,4,11,8,15,12];

meths = kmeths;
nmeths = numel(meths);
clear ER_k
for j = 1:nmeths; %nmeths;
    j
%     load(['~/work/Peron_crcns/example_data/',data_ID{j},'.mat'])
    load([data_path,data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    
    for c = 1:ncells
        c;
        this_c = data(c,:);

        % Threshold based on std of residual noise (remove smooth version of data
        % first)
        smooth_c = conv(this_c,ones(4,1),'same')/4;
        resid = this_c - smooth_c;
        
        % Clear zeros
        resid = resid(find(resid));
        
        sig = std(resid);
        level = mean(this_c) + 3*sig;
        
        % Plot
        % clf
        % plot(this_c)
        % hold all
        % plot([0,nt],level*ones(2,1))
        
        found_events = find(this_c>=level);
        
        % Event rate (Hz)
        ER_k(j,c) = 7 * numel(found_events)/nt;
    end
end

%% Event rate for spike inference methods
meths = emeths;
nmeths = numel(meths);
clear ER_e
for j = 1:nmeths;
    j
    load([data_path,data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    for c = 1:ncells;
        this_c = data(c,:);

        ER_e(j,c) = 7 * (numel(find(this_c))/nt);
    end
end

%% Plot (kernel methods)

cmap = distinguishable_colors(numel(methods));
figure(2); clf;
mx = max(ER_k(:));
mn = min(ER_k(:));

meths = kmeths;
nmeths = numel(meths);

bins = linspace(mn,mx,100);
clear ax
for j = 1:nmeths;
    ax(j) = subplot(nmeths,1,j)
    [F] = hist(ER_k(j,:),bins);
    
    bar(bins,F,'facecolor',cmap(meths(j),:),'edgecolor',cmap(meths(j),:))
    ylabel(meth_names{meths(j)})
    set(get(gca,'ylabel'),'rotation',0, 'HorizontalAlignment','right')
%     axis off
end
linkaxes(ax,'x')
xlabel('Event rate (Hz)')

%% Plot (event methods)
%% Plot

cmap = distinguishable_colors(numel(methods));
figure(3); clf
mx = max(ER_k(:));
mn = min(ER_k(:));

meths = emeths;
nmeths = numel(meths);

bins = linspace(mn,mx,100);
clear ax
for j = 1:nmeths;
    ax(j) = subplot(nmeths,1,j)
    [F] = hist(ER_e(j,:),bins);
    
    bar(bins,F,'facecolor',cmap(meths(j),:),'edgecolor',cmap(meths(j),:))
    ylabel(meth_names{meths(j)})
    set(get(gca,'ylabel'),'rotation',0, 'HorizontalAlignment','right')
%     axis off
end
linkaxes(ax,'x')
xlabel('Event rate (Hz)')

%% Reduce list to 'best' methods

%% Event rate of kernel/orig version of data
kmeths = [1,2,6,3,4,5,11,13,14,16];
emeths = [7,8,9,12,15,17];

% Good methods to use in the paper
gmeths = [1,2,6,13,4,11,8,15,12];

meths = gmeths;
nmeths = numel(meths);
clear ER_g
for j = 1:nmeths; %nmeths;
    j
%     load(['~/work/Peron_crcns/example_data/',data_ID{j},'.mat'])
    load([data_path,data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    
    
    if ismember(meths(j),kmeths); % deconvolved calcium methods
        for c = 1:ncells
            c;
            this_c = data(c,:);
            
            % Threshold based on std of residual noise (remove smooth version of data
            % first)
            smooth_c = conv(this_c,ones(4,1),'same')/4;
            resid = this_c - smooth_c;
            
            % Clear zeros
            resid = resid(find(resid));
            
            sig = std(resid);
            level = mean(this_c) + 3*sig;
            
            % Plot
            % clf
            % plot(this_c)
            % hold all
            % plot([0,nt],level*ones(2,1))
            
            found_events = find(this_c>=level);
            
            % Event rate (Hz)
            ER_g(j,c) = 7 * numel(found_events)/nt;
        end
        
    elseif ismember(meths(j),emeths); % Spike inference methods
        
        for c = 1:ncells;
            this_c = data(c,:);
            
            ER_g(j,c) = 7 * (numel(find(this_c))/nt);
        end
    end
end


% Fix buggy LZero result where all elements are 1 instead of 0 for silent cells
silly_cells = find(ER_g(9,:)>3.5);
ER_g(9,silly_cells) = 0;

% save example_data/event_rates_best.mat ER_g
%% Plot results for best methods

cmap = distinguishable_colors(numel(methods));
figure(4); clf
mx = max(ER_k(:));
mn = min(ER_k(:));

meths = gmeths;
nmeths = numel(meths);

bins = linspace(mn,mx,100);
clear ax
for j = 1:nmeths;
    ax(j) = subplot(nmeths,1,j)
    [F] = hist(ER_g(j,:),bins);
    
    bar(bins,F,'facecolor',cmap(meths(j),:),'edgecolor',cmap(meths(j),:))
    ylabel(meth_names{meths(j)})
    set(get(gca,'ylabel'),'rotation',0, 'HorizontalAlignment','right')
%     axis off
end
linkaxes(ax,'x')
xlabel('Event rate (Hz)')

%% Event rate boxplot for all methods
figure(5);clf;
% boxplot(ER_g')
h = boxplot(ER_g','boxstyle','outline','colors',cmap(meths,:),'orientation','horizontal','symbol','o','notch','on')
set(gca, 'YTick', [1:numel(gmeths)]);
set(gca,'Ydir','reverse')
set(h,{'linew'},{2})
% set(gca, 'XTickLabels',meth_names(gmeths),'XTickLabelrotation',45);% {'Calcium','Peron','LZero','MLSpike','Suite2P','Yaksi','ML_evnt','S2P_evnt','LZ_evnt'})
xlabel('Event rate (Hz)')

% Proper tick labels
% Generate figure and remove ticklabels
% close all;
% plot(1:numel(gmeths));
set(gca, 'yticklabel', []) %Remove tick labels
% Get tick mark positions
yTicks = get(gca, 'ytick');
xTicks = get(gca, 'xtick');
ax = axis; %Get left most x-position

HorizontalOffset = 0.1;
% Reset the ytick labels in desired font
for i = 1:length(yTicks)
%Create text box and set appropriate properties
     text(ax(1) - HorizontalOffset,yTicks(i),meth_names(gmeths(i)),...
         'HorizontalAlignment','Right');   
end

% % Reset the xtick labels in desired format 
% minX = min(xTicks);
% minY = min(yTicks);
% verticalOffset = 0.2;
% for xx = 1:length(yTicks)
% %Create text box and set appropriate properties
% %      text(xTicks(xx), minY - verticalOffset, meth_names(gmeths(xx)),...
% %          'HorizontalAlignment','right','Rotation',45);   
%     text(yTicks(xx),minX - verticalOffset,meth_names(gmeths(xx)),...
%         'VerticalAlignment','top')
% end

%% Proportion of silent cells (using Peron 0.0083Hz definition) 
% - O'Connor 2010 found ~13% silent cells across cortical layers. 
% 44 of cells had FRs < 1 Hz. L2/3 this figure was closer to 60%
% L2/3 - mean 3Hz, median 0.18Hz. ~26% silent cells
figure(6);
clf;
for j = 1:nmeths
%     subplot(2,ceil(nmeths/2),j)
    subplot(3,3,j)
    not_silent = zeros(ncells,1);
    not_silent(find(ER_g(j,:)>0.0083)) = 1;
    
    pie(hist(not_silent,2))
    colormap([0,0,0;.5,.5,.5])
    title(meth_names{gmeths(j)})
end
axis off
 suptitle('Silent (black) vs active (gray) cells')

%% Touch tuning - touch triggered average + stats test on data distribution vs shuffle

% Shift Yaksi to address causality problem
%Y_ey = circshift(Y_ey,[0,1]);

x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

% Find the first touch on every trial
first_touch = [];
for i = 1:numel(trials)
    tt = [];
    x1 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{1} ==i,1,'first');
    tt1 = dat.eventSeriesArrayHash.value{2}.eventTimes{1}(x1);
    x2 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{2} ==i,1,'first');
    tt2 = dat.eventSeriesArrayHash.value{2}.eventTimes{2}(x2);
    tt = [tt1,tt2];
    if find(tt);
        first_touch = [first_touch, min(tt(find(tt)))];
    end
    
end

% Closest frame to touch time (as touch is measured in milliseconds)
ftf = [];
for i = 1:numel(first_touch)
    [mn,im] = min(abs(dat.timeSeriesArrayHash.value{2}.time - first_touch(i)));
    ftf(i) = im; % dat.timeSeriesArrayHash.value{2}.time(im);
end

TTA = {};
TT_tuned = {};
for j = 1:nmeths
    load([data_path,data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    tth = zeros(ncells,15);
    tt_t = zeros(ncells,1);
    tt_P = zeros(ncells,1);
    
    for c = 1:ncells
        this_c = data(c,:);
        tt_sample = zeros(numel(ftf),15);
        
        % Equivalently sized array of random segments
        tt_shuff = zeros(numel(ftf),15);
        rsamp = randsample(1:nt-15,numel(ftf));
        
        for i = 1:numel(ftf);
            tt_sample(i,:) = this_c(ftf(i)-7 : ftf(i)+7);
            % Random sample
            tt_shuff(i,:) = this_c(rsamp(i):rsamp(i)+14);
        end
        tth(c,:) = nanmean(tt_sample,1);
        
        % Find index of max
        [~,mx] = max(mean(tt_sample,1));
        
        % Determine whether data and shuffled distributions are different
        ranksum(tt_sample(:,mx),tt_shuff(:,mx));
        
        % Bonferroni corrected ranksum test
        [P,H] = ranksum(tt_sample(:,mx),tt_shuff(:,mx),'alpha',0.05/1552);
        
        tt_t(c) = H;
        tt_P(c) = P;
    end
    TTA{j} = tth;
    TT_tuned{j,1} = tt_t;
    TT_tuned{j,2} = tt_P;
end


%% Recompute touch tuning based on Benjamimi-Hochberg correction of p values
for j = 1:9;
    p_values = TT_tuned{j,2};
    
    [H,T] = benjaminihochberg(p_values,0.05); % Mark's code
    TT_tuned{j,3} = H;
end

% 95% confidence interval (Jeffreys Interval) on number of tuned cells
for j = 1:9; 
    N_tuned(j) = numel(find(TT_tuned{j,3}));
end

JIs(1,:) = betainv(0.025,N_tuned+0.5, 1552 - N_tuned +0.5);
JIs(2,:) = betainv(0.975,N_tuned+0.5, 1552 - N_tuned +0.5);

%% Bar graph of N_touch tuned
clf
bar(N_tuned,'k')
set(gca, 'XTickLabels', meth_names(gmeths),'XTickLabelRotation',45)
ylabel('N touch tuned')
hold all
errorbar(1:numel(N_tuned),N_tuned,N_tuned-ncells*JIs(1,:),ncells*JIs(2,:)-N_tuned,'color',[.5,.5,.5],'linewidth',2,'Linestyle','none')

%% Plot cell 495 as it is clearly touch tuned
cmap = distinguishable_colors(numel(methods));
clf;
for i = 1:9; 
    plot(zscore(TTA{i}(495,:)),'color',cmap(meths(i),:),'linewidth',2); hold all;
end
xlabel('Time (frames)')
ylabel('dF/F (z-scored)')
title('Touch triggered average (Cell 495)')
plot([7.5,7.5],[-2,4],'k--')
legend({meth_names{gmeths},'Touch'},'location','best')

