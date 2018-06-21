%% plot_noise_reject_all.m
% load Noise rejection table and plot pairs of variables

clear all
load('/Users/mathew/work/PopWorld/Results_batch1/Network_Rejection_Table_wStats.mat')
%% Identify learning subvolumes to colour separately
% NB: THERE IS SOMETHING WRONG WITH EVENTS FOR an197522 (animal 4)
a = 2:5;
SV  = 1;
Sess = {[2:10]; [2:14]; [2:13,15:17]; [2:11,13]}; % {[1:10]; [1:14]; [2:13,15:17]; [1:11,13]};
L_sess = [];
L_sess_index = [];
Network_Rejection_Table.Learning = zeros(height(Network_Rejection_Table),1);

for i = 1:4;
    sessions = Sess{i};
    this_m = find(Network_Rejection_Table.Animal == a(i));
    these_sv = find(Network_Rejection_Table.Subvolume(this_m) == 1);
    these_sessions = find(ismember(Network_Rejection_Table.Session(this_m(these_sv)),sessions));
    L_sess_index = this_m(these_sv(these_sessions));
    L_sess = [L_sess; Network_Rejection_Table(L_sess_index,:)];
    % Add learning session variable to Network_Rejection_Table
    Network_Rejection_Table.Learning(L_sess_index) = 1;
end

%% Adding two variables (3 columns each) with plotting colours (for dots and edges):
% Dots: black/grey depending on method
% Edges - Not learning: black or grey based on method.
%         Learning: 4 varycolor colours

colours = varycolor(4);

% Dot colour based on method
dotcolour = zeros(height(Network_Rejection_Table),3);
edgecolour = zeros(height(Network_Rejection_Table),3);
colour_ID = ones(height(Network_Rejection_Table),1);
 
for n = 1:height(Network_Rejection_Table)
    if strcmp(Network_Rejection_Table.method{n},'calcium')
        dotcolour(n,:) = [.5,.5,.5];
        edgecolour(n,:) = [.5,.5,.5];
        colour_ID(n) = 2;
        
        if Network_Rejection_Table.Learning(n) == 1;
            edgecolour(n,:) = colours(Network_Rejection_Table.Animal(n) - 1,:);
            colour_ID(n) = Network_Rejection_Table.Animal(n)+1;
        end
    else
        if Network_Rejection_Table.Learning(n) == 1;
            edgecolour(n,:) = colours(Network_Rejection_Table.Animal(n) - 1,:);
            colour_ID(n) = Network_Rejection_Table.Animal(n)+5;
        end
        
    end
end
%% First, actual network size vs signal network size
% (var1,var2,labels,colours)
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.Network_Size;
labels = {'N';'Signal network size'};

% plot_rejection_pairs(va1,va2,labels);
figure(1); clf;
subplot(1,2,1)
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
var2 = var2./var1;
labels = {'N','Proportion of retained neurons'};
h(2) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

%% network size vs dimensions
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;
labels = {'N';'WCM_{RejectionDn}'};

% plot_rejection_pairs(va1,va2,labels);
figure(2); clf;
subplot(1,2,1)
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'N';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

%% Significantly contributing neurons vs dimensionality
var1 = Network_Rejection_Table.WCM_Dn;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_Dn;
var4 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(6); clf;
subplot(1,2,1)
labels = {'WCM_{Dn}';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'Config_{Dn}';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)

%% Correlates with recording duration
var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(7); clf;
subplot(1,2,1)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'Time (samples)';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

%% Same but split to learning sessions vs non learning sessions

var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;

l = find(Network_Rejection_Table.Learning == 1);
nl = 1:numel(var1);
nl(l) = [];

% plot_rejection_pairs(va1,va2,labels);
figure(8); clf;
subplot(2,2,1)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1(l),var2(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(2,2,2)
labels = {'Time (samples)';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1(l),var3(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(2,2,3)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(3) = plot_rejection_pairs(var1(nl),var2(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

subplot(2,2,4)
labels = {'Time (samples)';'Config_{RejectionDn}'};
h(4) = plot_rejection_pairs(var1(nl),var3(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

%%

for c = 1:numel(unique(colour_ID));
    these_d = find(colour_ID == c);
    plot(var1(these_d),var2(these_d),'o','markeredgecolor',edgecolour(these_d(1),:),'markerfacecolor',dotcolour(these_d(1),:),'markersize',5)
end
axis square
xlabel(labels{1});
ylabel(labels{2});

% Normalizing one by the other - should be flat
subplot(1,2,2)
hold all
for c = 1:numel(unique(colour_ID));
    these_d = find(colour_ID == c);
    plot(var1(these_d),var2(these_d)./var1(these_d),'o','markeredgecolor',edgecolour(these_d(1),:),'markerfacecolor',dotcolour(these_d(1),:),'markersize',5)
end
axis square
xlabel(labels{1});
ylabel('Proportion of retained neurons');

legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','best')


%% Network size vs WCM_Dn/Config_Dn
figure(1); clf;
methods = {'Peron';'calcium'};
dotcolours = [0,0,0;0.5,0.5,0.5];%varycolor(2);
subplot(1,3,3)
plot(0,0,'o','markeredgecolor',dotcolours(1,:),'markerfacecolor',dotcolours(1,:),'markersize',5); hold all
plot(0,0,'o','markeredgecolor',dotcolours(2,:),'markerfacecolor',dotcolours(2,:),'markersize',5)
for m = 1:2
    clear these_m
    for n = 1:height(Network_Rejection_Table)
        these_m(n) = strcmp(Network_Rejection_Table.method{n},methods{m});
    end
    
    m_array = find(these_m);
    subplot(1,3,1)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.WCM_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.Network_Size,Network_Rejection_Table.WCM_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,2000],[0,2000],'k--')
    xlim([0,2000])
    ylim([0,1000])
    xlabel('Network Size')
    ylabel(['WCM_{Dn}'])
    axis square
    
    subplot(1,3,2)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.Config_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.Network_Size,Network_Rejection_Table.Config_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,2000],[0,2000],'k--')
    xlim([0,2000])
    ylim([0,800])
    xlabel('Network Size')
    ylabel(['Config_{Dn}'])
    axis square
    
    subplot(1,3,3)
    plot(Network_Rejection_Table.WCM_Dn(m_array),Network_Rejection_Table.Config_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.WCM_Dn,Network_Rejection_Table.Config_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,1000],[0,1000],'k--')
    xlim([0,1000])
    ylim([0,800])
    xlabel(['WCM_{Dn}'])
    ylabel(['Config_{Dn}'])
    axis square
end


%% Plot learning subvolumes and in separate colours (coloured rings around dots)
figure(1);
methods = {'Peron';'calcium'};
colours = varycolor(10);
edgecolours = colours([2,4,6,8],:);
subplot(1,3,3)
for i = 1:4;
    plot(-5,-5,'o','markeredgecolor',edgecolours(i,:),'markerfacecolor',dotcolours(1,:),'markersize',5);
end
for m = 1:2;
    for animal = 2:5;
        this_a = find(L_sess.Animal == animal);
        clear these_m
        for n = 1:numel(this_a)
            these_m(n) = strcmp(L_sess.method{this_a(n)},methods{m});
        end
        
        m_array = this_a(find(these_m));
        subplot(1,3,1)
        plot(L_sess.Network_Size(m_array),L_sess.WCM_Dn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        
        subplot(1,3,2)
        plot(L_sess.Network_Size(m_array),L_sess.Config_Dn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        
        subplot(1,3,3)
        plot(L_sess.WCM_Dn(m_array),L_sess.Config_Dn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        
        
    end
end


legend('Peron','Calcium','location','best');

%% Network size vs WCM_RejectionDn/Config_RejectionDn
figure(2); clf
subplot(1,3,3)
dotcolours = [0,0,0;0.5,0.5,0.5];%varycolor(2);
subplot(1,3,3)
% plot(0,0,'o','markeredgecolor',dotcolours(1,:),'markerfacecolor',dotcolours(1,:),'markersize',5); hold all
% plot(0,0,'o','markeredgecolor',dotcolours(2,:),'markerfacecolor',dotcolours(2,:),'markersize',5)

for m = 1:2;
    for n = 1:height(Network_Rejection_Table)
        these_m(n) = strcmp(Network_Rejection_Table.method{n},methods{m});
    end
    
    m_array = find(these_m);
    subplot(1,3,1)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.WCM_RejectionDn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    hold all
    % plot([0,2000],[0,2000],'k--')
    % xlim([0,2000])
    % ylim([0,1000])
    xlabel('Network Size')
    ylabel(['WCM_{RejectionDn}'])
    axis square
    
    subplot(1,3,2)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.Config_RejectionDn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    hold all
    % plot([0,2000],[0,2000],'k--')
    % xlim([0,2000])
    % ylim([0,800])
    xlabel('Network Size')
    ylabel(['Config_{RejectionDn}'])
    axis square
    
    subplot(1,3,3)
    plot(Network_Rejection_Table.WCM_RejectionDn(m_array),Network_Rejection_Table.Config_RejectionDn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    hold all
    % plot([0,1000],[0,1000],'k--')
    % xlim([0,1000])
    % ylim([0,800])
    xlabel(['WCM_{RejectionDn}'])
    ylabel(['Config_{RejectionDn}'])
    axis square
end

% legend('Peron','Calcium','location','northeast');
%% Plot learning subvolumes and in separate colours (coloured rings around dots)
figure(2)
colours = varycolor(10);
edgecolours = colours([2,4,6,8],:);
subplot(1,3,3)
for i = 1:4;
    plot(-5,-5,'o','markeredgecolor',edgecolours(i,:),'markerfacecolor',dotcolours(1,:),'markersize',5);
end
for m = 1:2;
    for animal = 2:5;
        this_a = find(L_sess.Animal == animal);
        clear these_m
        for n = 1:numel(this_a)
            these_m(n) = strcmp(L_sess.method{this_a(n)},methods{m});
        end
        
        m_array = this_a(find(these_m));
        subplot(1,3,1);
        plot(L_sess.Network_Size(m_array),L_sess.WCM_RejectionDn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        hold all
        xlabel('Network Size')
        ylabel(['WCM_{RejectionDn}'])
        axis square
        
        subplot(1,3,2)
        plot(L_sess.Network_Size(m_array),L_sess.Config_RejectionDn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        hold all
        xlabel('Network Size')
        ylabel(['Config_{RejectionDn}'])
        axis square
        
        subplot(1,3,3)
        plot(L_sess.WCM_RejectionDn(m_array),L_sess.Config_RejectionDn(m_array),'o','markeredgecolor',edgecolours(animal-1,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
        hold all
        xlabel(['WCM_{RejectionDn}'])
        ylabel(['Config_{RejectionDn}'])
        axis square
    end
end

legend('Peron','Calcium','Learning SV 1','Learning SV 2','Learning SV 3','Learning SV 4','location','best')

%% ksdensity and scatter of various variables
% TO Do plot median with label
figure(3); clf;
subplot(2,2,1);
% [y,x] = ksdensity(Network_Rejection_Table.WCM_Dn);
% plot(x,y,'k','linewidth',2)
[y,x] = hist(Network_Rejection_Table.WCM_Dn,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(450,43,['Peak = ',num2str(x(mx),'%.3g')])
text(450,35,['Median = ',num2str(median(Network_Rejection_Table.WCM_Dn),'%.3g')])
xlabel(['WCM_{Dn}'])
ylabel('Probability density')
% title('Nretain')

subplot(2,2,2);
% [y,x] = ksdensity(Network_Rejection_Table.WCM_RejectionDn);
[y,x] = hist(Network_Rejection_Table.WCM_RejectionDn,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(60,50,['Peak = ',num2str(x(mx),'%.3g')])
text(60,43,['Median = ',num2str(median(Network_Rejection_Table.WCM_RejectionDn),'%.3g')])
xlabel(['WCM_{RejectionDn}'])
ylabel('Probability density')
% title('Npos')

subplot(2,2,3);
data = Network_Rejection_Table.Network_Size./Network_Rejection_Table.WCM_RejectionDn;
% data(find(isinf(data))) = [];
% [y,x] = ksdensity(data);
% plot(x,y,'k','linewidth',2)
[y,x] = hist(data,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(60,50,['Peak = ',num2str(x(mx),'%.3g')])
text(60,43,['Median = ',num2str(median(data),'%.3g')])
xlabel(['Network Size / WCM_{RejectionDn}'])
ylabel('Probability density')
% title('Pretain')

subplot(2,2,4);
data = Network_Rejection_Table.WCM_Dn./Network_Rejection_Table.WCM_RejectionDn;
[y,x] = hist(data,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(25,67,['Peak = ',num2str(x(mx),'%.3g')])
text(25,58,['Median = ',num2str(median(data),'%.3g')])
xlabel(['WCM_{Dn} / WCM_{RejectionDn}'])
ylabel('Probability density')
% title('Normalized dimensionality')

suptitle('WCM results overview')
% print(gcf,'-depsc','-painters','noiserejection/figs/NRksdensity.eps')

%% Repeat with Config version
% TO Do plot median with label
figure(4); clf;
subplot(2,2,1);
% [y,x] = ksdensity(Network_Rejection_Table.WCM_Dn);
% plot(x,y,'k','linewidth',2)
[y,x] = hist(Network_Rejection_Table.Config_Dn,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(450,35,['Peak = ',num2str(x(mx),'%.3g')])
text(450,30,['Median = ',num2str(median(Network_Rejection_Table.Config_Dn),'%.3g')])
xlabel(['Config_{Dn}'])
ylabel('Probability density')
% title('Nretain')

subplot(2,2,2);
% [y,x] = ksdensity(Network_Rejection_Table.WCM_RejectionDn);
[y,x] = hist(Network_Rejection_Table.Config_RejectionDn,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(25,165,['Peak = ',num2str(x(mx),'%.3g')])
text(25,140,['Median = ',num2str(median(Network_Rejection_Table.Config_RejectionDn),'%.3g')])
xlabel(['Config_{RejectionDn}'])
ylabel('Probability density')
% title('Npos')

subplot(2,2,3);
data = Network_Rejection_Table.Network_Size./Network_Rejection_Table.Config_RejectionDn;
% data(find(isinf(data))) = [];
% [y,x] = ksdensity(data);
% plot(x,y,'k','linewidth',2)
[y,x] = hist(data,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(800,160,['Peak = ',num2str(x(mx),'%.3g')])
text(800,140,['Median = ',num2str(median(data),'%.3g')])
xlabel(['Network Size / Config_{RejectionDn}'])
ylabel('Probability density')
% title('Pretain')

subplot(2,2,4);
data = Network_Rejection_Table.WCM_Dn./Network_Rejection_Table.Config_RejectionDn;
[y,x] = hist(data,100);
bar(x,y,'k','linewidth',2)
hold all; xlim([min(x),max(x)]); [~,mx] = max(y);
text(220,115,['Peak = ',num2str(x(mx),'%.3g')])
text(220,100,['Median = ',num2str(median(data),'%.3g')])
xlabel(['Config_{Dn} / Config_{RejectionDn}'])
ylabel('Probability density')
% title('Normalized dimensionality')

suptitle('Default config model results overview')
% print(gcf,'-depsc','-painters','noiserejection/figs/NRksdensity.eps')

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


%% Local plotting functions (required Matlab 2016a or later)

figure;
methods = labels;
dotcolours = colours;
subplot(1,3,3);
plot(0,0,'o','markeredgecolor',dotcolours(1,:),'markerfacecolor',dotcolours(1,:),'markersize',5); hold all
plot(0,0,'o','markeredgecolor',dotcolours(2,:),'markerfacecolor',dotcolours(2,:),'markersize',5)
for m = 1:2
    clear these_m
    for n = 1:height(Network_Rejection_Table)
        these_m(n) = strcmp(Network_Rejection_Table.method{n},methods{m});
    end
    
    m_array = find(these_m);
    subplot(1,3,1)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.WCM_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.Network_Size,Network_Rejection_Table.WCM_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,2000],[0,2000],'k--')
    xlim([0,2000])
    ylim([0,1000])
    xlabel('Network Size')
    ylabel(['WCM_{Dn}'])
    axis square
    
    subplot(1,3,2)
    plot(Network_Rejection_Table.Network_Size(m_array),Network_Rejection_Table.Config_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.Network_Size,Network_Rejection_Table.Config_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,2000],[0,2000],'k--')
    xlim([0,2000])
    ylim([0,800])
    xlabel('Network Size')
    ylabel(['Config_{Dn}'])
    axis square
    
    subplot(1,3,3)
    plot(Network_Rejection_Table.WCM_Dn(m_array),Network_Rejection_Table.Config_Dn(m_array),'o','markeredgecolor',dotcolours(m,:),'markerfacecolor',dotcolours(m,:),'markersize',5)
    % plot(Network_Rejection_Table.WCM_Dn,Network_Rejection_Table.Config_Dn,'.','color',[.5,.5,.5])
    hold all
    plot([0,1000],[0,1000],'k--')
    xlim([0,1000])
    ylim([0,800])
    xlabel(['WCM_{Dn}'])
    ylabel(['Config_{Dn}'])
    axis square
end
% end

%% Stripped down single plot figure
function h = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)
cla; hold all
C = unique(colour_ID);
for c = 1:numel(C)
    these_d = find(colour_ID == C(c));
    h = plot(var1(these_d),var2(these_d),'o','markeredgecolor',edgecolour(these_d(1),:),'markerfacecolor',dotcolour(these_d(1),:),'markersize',5)
end
axis square
xlabel(labels{1});
ylabel(labels{2});
end