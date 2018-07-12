%% plot_noise_reject_all.m
% load Noise rejection table and plot pairs of variables

clear all
load('/Users/mathew/work/PopWorld/Results_batch1/Network_Rejection_Table_wStats.mat')
%% Identify learning subvolumes to colour separately
% NB: THERE IS SOMETHING WRONG WITH EVENTS FOR an197522 (animal 4)
a = 2:5;
SV  = 1;
Sess =  {[1:10]; [1:14]; [2:13,15:17]; [1:11,13]}; % {[2:10]; [2:14]; [2:13,15:17]; [2:11,13]}; %
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
plot([0,2500],[0,2500],'k--')

subplot(1,2,2)
var2 = var2./var1;
labels = {'N','Proportion of retained neurons'};
h(2) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)
legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','southeast')


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
legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')

%% Significantly contributing neurons vs dimensionality
var1 = Network_Rejection_Table.WCM_Dn;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_Dn;
var4 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(3); clf;
subplot(1,2,1)
labels = {'WCM_{Dn}';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'Config_{Dn}';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)
legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')

%% Significantly contributing neurons vs N
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.WCM_Dn;
var3 = Network_Rejection_Table.N;
var4 = Network_Rejection_Table.Config_Dn;


% plot_rejection_pairs(va1,va2,labels);
figure(4); clf;
subplot(1,2,1)
labels = {'N';'WCM_{Dn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'N';'Config_{Dn}'};
h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)

legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northwest')

%% Correlates with recording duration
var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(5); clf;
subplot(1,2,1)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'Time (samples)';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')

%% Same but split to learning sessions vs non learning sessions

var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;

l = find(Network_Rejection_Table.Learning == 1);
nl = 1:numel(var1);
nl(l) = [];

% plot_rejection_pairs(va1,va2,labels);
figure(6); clf;
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

%% Correlates with recording duration * N (i.e. number of datapoints)
var1 = Network_Rejection_Table.N .* Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(7); clf;
subplot(1,2,1)
labels = {'N * T';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'N * T';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

%% Same but split to learning sessions vs non learning sessions

var1 = Network_Rejection_Table.N .* Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Config_RejectionDn;

l = find(Network_Rejection_Table.Learning == 1);
nl = 1:numel(var1);
nl(l) = [];

% plot_rejection_pairs(va1,va2,labels);
figure(8); clf;
subplot(2,2,1)
labels = {'N * T';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1(l),var2(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(2,2,2)
labels = {'N * T';'Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var1(l),var3(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(2,2,3)
labels = {'N * T';'WCM_{RejectionDn}'};
h(3) = plot_rejection_pairs(var1(nl),var2(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

subplot(2,2,4)
labels = {'N * T';'Config_{RejectionDn}'};
h(4) = plot_rejection_pairs(var1(nl),var3(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

    
legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','best')

%% ksdensity (and scatter ?) of various variables
% TO Do plot median with label
var1 = Network_Rejection_Table.Config_Dn;
var2 = Network_Rejection_Table.Config_RejectionDn;
var3 = Network_Rejection_Table.Network_Size ./ Network_Rejection_Table.N;
var4 = Network_Rejection_Table.Config_RejectionDn ./Network_Rejection_Table.Network_Size;

labels = {'Config_{Dn}';'Config_{RejectionDn}';'Network Size/ N [P(retained)]';'Config_{RejectionDn}/Network Size [Norm dimensionality]'};

nbins = 50; 

figure(9); clf;
for c = 1:2; % Just the two methods
    these_d = find(colour_ID == c);
    for i = 1:4;
        bins = linspace(min(eval(['var',num2str(i)])),max(eval(['var',num2str(i)])),nbins);
        subplot(2,2,i);
        h = hist(eval(['var',num2str(i),'(these_d)']),bins);
%         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
        stairs(bins,h./sum(h),'color',dotcolour(these_d(1),:),'linewidth',2);
        xlabel(labels{i})
        hold all
    end
end

subplot(2,2,1)
legend('Peron','Calcium')

%% Same again but with the config model
%% ksdensity and scatter of various variables
% TO Do plot median with label
var1 = Network_Rejection_Table.WCM_Dn;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Network_Size ./ Network_Rejection_Table.N;
var4 = Network_Rejection_Table.WCM_RejectionDn ./Network_Rejection_Table.Network_Size;

labels = {'WCM_{Dn}';'WCM_{RejectionDn}';'Network Size/ N [P(retained)]';'WCM_{RejectionDn}/Network Size [Norm dimensionality]'};

nbins = 50; 

figure(10); clf;
for c = 1:2; % Just the two methods
    these_d = find(colour_ID == c);
    for i = 1:4;
        bins = linspace(min(eval(['var',num2str(i)])),max(eval(['var',num2str(i)])),nbins);
        subplot(2,2,i);
        h = hist(eval(['var',num2str(i),'(these_d)']),bins);
%         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
        stairs(bins,h./sum(h),'color',dotcolour(these_d(1),:),'linewidth',2);
        xlabel(labels{i})
        hold all
    end
end

subplot(2,2,1)
legend('Peron','Calcium')

%% 'Normalized dimensionality' by N or time
%% Correlates with recording duration * N (i.e. number of datapoints)

var1 = Network_Rejection_Table.WCM_RejectionDn ./ Network_Rejection_Table.WCM_Dn;
var2 = Network_Rejection_Table.Config_RejectionDn ./ Network_Rejection_Table.Config_Dn;
var3 = Network_Rejection_Table.N; 
var4 = Network_Rejection_Table.T;

% plot_rejection_pairs(va1,va2,labels);
figure(11); clf;
subplot(1,2,1)
labels = {'N','WCM_{RejectionDn} ./ WCM_{Dn}'};
h(1) = plot_rejection_pairs(var3,var1,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'N','Config_{RejectionDn} ./ Config_{Dn}'};
h(2) = plot_rejection_pairs(var3,var2,labels,colour_ID,edgecolour,dotcolour)

figure(12); clf;
subplot(1,2,1)
labels = {'T','WCM_{RejectionDn} ./ WCM_{Dn}'};
h(1) = plot_rejection_pairs(var4,var1,labels,colour_ID,edgecolour,dotcolour)

subplot(1,2,2)
labels = {'T','Config_{RejectionDn} ./ Config_{Dn}'};
h(1) = plot_rejection_pairs(var4,var2,labels,colour_ID,edgecolour,dotcolour)


%% Correlates with performance - Learning sessions only

L_sess = find(colour_ID>2);
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);
var6 = Network_Rejection_Table.Network_Size(L_sess) .* Network_Rejection_Table.T(L_sess);
var7 = Network_Rejection_Table.WCM_RejectionDn(L_sess) ./ Network_Rejection_Table.Network_Size(L_sess);
var8 = Network_Rejection_Table.eig90(L_sess);
% plot_rejection_pairs(va1,va2,labels);
figure(13); clf;
subplot(3,2,1)
labels = {'Session','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs_lines(var5,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,2,2)
labels = {'Session','Config_{RejectionDn}'};
h(2) = plot_rejection_pairs_lines(var5,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,2,3)
labels = {'Session','Pnolick'};
h(1) = plot_rejection_pairs_lines(var5,var3,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,2,4)
labels = {'Session','Pcorrect'};
h(1) = plot_rejection_pairs_lines(var5,var4,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,2,5)
labels = {'Session','Datapoints'};
h(1) = plot_rejection_pairs_lines(var5,var6,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,2,6)
% labels = {'Session','WCM_{RejectionDn}./Network Size'};
labels = {'Session','Eig90'};
h(1) = plot_rejection_pairs_lines(var5,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

%% Eigenvectors to reach 90% explained variance
% NL_sess = find(colour_ID<=2);
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.eig90;

% plot_rejection_pairs(va1,va2,labels);
figure(14); clf;
labels = {'N','Eig 90'};
plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

legend('Peron','Calcium')


%% Correlates with performance - 3 by 3

L_sess = find(colour_ID>2); % find(colour_ID>2 & colour_ID<7); % find(colour_ID > 6); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);
var6 = Network_Rejection_Table.Network_Size(L_sess) .* Network_Rejection_Table.T(L_sess);
var7 = Network_Rejection_Table.WCM_RejectionDn(L_sess) ./ Network_Rejection_Table.Network_Size(L_sess);
var8 = Network_Rejection_Table.eig90(L_sess);


% plot_rejection_pairs(va1,va2,labels);
figure(15); clf;

subplot(3,3,1)
labels = {'Pnolick','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var3,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,2)
labels = {'Pcorrect','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var4,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,3)
labels = {'Session','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs_lines(var5,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,4)
labels = {'Pnolick','Config_{RejectionDn}'};
h(2) = plot_rejection_pairs(var3,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,5)
labels = {'Pcorrect','Config_{RejectionDn}'};
h(1) = plot_rejection_pairs(var4,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,6)
labels = {'Session','Config_{RejectionDn}'};
h(1) = plot_rejection_pairs_lines(var5,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,7)
labels = {'Pnolick','Eig90'};
h(1) = plot_rejection_pairs(var3,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,8)
labels = {'Pcorrect','Eig90'};
h(1) = plot_rejection_pairs(var4,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(3,3,9)
% labels = {'Session','WCM_{RejectionDn}./Network Size'};
labels = {'Session','Eig90'};
h(1) = plot_rejection_pairs_lines(var5,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

%% Loop version of above

L_sess = find(colour_ID>2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

labelA = {'Pnolick';'Pcorrect';'Session'};
labelB = {'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'};
VARA = {var3;var4;var5}
VARB = {var1;var2;var8};

figure(16); clf;
hold all
n = 0;
for i = 1:3;
    for j = 1:3;
        n = n + 1;
        subplot(3,3,n);
        h(1) = plot_rejection_pairs(VARA{j},VARB{i},{labelA{j},labelB{i}},colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))
    end
end
% suptitle('Peron')

%% Need a new variable that counts the sessions that a given subvolume was 
% present in. I.e. if the first time a cell was imaged on was day 8, it's
% 'imaging day' 8 would be marked as 1.
% Basing this on N_unique (defined in peron_noise_reject_all.m)

% First, we need to process the methods separately
meths = {'calcium';'Peron'}
for m = 1:2
    this_m = find(strcmp(Network_Rejection_Table.method,meths{m}));
    all_svs = unique(Network_Rejection_Table.N_unique);
    for i = 1:numel(all_svs)
        this_sv = find(Network_Rejection_Table.N_unique(this_m)==all_svs(i));
        sv_dates = [];
        for s = 1:numel(this_sv)
            Y = Network_Rejection_Table.NetworkName{this_m(this_sv(s))}(10:13);
            M = Network_Rejection_Table.NetworkName{this_m(this_sv(s))}(15:16);
            D = Network_Rejection_Table.NetworkName{this_m(this_sv(s))}(18:19);
            sv_dates = [sv_dates, datenum([Y,'-',M,'-',D])];
        end
        Network_Rejection_Table.datenum(this_m(this_sv)) = sv_dates;
        % Sort dates by datenum and assign increasing values
        [~,sv_order] = sort(sv_dates);
        Network_Rejection_Table.sv_order(this_m(this_sv)) = sv_order;
    end
end

%% Plot N by datenum - quick way of looking at the experiment
clf
for a = 1:8
    this_a = find(Network_Rejection_Table.Animal(this_m)==a);
    plot(Network_Rejection_Table.datenum(this_m(this_a)),Network_Rejection_Table.N(this_m(this_a)),'.');
    hold all;
end
        
        
        
this_a = find(Network_Rejection_Table.Animal == 2);

%% Plot non-learning sessions vs session order
L_sess = find(colour_ID==2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.sv_order(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

subvols = Network_Rejection_Table.N_unique(L_sess);

labelA = {'Pnolick';'Pcorrect';'SV_{order}'};
labelB = {'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'};
VARA = {var3;var4;var5};
VARB = {var1;var2;var8};

figure(20); clf;
hold all
n = 0;
for i = 1:3
    for j = 1:3
        n = n + 1;
        subplot(3,3,n);
        h(1) = plot_rejection_pairs_lines(VARA{j},VARB{i},{labelA{j},labelB{i}},subvols,edgecolour(L_sess,:),dotcolour(L_sess,:));        
        
    end
end
suptitle('Calcium')

%% Add distributions of each 1D variable

L_sess = find(colour_ID<=2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
% var5 = Network_Rejection_Table.Session(L_sess);
var5 = Network_Rejection_Table.sv_order(L_sess);
subvols = Network_Rejection_Table.N_unique(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

labelA = {'Pnolick';'Pcorrect';'Session'};
labelB = {'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'};
VARA = {var3;var4;var5};
VARB = {var1;var2;var8};

figure(17); clf;
hold all
nbins = 50; 

for i = 1:3
    bins = linspace(min(VARA{i}),max(VARA{i}),nbins);
    subplot(4,4,i);
    h = hist(VARA{i},bins);
    %         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
    stairs(bins,h./sum(h),'k','linewidth',2);
    xlabel(labelA{i});
    axis square
    hold all
end


for i = 1:3;
    for j = 1:3;
        n = 4*i + j;
        subplot(4,4,n);
        h(1) = plot_rejection_pairs(VARA{j},VARB{i},{labelA{j},labelB{i}},colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:));
        axis square
    end
end

for c = 1:2;
    these_d = find(colour_ID(L_sess)==c);
    for i = 1:3
        bins = linspace(min(VARB{i}(these_d)),max(VARB{i}(these_d)),nbins);
        subplot(4,4,4+i*4);
        h = hist(VARB{i}(these_d),bins);
        %         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
        stairs(h./sum(h),bins,'color',dotcolour(L_sess(these_d(1)),:),'linewidth',2);
        ylabel(labelB{i});
        axis square
        hold all
    end
end
% suptitle('Peron')

%% Add distributions of each 1D variable

L_sess = find(colour_ID==1); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess); 
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

labelA = {'Pnolick';'Pcorrect';'Session'};
labelB = {'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'};
VARA = {var3;var4;var5};
VARB = {var1;var2;var8};

figure(18); clf;
hold all
nbins = 40; 

for i = 1:3
    bins = linspace(min(VARA{i}),max(VARA{i}),nbins);
    subplot(4,4,i);
    h = hist(VARA{i},bins);
    %         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
    stairs(bins,h./sum(h),'k','linewidth',2);
    xlabel(labelA{i});
    axis square
    hold all
end


for i = 1:3
    for j = 1:3
        n = 4*i + j;
        subplot(4,4,n);
        [X,imax] = plot_rejection_image(VARA{j},VARB{i},{labelA{j},labelB{i}},nbins/2);
        axis square

%         if i == 1
%             subplot(4,4,j)
%             stairs(mean(X)); xlim([1,nbins+1]);
%             xticks([5,10,15,20]); xtickformat('%2.2f');
%             xticklabels(c{1}([5,10,15,20]));
%             
%         end
        
%         if j == 1
%             subplot(4,4,4*i+4)
%             stairs(mean(X'),linspace(nbins,1,nbins)); ylim([1,nbins+1])
%         end
    end
end

for c = 1:2
    these_d = find(colour_ID(L_sess)==c);
    if these_d
        for i = 1:3
            bins = linspace(min(VARB{i}(these_d)),max(VARB{i}(these_d)),nbins);
            subplot(4,4,4+i*4);
            h = hist(VARB{i}(these_d),bins);
            %         bar(bins,h./sum(h),'facecolor','none','edgecolor',dotcolour(these_d(1),:),'linewidth',2);
            stairs(h./sum(h),bins,'color',dotcolour(L_sess(these_d(1)),:),'linewidth',2);
            ylabel(labelB{i});
            axis square
            hold all
        end
    end
end


suptitle('Peron')

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

%% Stripped down single plot figure with lines
function h = plot_rejection_pairs_lines(var1,var2,labels,colour_ID,edgecolour,dotcolour)
cla; hold all
C = unique(colour_ID);
for c = 1:numel(C)
    these_d = find(colour_ID == C(c));
    
    h = plot(var1(these_d),var2(these_d),'color',dotcolour(these_d(1),:),'linewidth',2);
    h = plot(var1(these_d),var2(these_d),'o','markeredgecolor',edgecolour(these_d(1),:),'markerfacecolor',dotcolour(these_d(1),:),'markersize',5);
end
axis square
xlabel(labels{1});
ylabel(labels{2});
end

%% Shaded error bar version of above
function [X,c] = plot_rejection_image(var1,var2,labels,nbins)
cla; hold all

binsA = linspace(min(var1),max(var1),nbins);
binsB = linspace(min(var2),max(var2),nbins);

[X,c] = hist3([var1,var2],'Ctrs',{binsA,binsB}); 
X = rot90(X);
X = flipud(X)
h = imagesc(c{1},c{2},X);
% xlim([floor(min(c{1})),ceil(max(c{1}))])
ylim([floor(min(c{2})),ceil(max(c{2}))])
% set(gca,'YDir','reverse')
ticklabs = get(gca,'xticklabels')
axis square
xlabel(labels{1});
ylabel(labels{2});
end