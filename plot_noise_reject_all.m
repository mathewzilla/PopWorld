%% plot_noise_reject_all.m
% load Noise rejection table and plot pairs of variables

clear all
% load('/Users/mathew/work/PopWorld/Results_batch1/Network_Rejection_Table_wStats.mat')
load('Results_reject_preround/Network_Rejection_Table_wStats_preround3.mat') % Includes learning sessions
%% Identify learning subvolumes to colour separately
% NB: THERE IS SOMETHING WRONG WITH EVENTS FOR an197522 (animal 4)
a = 2:5;
SV  = 1;
Sess =  {[1:10]; [1:14]; [2:13,15:17]; [1:11,13]}; % {[2:10]; [2:14]; [2:13,15:17]; [2:11,13]}; %
L_sess = [];
L_sess_index = [];
Network_Rejection_Table.Learning = zeros(height(Network_Rejection_Table),1);

for i = 1:4
    sessions = Sess{i};
    this_m = find(Network_Rejection_Table.Animal == a(i));
    these_sv = find(Network_Rejection_Table.Subvolume(this_m) == 1);
    these_sessions = find(ismember(Network_Rejection_Table.Session(this_m(these_sv)),sessions));
    L_sess_index = this_m(these_sv(these_sessions));
    L_sess = [L_sess; Network_Rejection_Table(L_sess_index,:)];
    % Add learning session variable to Network_Rejection_Table
    Network_Rejection_Table.Learning(L_sess_index) = 1;
end

%% Load events version and integrate with Network Rejection Table
% temp = load('/Users/mathew/work/PopWorld/Results_reject_events/Network_Rejection_Table_events.mat');
temp = load('Results_reject_events/Network_Rejection_Table_events.mat');
events_table = temp.Network_Rejection_Table; clear temp
for i =1:height(events_table)
    Network_Rejection_Table.Event_Network_Size = events_table.Network_Size;
    Network_Rejection_Table.Event_Signal_Size_WCM = events_table.Signal_Size_WCM;
    Network_Rejection_Table.Event_WCM_Dn = events_table.WCM_Dn;
    Network_Rejection_Table.Event_WCM_RejectionDn = events_table.WCM_RejectionDn;
    Network_Rejection_Table.Event_Signal_Components = events_table.Signal_Components;
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

%% Simpler animal based colours, with not colour indicating learning session
% Dots: black/grey based on learning/not learning
% Edges - 8 animals

colours = varycolor(8);

% Dot colour based on method
dotcolour = zeros(height(Network_Rejection_Table),3);
edgecolour = zeros(height(Network_Rejection_Table),3);
colour_ID = ones(height(Network_Rejection_Table),1);

for n = 1:height(Network_Rejection_Table)
    %     if Network_Rejection_Table.Learning(n) == 1
    dotcolour(n,:) = [.5,.5,.5];
    edgecolour(n,:) = colours(Network_Rejection_Table.Animal(n),:);
    colour_ID(n) = Network_Rejection_Table.Animal(n);
    %     else
    %         dotcolour(n,:) = [0,0,0];
    %         edgecolour(n,:) = colours(Network_Rejection_Table.Animal(n),:);
    %         colour_ID(n) = Network_Rejection_Table.Animal(n);
    %     end
end

%% First, actual network size vs signal network size
% (var1,var2,labels,colours)
var1 = Network_Rejection_Table.N;
% var2 = Network_Rejection_Table.Network_Size;
var2 = Network_Rejection_Table.Event_Network_Size;
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
% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','southeast')
legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5','Mouse 6','Mouse 7','Mouse 8','location','southeast')


%% N/network size vs dimensions
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.WCM_RejectionDn;
var3 = Network_Rejection_Table.Network_Size;

var4 = Network_Rejection_Table.eig90;
var4 = Network_Rejection_Table.Event_WCM_RejectionDn;
% labels = {'N';'WCM_{RejectionDn}'};

% plot_rejection_pairs(va1,va2,labels);
figure(2); clf;
subplot(1,2,1)
% h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

% subplot(1,2,2)
% labels = {'Retained N';'Rejection Dn'};
% labels = {'Retained N';'Eig 90'};
labels = {'Retained N';'Events Rejection Dn'};
h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)
% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')
% legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5','Mouse 6','Mouse 7','Mouse 8','location','southeast')

%% Plot eig90 vs D_retained for all animals, then each individually
NL_sess = find(~Network_Rejection_Table.Learning);
% Also strip out very long sessions
short_sess = find(Network_Rejection_Table.T<=5000);
keepers = NL_sess(ismember(NL_sess,short_sess));

var2 = Network_Rejection_Table.WCM_RejectionDn;
var4 = Network_Rejection_Table.eig90;

% Normalise vs zscore
var2 = var2./max(var2(keepers));
var4 = var4./max(var4(keepers));

% var2 = local_z(var2,keepers);
% var4 = local_z(var4,keepers);


figure(2);clf;
subplot(1,2,1)
hold all
labels = {'Calcium Spectral D (%)';'Eig_{90} dimensions (%)'};
h(2) = plot_rejection_pairs(var2(keepers)*100,var4(keepers)*100,labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))
plot([0,100],[0,100],'linewidth',2,'color',[.7,.7,.7]); grid on


var2 = Network_Rejection_Table.Event_WCM_RejectionDn;
var2 = var2./max(var2(keepers));
subplot(1,2,2)
hold all
labels = {'Events Spectral D (%)';'Eig_{90} dimensions (%)'};
h(2) = plot_rejection_pairs(var2(keepers)*100,var4(keepers)*100,labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))
plot([0,100],[0,100],'linewidth',2,'color',[.7,.7,.7]); grid on
%% title('Proportion')

figure(101); clf;
ax(1) = subplot(3,3,1);

h(2) = plot_rejection_pairs(var2(keepers),var4(keepers),labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))
% grid on
animals = Network_Rejection_Table.Animal;

for i = 1:numel(unique(colour_ID))
    ax(i+1) = subplot(3,3,i+1);
    these_As = find(animals(keepers) == i);
    this_A = keepers(these_As);
    h(2) = plot_rejection_pairs(var2(this_A),var4(this_A),labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
%     grid on
end
suptitle ('Proportion')

%% Same but not plotting learning sessions
NL_sess = find(~Network_Rejection_Table.Learning);

figure(102); clf;
ax(1) = subplot(3,3,1);
labels = {'D_{rejection}';'Eig_{90}'};
h(1) = plot_rejection_pairs(var2(NL_sess),var4(NL_sess),labels,colour_ID(NL_sess),edgecolour(NL_sess,:),dotcolour(NL_sess,:))

animals = Network_Rejection_Table.Animal;


for i = 1:numel(unique(colour_ID))
    ax(i+1) = subplot(3,3,i+1);
    these_As = find(animals(NL_sess) == i);
    this_A = NL_sess(these_As);
    h(i+1) = plot_rejection_pairs(var2(this_A),var4(this_A),labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
end
linkaxes(ax)

%% Fit linear regression line and extrapolate to ~10000 barrel column population (this number might contain other layers)B = x\var2;
[x,ix] = sort(var3);
y = var2(ix);
B = x\y;
yEst1 = x*B;
plot(x,yEst1,'r','linewidth',2)

% Add intercept
X = [ones(length(x),1) x];
B2 = X\y;
yEst2 = X*B2;
plot(x,yEst2,'b','linewidth',2)

% Use polyfit for higher order regression
p = polyfit(x,y,2);
yFit = polyval(p,x);
plot(x,yFit,'g','linewidth',2)

p2 = polyfit(x,y,3);
yFit2 = polyval(p2,x);
plot(x,yFit2,'c','linewidth',2)


%% Extrapolating out to 10K cells
xrange = linspace(1,1e4,100)';
yExt1 = xrange*B;
plot(xrange,yExt1,'r','linewidth',2)

yExt2 = [ones(100,1), xrange]*B2;
plot(xrange,yExt2,'b','linewidth',2)

yExt3 = polyval(p,xrange);
plot(xrange,yExt3,'g','linewidth',2)

yExt4 = polyval(p2,xrange);
plot(xrange,yExt4,'c','linewidth',2)

%% GP version
modelgp = fitrgp(x,y);
[ypred,ysd,yint] = predict(modelgp,xrange);
% myeb(xrange',ypred,ysd);
plot(xrange,ypred,'y','linewidth',2)

%%
xlim([0,2500])
% ylim([0,60])
ylim([0,1400])

%%
% ylim([0,200])
ylim([0,7000])
xlim([0,10000])
% plot([0;1e4],yExt1(end).*ones(2,1),'r--')
% plot([0;1e4],yExt2(end).*ones(2,1),'b--')
plot([0;1e4],yExt1(end).*ones(2,1),'r--')
plot([0;1e4],yExt3(end).*ones(2,1),'g--')

%% Non-learning data. Shorter sessions only
NL_sess = find(~Network_Rejection_Table.Learning);
% Also strip out very long sessions
short_sess = find(Network_Rejection_Table.T<=5000);

keepers = NL_sess(ismember(NL_sess,short_sess));

% var1 = Network_Rejection_Table.Signal_Size_WCM;
var1 = Network_Rejection_Table.Event_Signal_Size_WCM;
% var2 = Network_Rejection_Table.WCM_RejectionDn;
var2 = Network_Rejection_Table.Event_WCM_RejectionDn;
var3 = Network_Rejection_Table.Network_Size;
var4 = Network_Rejection_Table.eig90;

figure(2); clf;
labels = {'Retained N';'D_{Rejection}'};
h(1) = plot_rejection_pairs(var1(keepers),var2(keepers),labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))


figure(103); clf;
ax(1) = subplot(3,3,1);
h(1) = plot_rejection_pairs(var1(keepers),var2(keepers),labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))

animals = Network_Rejection_Table.Animal;


for i = 1:numel(unique(colour_ID))
    ax(i+1) = subplot(3,3,i+1);
    these_As = find(animals(keepers) == i);
    this_A = keepers(these_As);
    h(i+1) = plot_rejection_pairs(var1(this_A),var2(this_A),labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
end
linkaxes(ax)

%% Fitting linear and rise to max models (Mark's code). 
clear AICs BICs
% var1 = Network_Rejection_Table.Signal_Size_WCM;
% var2 = Network_Rejection_Table.WCM_RejectionDn;
% var2 = Network_Rejection_Table.eig90;
var1 = Network_Rejection_Table.Event_Signal_Size_WCM;
var2 = Network_Rejection_Table.Event_WCM_RejectionDn;
xdata = var1(keepers);
ydata = var2(keepers);

SS = zeros(2,1);
coeffs = cell(2,1);
residuals = cell(2,1);
IV = rand(2,1);

% linear fit
LB = -1e6;
UB = 1e6;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
options.MaxFunctionEvaluations = 1000;
    
linfun = inline('x(1) + x(2) .* xdata','x','xdata');
[coeffs{1},SS(1),residuals{1}] = lsqcurvefit(linfun,[1;0], xdata, ydata,LB,UB,options);   

% exponential rise-to-max (2par) a*(1-exp(-b*x))
expmaxAfun = inline('x(1) .* (1-exp(-x(2).*xdata))','x','xdata');
[coeffs{2},SS(2),residuals{2}] = lsqcurvefit(expmaxAfun,[100;0], xdata, ydata,LB,UB,options);              

% Model selection stuff
for j = 1:2
    AICs(j) = AICSS(SS(j),length(ydata),length(coeffs{j}));
    BICs(j) = BICSS(SS(j),length(ydata),length(coeffs{j}));
end

%% Linear through the origin
B = xdata\ydata;
yEst1 = xdata*B;
% plot(xdata,yEst1,'g','linewidth',2)

%% Plotting
% figure(104); clf;
% subplot(1,2,1);
figure(2);
% labels = {'Retained N';'D_{Rejection}'};
labels = {'N';'Spectral rejection dimensions'};
% labels = {'N';'Eig_{90} dimensions'};
h(1) = plot_rejection_pairs(xdata,ydata,labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:));

xrange = linspace(0,6e3,100);

linpred = linfun(coeffs{1},xrange);
yEst2 = xrange*B;

r2maxpred = expmaxAfun(coeffs{2},xrange);

% plot(xrange,linpred,'r','linewidth',2)
plot(xrange,yEst2,'r','linewidth',2)
plot(xrange,r2maxpred,'b','linewidth',2)

% Plot extrapolated dims on top
% plot([0;6e3],linpred(end).*ones(2,1),'r--','linewidth',2)
% plot([0;6e3],r2maxpred(end).*ones(2,1),'b--','linewidth',2)

% plot([0;xrange(31)],linpred(31).*ones(2,1),'r:','linewidth',2)
plot([0;xrange(31)],yEst2(31).*ones(2,1),'r:','linewidth',2)
plot([0;xrange(31)],r2maxpred(31).*ones(2,1),'b:','linewidth',2)

% plot([xrange(31);xrange(31)],[0,linpred(31)],'r:','linewidth',2)
plot([xrange(31);xrange(31)],[0,yEst2(31)],'r:','linewidth',2)
plot([xrange(31);xrange(31)],[0,r2maxpred(31)],'b:','linewidth',2)

xlim([0,2050])
%% Refit to individual mice
% var1 = Network_Rejection_Table.Signal_Size_WCM;
% var2 = Network_Rejection_Table.WCM_RejectionDn;
var1 = Network_Rejection_Table.Event_Signal_Size_WCM;
var2 = Network_Rejection_Table.Event_WCM_RejectionDn;
xdata = var1(keepers);
ydata = var2(keepers);

clear ax h 
xrange = linspace(0,6e3,100);
figure(105); clf;
ax(i) = subplot(3,3,1);
labels = {'N';'D'};
plot_rejection_pairs(xdata,ydata,labels,colour_ID(keepers),edgecolour(keepers,:),dotcolour(keepers,:))
% plot(xrange,linpred,'r','linewidth',2)
plot(xrange,yEst2,'r','linewidth',2)
plot(xrange,r2maxpred,'b','linewidth',2)
% title(['D1.8K:{\color{red}',num2str(round(linpred(31))),'}|{\color{blue}',num2str(round(r2maxpred(31))), '}. D6K:{\color{red}',num2str(round(linpred(end))),'}|{\color{blue}',num2str(round(r2maxpred(end))),'}'])

clear lindims_L23 lindims_B expdims_L23 expdims_B
lindims_L23(1,1) = round(yEst2(31)); %round(linpred(31));
expdims_L23(1,1) = round(r2maxpred(31));
lindims_B(1,1) = round(yEst2(end)); %round(linpred(end));
expdims_B(1,1) = round(r2maxpred(end));
%%
clear AICs BICs coeffs SS residuals
% IV = [0.5,0];  % [0.6557;0.0357]
for i = 1:numel(unique(colour_ID))
    
    these_As = find(animals(keepers) == i);
    this_A = keepers(these_As);
    xdata = var1(this_A);
    ydata = var2(this_A);
    
    ax(i+1) = subplot(3,3,i+1);
    h(i+1) = plot_rejection_pairs(xdata,ydata,labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
    
    [coeffs{i,1},SS(i,1),residuals{i,1}] = lsqcurvefit(linfun,[1;0], xdata, ydata,LB,UB,options); 
    [coeffs{i,2},SS(i,2),residuals{i,2}] = lsqcurvefit(expmaxAfun,[100;0], xdata, ydata,LB,UB,options);
    
    linpred = linfun(coeffs{i,1},xrange);
    B = xdata\ydata;
    yEst2 = xrange*B;
    
    r2maxpred = expmaxAfun(coeffs{i,2},xrange);
    

%     plot(xrange,linpred,'r','linewidth',2)
    plot(xrange,yEst2,'r','linewidth',2)
    plot(xrange,r2maxpred,'b','linewidth',2)
%     title(['D1.8K:{\color{red}',num2str(round(linpred(31))),'}|{\color{blue}',num2str(round(r2maxpred(31))), '}. D6K:{\color{red}',num2str(round(linpred(end))),'}|{\color{blue}',num2str(round(r2maxpred(end))),'}'])

    lindims_L23(1,i+1) = round(yEst2(31)); %round(linpred(31));
    expdims_L23(1,i+1) = round(r2maxpred(31));
    lindims_B(1,i+1) = round(yEst2(end)); %round(linpred(end));
    expdims_B(1,i+1) = round(r2maxpred(end));
    
    for j = 1:2
        AICs(i,j) = AICSS(SS(i,j),length(ydata),length(coeffs{i,j}));
        BICs(i,j) = BICSS(SS(i,j),length(ydata),length(coeffs{i,j}));
    end
end

linkaxes(ax)

% Bootstrap error bars TO DO

%% plot extrapolation dim stats
figure(106); clf
subplot(1,2,1);
plot([ones(8,1),2*ones(8,1)]',[lindims_L23(2:9);expdims_L23(2:9)],'linewidth',2,'color',[.7,.7,.7])
hold all
plot(ones(8,1),lindims_L23(2:9),'r.','markersize',20)
plot(2*ones(8,1),expdims_L23(2:9),'b.','markersize',20)
plot([1,2],[lindims_L23(1);expdims_L23(1)],'linewidth',2,'color',[.2,.2,.2])
plot([1,2],[lindims_L23(1);expdims_L23(1)],'k.','markersize',20)
xlim([0.5,2.5]); axis square; ylabel('Dimensions')
set(gca,'Xtick', 1:2,'Xticklabel',{'Linear';'Exponential'}); title('L2/3')

subplot(1,2,2);
plot([ones(8,1),2*ones(8,1)]',[lindims_B(2:9);expdims_B(2:9)],'linewidth',2,'color',[.7,.7,.7])
hold all
plot(ones(8,1),lindims_B(2:9),'r.','markersize',20)
plot(2*ones(8,1),expdims_B(2:9),'b.','markersize',20)
plot([1,2],[lindims_B(1);expdims_B(1)],'linewidth',2,'color',[.2,.2,.2])
plot([1,2],[lindims_B(1);expdims_B(1)],'k.','markersize',20)
xlim([0.5,2.5]); axis square; ylabel('Dimensions')
set(gca,'Xtick', 1:2,'Xticklabel',{'Linear';'Exponential'}); title('Barrel column')
%% Separate figure for each mouse. Plot slope for each on top
% var2 = Network_Rejection_Table.WCM_RejectionDn;
var2 = Network_Rejection_Table.eig90;
var3 = Network_Rejection_Table.Network_Size;
xrange = linspace(1,1e4,100)';

animals = Network_Rejection_Table.Animal;

% labels = {'Retained N';'RejectionDn'};
labels = {'Retained N';'Eig 90'};

% plot_rejection_pairs(va1,va2,labels);
figure(3); clf;
clear Betas
for i = 1:numel(unique(colour_ID))
    subplot(3,3,i);
    this_A = find(animals == i);
    % labels = {'Retained N';'WCM_{RejectionDn}'};
    
    h(2) = plot_rejection_pairs(var3(this_A),var2(this_A),labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
    
    [x,ix] = sort(var3(this_A));
    y = var2(this_A(ix));
    X = [ones(length(x),1) x];
    B2 = X\y;
    Betas(i) = B2(2);
    yEst2 = X*B2;
    plot(x,yEst2,'color',edgecolour(this_A(1),:),'linewidth',2)
    
    yExt2 = [ones(100,1), xrange]*B2;
    title(['D10K:',num2str(round(yExt2(end))),' Slope:',num2str(B2(2),3)])
    xlim([0,2000])
%     ylim([0,60])
    ylim([0,1400])
    %     plot(xrange,yExt2,'color',edgecolour(this_A(1),:),'linewidth',2)
end

%% Plot D vs N for recordings within a session
% var2 = Network_Rejection_Table.WCM_RejectionDn;
var2 = Network_Rejection_Table.eig90;
var3 = Network_Rejection_Table.Network_Size;

Sess_Betas = [];
for A = 1:8
    figure(4); clf
    this_A = find(Network_Rejection_Table.Animal == A);
    sess_list = unique(Network_Rejection_Table.Session(this_A));
    n = 0;
    for S = 1:numel(sess_list)
   
        
        this_S = find(Network_Rejection_Table.Session(this_A) == sess_list(S));
        if numel(this_S)>=5
            n = n + 1;
            
            ax(n) = subplot(4,3,n)
            this_plot = this_A(this_S);
            labels = {'',''};
            h(4) = plot_rejection_pairs(var3(this_plot),var2(this_plot),labels,colour_ID(this_plot),edgecolour(this_plot,:),dotcolour(this_plot,:));
            %         title(Network_Rejection_Table.NetworkName{this_plot(1)}(1:19))
%             xlim([0,2000])
%             ylim([0,60])
            %         pause
            hold all
            
            
            y = var2(this_plot);
            X = [ones(length(this_plot),1) var3(this_plot)];
            B2 = X\y;
            Sess_Betas = [Sess_Betas,B2(2)];
            yEst2 = X*B2;
            plot(var3(this_plot),yEst2,'k','linewidth',2)
            
        end 
    end
    linkaxes(ax);
    enw = Network_Rejection_Table.NetworkName{this_plot(1)}(1:8);
    suptitle(enw)
    print(['Figures/noise_rejection/dimensionality_preround/Each_session_N_eig90_',enw],'-dpdf','-bestfit')
   
%     pause
end

%% model comparison plot
Dsets = {SS,AICs,BICs};
meth_names = {'Sum of squares';'AIC';'BIC'};
clf
for i = 1:3
    subplot(2,3,i);
    plot(Dsets{i}(:,1),'r','linewidth',2)
    hold all
    plot(Dsets{i}(:,2),'b','linewidth',2)
%     legend('linear','rise to max')
    xlabel('Mice')
    ylabel(meth_names{i})
    title(meth_names{i})
    axis square
    
    subplot(2,3,i+3);
    plot([ones(8,1),2*ones(8,1)]',[Dsets{i}(:,1),Dsets{i}(:,2)]','k','linewidth',2)
    hold all
    plot(ones(8,1),Dsets{i}(:,1),'r.','markersize',20)
    plot(2*ones(8,1),Dsets{i}(:,2),'b.','markersize',20)
    xlim([0,3])
    set(gca,'Xtick',[1,2],'Xticklabel',{'Linear';'Rise2Max'})
%     xtickangle(30)
%     legend('linear','rise to max')
%     xlabel('Mice')
    ylabel(meth_names{i})
    title(meth_names{i})
    axis square
    
end
suptitle('Red: Linear, Blue: Rise to max')
% legend('linear','rise to max','location','bestoutside')
%% Organised data so you get D vs N for recordings within a session, per mouse
% var2 = Network_Rejection_Table.WCM_RejectionDn;
var2 = Network_Rejection_Table.eig90;
var3 = Network_Rejection_Table.Network_Size;

All_betas = [];
All_animals = [];
All_intercepts = [];

figure(5); clf;
for A = 1:8
    Xdata = [];
    Ydata = [];
    Sdata = [];
    this_A = find(Network_Rejection_Table.Animal == A);
    sess_list = unique(Network_Rejection_Table.Session(this_A));
    for S = 1:numel(sess_list)
        
        this_S = find(Network_Rejection_Table.Session(this_A) == sess_list(S));
        this_plot = this_A(this_S);
        Sdata = [Sdata,S*ones(1,numel(this_S))];
        [sortedx,ix] = sort(var3(this_plot));
        Xdata = [Xdata;sortedx];
        Ydata = [Ydata;var2(this_plot(ix))];
        
    end
    
    
    %% Plot somehow
    ix = 1:numel(Xdata);
    % figure(5); clf;
    ax(A) = subplot_tight(8,1,A)
    cmap = lines(numel(unique(Sdata)));
    max_so_far = 0;
    Sess_Betas = [];
    Sess_intercepts = [];
    n = 0;
    for i = 1:numel(unique(Sdata))
        this_S = find(Sdata == i);
        if numel(this_S) >= 5
            %     plot(ix(this_S),Xdata(this_S),'o','markersize',5);
            plot(Xdata(this_S)+max_so_far,Ydata(this_S),'.','markersize',10,'color',colours(A,:));
            
            hold all
            
            
            y = Ydata(this_S);
            X = [ones(length(this_S),1) Xdata(this_S)];
            B2 = X\y;
            Sess_Betas = [Sess_Betas,B2(2)];
            Sess_intercepts = [Sess_intercepts,B2(1)];
            yEst2 = X*B2;
            plot(Xdata(this_S)+max_so_far,yEst2,'color',colours(A,:),'linewidth',2)
            
            
            max_so_far = max_so_far + Xdata(this_S(end));
            n = n +1;
        end
    end
    
    All_betas = [All_betas,Sess_Betas];
    All_intercepts = [All_intercepts,Sess_intercepts];
    All_animals = [All_animals,A*ones(1,n)];
    
    % ylim([0,60])
end
linkaxes(ax)

%% Plot session betas (sessions with 4 or more subvolumes)
figure(6); clf;
for i = 1:8
    this_A = find(All_animals == i);
    plot(i*ones(numel(this_A)),All_betas(this_A),'.','markersize',10,'color',colours(i,:))
    hold all
    
end
xlim([0,9])
% ylim([0,0.05])
ylabel('Regression Betas')


%% Plot regression slopes
figure(7); clf; hold all;
xrange = [0;1000];
for i = 1:numel(All_betas)
    yExt2 = [ones(2,1), xrange]*[All_intercepts(i);All_betas(i)];
    plot(xrange,yExt2,'color',colours(All_animals(i),:),'linewidth',2)
end



%% Significantly contributing neurons vs dimensionality
var1 = Network_Rejection_Table.WCM_Dn;
var2 = Network_Rejection_Table.WCM_RejectionDn;
% var3 = Network_Rejection_Table.Config_Dn;
% var4 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(3); clf;
subplot(1,2,1)
labels = {'WCM_{Dn}';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

% subplot(1,2,2)
% labels = {'Config_{Dn}';'Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)
% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')

%% Significantly contributing neurons vs N
var1 = Network_Rejection_Table.N;
var2 = Network_Rejection_Table.WCM_Dn;
% var3 = Network_Rejection_Table.N;
% var4 = Network_Rejection_Table.Config_Dn;


% plot_rejection_pairs(va1,va2,labels);
figure(4); clf;
subplot(1,2,1)
labels = {'N';'WCM_{Dn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

% subplot(1,2,2)
% labels = {'N';'Config_{Dn}'};
% h(2) = plot_rejection_pairs(var3,var4,labels,colour_ID,edgecolour,dotcolour)
% 
% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northwest')

%% Correlates with recording duration
var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
% var3 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(5); clf;
subplot(1,2,1)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)

% subplot(1,2,2)
% labels = {'Time (samples)';'Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','northeast')

%% Same but split to learning sessions vs non learning sessions

var1 = Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
% var3 = Network_Rejection_Table.Config_RejectionDn;

l = find(Network_Rejection_Table.Learning == 1);
nl = 1:numel(var1);
nl(l) = [];

% plot_rejection_pairs(va1,va2,labels);
figure(6); clf;
subplot(2,2,1)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1(l),var2(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

% subplot(2,2,2)
% labels = {'Time (samples)';'Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var1(l),var3(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(2,2,3)
labels = {'Time (samples)';'WCM_{RejectionDn}'};
h(3) = plot_rejection_pairs(var1(nl),var2(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

% subplot(2,2,4)
% labels = {'Time (samples)';'Config_{RejectionDn}'};
% h(4) = plot_rejection_pairs(var1(nl),var3(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))

%% Correlates with recording duration * N (i.e. number of datapoints)
var1 = Network_Rejection_Table.N .* Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
% var3 = Network_Rejection_Table.Config_RejectionDn;


% plot_rejection_pairs(va1,va2,labels);
figure(7); clf;
subplot(1,2,1)
labels = {'N * T';'WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour)
 
% subplot(1,2,2)
% labels = {'N * T';'Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var1,var3,labels,colour_ID,edgecolour,dotcolour)

%% Same but split to learning sessions vs non learning sessions

var1 = Network_Rejection_Table.N .* Network_Rejection_Table.T;
var2 = Network_Rejection_Table.WCM_RejectionDn;
% var3 = Network_Rejection_Table.Config_RejectionDn;

l = find(Network_Rejection_Table.Learning == 1);
nl = 1:numel(var1);
nl(l) = [];

% plot_rejection_pairs(va1,va2,labels);
figure(8); clf;
subplot(1,2,1)
labels = {'N * T';'Rejection Dn'};
h(1) = plot_rejection_pairs(var1(l),var2(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

% subplot(2,2,2)
% labels = {'N * T';'Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var1(l),var3(l),labels,colour_ID(l),edgecolour(l,:),dotcolour(l,:))

subplot(1,2,2)
labels = {'N * T';'Rejection Dn'};
h(3) = plot_rejection_pairs(var1(nl),var2(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))
 
% subplot(2,2,4)
% labels = {'N * T';'Config_{RejectionDn}'};
% h(4) = plot_rejection_pairs(var1(nl),var3(nl),labels,colour_ID(nl),edgecolour(nl,:),dotcolour(nl,:))


% legend('Peron','Calcium','Mouse 2','Mouse 3','Mouse 4','Mouse 5','location','best')
xlim([0,8e6])

%% Compare events D to calcium D
% var1 = Network_Rejection_Table.Network_Size;
% var2 = Network_Rejection_Table.Event_Network_Size;
var1 = Network_Rejection_Table.WCM_RejectionDn;
var2 = Network_Rejection_Table.Event_WCM_RejectionDn;

% labels = {'Calcium Network Size'; 'Event Network Size'};
labels = {'Calcium D'; 'Events D'};

% plot_rejection_pairs(va1,va2,labels);
figure(201); clf;
% subplot(1,2,1)
h(1) = plot_rejection_pairs(var1,var2,labels,colour_ID,edgecolour,dotcolour);
% plot([0,2500],[0,2500],'k--')

%% Feb 2019. Dimensionality of same neurons (SVIDs) over time

% Assigning new SVID to the non-learning subset via within mouse network
% size (happens to be unique in this dataset)
svid = 0;
clear this_A u_sv
for A = 1:8
    this_A = find(Network_Rejection_Table.Animal == A);
    u_sv = unique(Network_Rejection_Table.N(this_A),'stable');
    for i = 1:numel(u_sv)
        svid = svid + 1;
        these_svs = this_A(find(Network_Rejection_Table.N(this_A) == u_sv(i)));
        Network_Rejection_Table.SVID(these_svs) = svid;
    end
    
end

%% For each SVID, plot N by WCM_RejectionDn connected by a line
n_svid = hist(Network_Rejection_Table.SVID,svid);

rep_sv = find(n_svid>=4);

% unroll id of retained svids
kept_sv = [];
sv_unroll = [];
for i = 1:numel(rep_sv)
    x = find(Network_Rejection_Table.SVID ==rep_sv(i));
    kept_sv = [kept_sv;rep_sv(i)*ones(numel(x),1)];
    sv_unroll = [sv_unroll;x];
end

var1 = Network_Rejection_Table.Network_Size;
var2 = Network_Rejection_Table.WCM_RejectionDn;
labels = {'N','D'};
figure(301); clf; hold all
h(1) = plot_rejection_pairs_lines(var1(sv_unroll),var2(sv_unroll),labels,kept_sv,edgecolour(sv_unroll,:),dotcolour(sv_unroll,:));


%% OLDER (PRE COSYNE ABSTRACT) STUFF

%% NOT INTERESTING:
% % Repeat N*T but Separate figure for each mouse. Plot slope for each on top
% var1 = Network_Rejection_Table.N .* Network_Rejection_Table.T;
% var2 = Network_Rejection_Table.WCM_RejectionDn;
% % var2 = Network_Rejection_Table.eig90;
% xrange = linspace(1,1e4,100)';
% 
% animals = Network_Rejection_Table.Animal;
% 
% labels = {'N*T';'RejectionDn'};
% % labels = {'Retained N';'Eig 90'};
% 
% % plot_rejection_pairs(va1,va2,labels);
% figure(9); clf;
% clear Betas
% for i = 1:numel(unique(colour_ID))
%     subplot(3,3,i);
%     this_A = find(animals == i);
% 
% l = find(Network_Rejection_Table.Learning(this_A) == 1);
% nl = 1:numel(this_A);
% nl(l) = [];
% 
% labels = {'N * T';'Rejection Dn'};
% % h(2) = plot_rejection_pairs(var3(this_A),var2(this_A),labels,colour_ID(this_A),edgecolour(this_A,:),dotcolour(this_A,:))
%    
% h(3) = plot_rejection_pairs(var1(this_A(nl)),var2(this_A(nl)),labels,colour_ID(this_A(nl)),edgecolour(this_A(nl),:),dotcolour(this_A(nl),:))
% 
% % xlim([0,8e6])
% % ylim([0,60])
%     
%     
% %     [x,ix] = sort(var3(this_A));
% %     y = var2(this_A(ix));
% %     X = [ones(length(x),1) x];
% %     B2 = X\y;
% %     Betas(i) = B2(2);
% %     yEst2 = X*B2;
% %     plot(x,yEst2,'color',edgecolour(this_A(1),:),'linewidth',2)
% %     
% %     yExt2 = [ones(100,1), xrange]*B2;
% %     title(['D10K:',num2str(round(yExt2(end))),' Slope:',num2str(B2(2),3)])
% %     xlim([0,2000])
% % %     ylim([0,60])
% %     ylim([0,1400])
% end


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
% var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess);
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);
var6 = Network_Rejection_Table.Network_Size(L_sess) .* Network_Rejection_Table.T(L_sess);
var7 = Network_Rejection_Table.WCM_RejectionDn(L_sess) ./ Network_Rejection_Table.Network_Size(L_sess);
var8 = Network_Rejection_Table.eig90(L_sess);


% plot_rejection_pairs(va1,va2,labels);
figure(15); clf;

subplot(2,3,1)
labels = {'Pnolick','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var3,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(2,3,2)
labels = {'Pcorrect','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs(var4,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(2,3,3)
labels = {'Session','WCM_{RejectionDn}'};
h(1) = plot_rejection_pairs_lines(var5,var1,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

% subplot(3,3,4)
% labels = {'Pnolick','Config_{RejectionDn}'};
% h(2) = plot_rejection_pairs(var3,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))
% 
% subplot(3,3,5)
% labels = {'Pcorrect','Config_{RejectionDn}'};
% h(1) = plot_rejection_pairs(var4,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))
% 
% subplot(3,3,6)
% labels = {'Session','Config_{RejectionDn}'};
% h(1) = plot_rejection_pairs_lines(var5,var2,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(2,3,4)
labels = {'Pnolick','Eig90'};
h(1) = plot_rejection_pairs(var3,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(2,3,5)
labels = {'Pcorrect','Eig90'};
h(1) = plot_rejection_pairs(var4,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

subplot(2,3,6)
% labels = {'Session','WCM_{RejectionDn}./Network Size'};
labels = {'Session','Eig90'};
h(1) = plot_rejection_pairs_lines(var5,var8,labels,colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))

%% Loop version of above

L_sess = find(Network_Rejection_Table.Learning); %find(colour_ID>2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
% var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess);
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.Session(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

labelA = {'Pnolick';'Pcorrect';'Session'};
labelB = {'WCM_{RejectionDn}';'Eig90'}; % {'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'};
VARA = {var3;var4;var5}
VARB = {var1;var8}; % var1;var2;var8

figure(16); clf;
hold all
n = 0;
for i = 1:2;
    for j = 1:3;
        n = n + 1;
        subplot(2,3,n);
        h(1) = plot_rejection_pairs_lines(VARA{j},VARB{i},{labelA{j},labelB{i}},colour_ID(L_sess),edgecolour(L_sess,:),dotcolour(L_sess,:))
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

% Save learning and sv_order
save('Results_reject_preround/Network_Rejection_Table_wStats_preround3.mat','Network_Rejection_Table')

%% Plot N by datenum - quick way of looking at the experiment
figure(16);
clf
for a = 1:8
    this_a = find(Network_Rejection_Table.Animal(this_m)==a);
    plot(Network_Rejection_Table.datenum(this_m(this_a)),Network_Rejection_Table.N(this_m(this_a)),'.');
    hold all;
end



this_a = find(Network_Rejection_Table.Animal == 2);

%% Plot non-learning sessions vs session order
L_sess = find(Network_Rejection_Table.Learning); % find(colour_ID==2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
% var2 = Network_Rejection_Table.Config_RejectionDn(L_sess);

var3 = Network_Rejection_Table.Pnolick(L_sess);
var4 = Network_Rejection_Table.Pcorrect(L_sess);
var5 = Network_Rejection_Table.sv_order(L_sess);

var8 = Network_Rejection_Table.eig90(L_sess);

subvols = Network_Rejection_Table.N_unique(L_sess);

labelA = {'Pnolick';'Pcorrect';'SV_{order}'};
labelB = {'WCM_{RejectionDn}';'Eig90'}; % 'WCM_{RejectionDn}';'Config_{RejectionDn}';'Eig90'
VARA = {var3;var4;var5};
VARB = {var1;var8}; % var1;var2;var8

figure(20); clf;
hold all
n = 0;
for i = 1:2
    for j = 1:3
        n = n + 1;
        subplot(2,3,n);
        h(1) = plot_rejection_pairs_lines(VARA{j},VARB{i},{labelA{j},labelB{i}},subvols,edgecolour(L_sess,:),dotcolour(L_sess,:));
        
    end
end
suptitle('Calcium')

%% Look at each animal separately. Rejection Dn and Eig90 vs SV order and Pcorrect
L_sess = find(Network_Rejection_Table.Learning); % find(colour_ID==2); % find(colour_ID > 6); % find(colour_ID>2 & colour_ID<7); %
var1 = Network_Rejection_Table.WCM_RejectionDn(L_sess);
var2 = Network_Rejection_Table.eig90(L_sess);
var3 = Network_Rejection_Table.Pcorrect(L_sess);
var4 = Network_Rejection_Table.sv_order(L_sess);
labelA = {'Pcorrect';'SV_{order}'};
labelB = {'Rejection Dn';'Eig90'}; 
VARA = {var3;var4};
VARB = {var1;var2}; % var1;var2;var8

subvols = Network_Rejection_Table.N_unique(L_sess);
unique_Lsv = unique(subvols); 

figure(21); clf;
hold all
n = 0;
for i = 1:2
    for j = 1:2
        for A = 1:4
            this_LS = find(subvols == unique_Lsv(A));
            n = n + 1;
            subplot(4,4,n);
            vA = VARA{j}(this_LS);
            vB = VARB{i}(this_LS);
            h(1) = plot_rejection_pairs(vA,vB,{labelA{j},labelB{i}},colour_ID(L_sess(this_LS)),edgecolour(L_sess(this_LS),:),dotcolour(L_sess(this_LS),:));
            hold all
                        
            y = vB;
            X = [ones(length(vA),1) vA];
            B2 = X\y;

            yEst = X*B2;
            plot(vA,yEst,'color',edgecolour(L_sess(this_LS(1)),:),'linewidth',2)
            
            % Stats
            mdl = fitlm(vA,vB);
            title(['R:',num2str(mdl.Rsquared.Adjusted,3),' p:',num2str(mdl.Coefficients.pValue(2),3)])
            
            
            
        end
    end
end



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
for i = 2:5
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
    
%     h = plot(var1(these_d),var2(these_d),'color',dotcolour(these_d(1),:),'linewidth',2);
    h = plot(var1(these_d),var2(these_d),'color',edgecolour(these_d(1),:),'linewidth',2);
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

%% zscore function that normalises the whole dataset to values of a subset
function z = local_z(x,keepers)
mu = mean(x(keepers));
sigma = std(x(keepers));
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end