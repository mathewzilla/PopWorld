% conversion_test.m
% Script to load example data and run different versions of conversion in
% poissonSparseWCM


%% Compare changing the conversion factor to 1000, vs thresholding the
% matrix to remove links weaker than the lowest conversion factor
%
% Another option is to use ceil instead of round in the conversion code
% (but only looking at non-zero links):
%
% A_int2 = zeros(n,n);
% ixLinks = find(A>0);
% A_int2(ixLinks) = ceil(A(ixLinks)*conversion);

load('/Volumes/Extras/197522_rejection/Rejected_an197522_2013_02_20_data_s_sv_1.mat')
Data_A = Data;
load('/Volumes/Extras/197522_rejection/Rejected_an197522_2013_02_21_data_s_sv_1.mat')
Data_B = Data;
%% run poissonSparseWCM with conversion factor of 100 and 1000, on Data_A and Data_B
for run = 1
    pars.C = 100;
    % A = Data_B.A;
    A = Data_A.A;
    
    [PSWCM_100.E,PSWCM_100.D,PSWCM_100.V,PSWCM_100.ExpA,PSWCM_100.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);
    
    pars.C = 1000;
    
    [PSWCM_1K.E,PSWCM_1K.D,PSWCM_1K.V,PSWCM_1K.ExpA,PSWCM_1K.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);
    
    %% Plot the two ExpAs
    figure(5); clf;
    subplot(2,2,1);
    imagesc(PSWCM_100.ExpA);
    colorbar
    title('Conversion = 100')
    subplot(2,2,3);
    plot(sum(PSWCM_100.ExpA));
    
    
    subplot(2,2,2);
    imagesc(PSWCM_1K.ExpA);
    colorbar
    title('Conversion = 1000')
    subplot(2,2,4);
    plot(sum(PSWCM_1K.ExpA));
    
    
    suptitle('Feb 20 ExpA A_{int}')
end
%% Looking at distribution of links. How many zeros vs effectively zeros?
conversion = 100;
n = size(A,1);
idx = find(triu(ones(size(A)),1));
figure(6); clf;
plot(sort(A(idx),'descend'),'linewidth',2)
hold all
% ixLinks = find(A(idx) > 0);
% plot(sort(A(idx(ixLinks)),'descend'))

% Version 1 of A_int, then divided by 100
A_int = round(A*conversion);
plot(sort(A_int(idx),'descend')/100,'linewidth',2)

% Version 2 of A_int (preserving all links), then divided by 100
A_int2 = zeros(n,n);
ixLinks = find(A>0);
A_int2(ixLinks) = ceil(A(ixLinks)*conversion);
plot(sort(A_int2(idx),'descend')/100,'linewidth',2)


% Version 3 of A_int (setting links below 1 to 0), then divided by 100
A_int3 = zeros(n,n);
ixLinks = find(A>0);
A_int3(ixLinks) = floor(A(ixLinks)*conversion);
plot(sort(A_int3(idx),'descend')/100,'linewidth',2)
xlabel('Sorted links')
ylabel('Link strength (Corr coef)')

% V4, after round(A) set all values of A that are zero in A_int to zero
tiny_links2 = find(A_int(ixLinks)==0);
A_z = A;
A_z(ixLinks(tiny_links2)) = 0;
plot(sort(A_z(idx),'descend'),'linewidth',2)

% % V5, where A is thresholded
% ixLinks = find(A>0);
% C_thresh = max(A(:))/conversion;
% tiny_links = find(A(ixLinks)<(C_thresh));
% L_z = ixLinks(tiny_links);
% A_z = A;
% A_z(L_z) = 0;
% plot(sort(A_z(idx),'descend'),'linewidth',2)
%
% % V6, A is thresholded then converted
% A_z_int = round(A_z*conversion);
% plot(sort(A_z_int(idx),'descend')/100,'linewidth',2)

legend('A','A_{round}','A_{ceil}','A_{floor}','A_z')

axis square


% xlim([0,1e6]); ylim([0,0.1])
% xlim([0,100]); ylim([0.3,0.7])

%% Bar plot of link histograms
figure(7); clf;
cmap = colormap(lines);

[h,b] = hist(A(idx),100);
stairs(b,h,'linewidth',2)

hold all

[h,b] = hist(A_int(idx)/100,100);
stairs(b,h,'linewidth',2)

[h,b] = hist(A_z(idx),100);
stairs(b,h,'linewidth',2,'color',cmap(5,:))

legend('A','A_{round}','A_z')

%% repeat conversion = 100/1000 run but with A_int2 code (change round to ceil in poissonSparseWCM)
pars.C = 100;
% A = Data_B.A;
% A = Data_A.A;

[PSWCM_100two.E,PSWCM_100two.D,PSWCM_100two.V,PSWCM_100two.ExpA,PSWCM_100two.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);

pars.C = 1000;

[PSWCM_1Ktwo.E,PSWCM_1Ktwo.D,PSWCM_1Ktwo.V,PSWCM_1Ktwo.ExpA,PSWCM_1Ktwo.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);

%% Plot the two ExpAs
figure(7); clf;
subplot(2,2,1);
imagesc(PSWCM_100two.ExpA);
colorbar
title('Conversion = 100')
subplot(2,2,3);
plot(sum(PSWCM_100two.ExpA));


subplot(2,2,2);
imagesc(PSWCM_1Ktwo.ExpA);
colorbar
title('Conversion = 1000')
subplot(2,2,4);
plot(sum(PSWCM_1Ktwo.ExpA));


suptitle('Feb 20 ExpA A_{int2}')


%% repeat conversion = 100/1000 run but with A_int3 code (change round to floor in poissonSparseWCM)
pars.C = 100;
% A = Data_B.A;
% A = Data_A.A;

[PSWCM_100three.E,PSWCM_100three.D,PSWCM_100three.V,PSWCM_100three.ExpA,PSWCM_100three.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);

pars.C = 1000;

[PSWCM_1Kthree.E,PSWCM_1Kthree.D,PSWCM_1Kthree.V,PSWCM_1Kthree.ExpA,PSWCM_1Kthree.ALL] = poissonSparseWCM(A,pars.N,pars.C,optionsModel);

%%
save(['/Volumes/Extras/197522_rejection/Rejection_tests/conv_test_a_02_21.mat'],'PSWCM_100','PSWCM_1K','PSWCM_100two','PSWCM_1Ktwo','PSWCM_100three','PSWCM_1Kthree')
%% Plot the two ExpAs
figure(9); clf;
subplot(2,2,1);
imagesc(PSWCM_100three.ExpA);
colorbar
title('Conversion = 100')
subplot(2,2,3);
plot(sum(PSWCM_100three.ExpA));


subplot(2,2,2);
imagesc(PSWCM_1Kthree.ExpA);
colorbar
title('Conversion = 1000')
subplot(2,2,4);
plot(sum(PSWCM_1Kthree.ExpA));


suptitle('Feb 20 ExpA A_{int3}')

%% using poissonSparseWCMReal.m
% to run code again but set links < 0.01 to 0 ahead of analysis
% load('/Volumes/Extras/197522_rejection/Rejected_an197522_2013_02_20_data_s_sv_1.mat')
load('/Volumes/Extras/197522_rejection/Rejected_an197522_2013_02_21_data_s_sv_1.mat')
Data_A = Data;

pars.C = 100;
% A = Data_B.A;
A = Data_A.A;

[Pre_conv.E,Pre_conv.D,Pre_conv.V,Pre_conv.ExpA,Pre_conv.ALL,Pre_conv.A_z,Pre_conv.L_z] = poissonSparseWCMReal(A,pars.N,pars.C,optionsModel);

pars.C = 1000;

[Pre_conv_1K.E,Pre_conv_1K.D,Pre_conv_1K.V,Pre_conv_1K.ExpA,Pre_conv_1K.ALL,Pre_conv_1K.A_z,Pre_conv_1K.L_z] = poissonSparseWCMReal(A,pars.N,pars.C,optionsModel);

save(['/Volumes/Extras/197522_rejection/Rejection_tests/preconv_test_a_02_21.mat'],'Pre_conv','Pre_conv_1K')

%% Plot the two ExpAs with poissonSparseWCMReal
figure(9); clf;
subplot(2,2,1);
imagesc(Pre_conv.ExpA);
colorbar
title('Conversion = 100')
subplot(2,2,3);
plot(sum(Pre_conv.ExpA));


subplot(2,2,2);
imagesc(Pre_conv_1K.ExpA);
colorbar
title('Conversion = 1000')
subplot(2,2,4);
plot(sum(Pre_conv_1K.ExpA));


suptitle('Feb 21 ExpA A_z')

%% Compute dimensionality from the different versions.
Data = Pre_conv_1K; % Pre_conv
% Data.ALL = Data.A;
Data.A = Data.A_z;
pars.C = 100;

fname = 'Pre_conv_1K_21_02'; % 'Pre_conv_1K_20_02';

%% Load ceil version etc if needed
load('/Volumes/Extras/197522_rejection/Rejection_tests/conv_test_a_02_20.mat')
load('/Volumes/Extras/197522_rejection/Rejected_an197522_2013_02_20_data_s_sv_1.mat')

fname = 'Floor_20_02'; %'Ceil_20_02'; %'Ceil_1K_20_02'; % 'Ceil_20_02' % 'Round_21_02' % 'Floor_21_02'

% Order is ..:round, two:ceil, three:floor
Data.ExpA = PSWCM_100three.ExpA;
Data.E = PSWCM_100three.E;
Data.V = PSWCM_100three.V;
%% decompose nodes into signal and noise
B = Data.A - Data.ExpA;  % modularity matrix using chosen null model

% find low-dimensional projection
[Data.Dspace,Data.ixpos,Data.Dn,Data.EigEst,Data.Nspace,Data.ixneg,Data.Dneg,Data.NEigEst] = LowDSpace(B,Data.E,pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors

% compute dimensions based on just positive eigenvalues
egs = eig(B);  % eigenspectra of data modularity matrix
egs = sort(egs,'descend'); % sort eigenvalues into descending order
Data.PosDn = sum(egs > pars.eg_min);

% node rejection within low-dimensional projection
Rejection = NodeRejection(B,Data.E,pars.I,Data.V,optionsReject); % N.B. also calls LowDSpace function to find projections

% new signal matrix
Data.Asignal = Data.A(Rejection.ixSignal,Rejection.ixSignal);

% connected signal matrix: find largest component, and use that - store
% others
[Data.Asignal_comp,Data.ixRetain,Data.SignalComps,Data.SignalComp_sizes] = prep_A(Data.Asignal);
Data.ixSignal_comp = Rejection.ixSignal(Data.ixRetain);  % original node indices

% and then strip out leaves - nodes with single links
K = sum(Data.Asignal_comp);
ixLeaves = find(K==1); ixKeep = find(K > 1);

Data.ixSignal_Final = Data.ixSignal_comp(ixKeep);
Data.ixSignal_Leaves = Data.ixSignal_comp(ixLeaves);
Data.Asignal_final = Data.Asignal_comp(ixKeep,ixKeep);

%%

save(['/Volumes/Extras/197522_rejection/Rejection_tests/Rejected_' fname],'Rejection','Data','pars','optionsModel','optionsReject')

%% Summarise into a table

fnames = dir('/Volumes/Extras/197522_rejection/Rejection_tests/Rejected_*');

nF = numel(fnames);

netCtr = 0;
for iF = 1:nF
    
    if any(strfind(fnames(iF).name,'Rejected'))
        netCtr = netCtr + 1;
        result(netCtr).NetworkName = fnames(iF).name(10:end-4); % strip out 'Rejected' and .mat
        Data = {};
        load(['/Volumes/Extras/197522_rejection/Rejection_tests/',fnames(iF).name]);
        result(netCtr).Network_Size = numel(Data.ixRetain);
        result(netCtr).Signal_Size_WCM = numel(Data.ixSignal_Final);
        result(netCtr).WCM_Dn = Data.PosDn;
        result(netCtr).WCM_RejectionDn = Data.Dn;
%         result(netCtr).Config_Dn = Control.PosDn;
%         result(netCtr).Config_RejectionDn = Control.Dn;
        result(netCtr).Signal_Components = numel(Data.SignalComp_sizes);
    end
end

Network_Rejection_Table = struct2table(result);
save('/Volumes/Extras/197522_rejection/Rejection_tests/Scaling_Network_Rejection_Table','Network_Rejection_Table');

%% Make a plot
clear all
load('/Volumes/Extras/197522_rejection/Rejection_tests/Scaling_Network_Rejection_Table');

% Parse names and add variables for easier plotting
meth = {'Roun','Ceil','Floo','Pre_'};
for i = 1:height(Network_Rejection_Table)
    Name = Network_Rejection_Table.NetworkName{i};
    for j = 1:4
        if strcmp(Name(1:4),meth{j})
             Network_Rejection_Table.conv_meth(i) = j;
        end
    end
    
    D = Name(end-4:end-3);
    if strcmp(D,'20')
        Network_Rejection_Table.day(i) = 1;
    else
        Network_Rejection_Table.day(i) = 2;
    end
    
    C = Name(end-6);
    if strcmp(C,'K')
        Network_Rejection_Table.C(i) = 1000;
    else
        Network_Rejection_Table.C(i) = 100;
    end
    
end


%% Plot WCM_RejectionDn for each method
clf;
% cmap = colormap(lines);
% plot(ones(4,4),zeros(4,4),'.','markersize',20);
% hold all
cmap1 = colormap(redblue(7));
cmap = cmap1([1,3,7,5],:);

for i = 1:4
plot(1,0,'.','markersize',20,'color',cmap(i,:));
hold all
end

C_100 = find(Network_Rejection_Table.C == 100);
C_1000 = find(Network_Rejection_Table.C == 1000);
D_1 = find(Network_Rejection_Table.day(C_100) == 1);
D_2 = find(Network_Rejection_Table.day(C_100) == 2);

plotting = C_100(D_1);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(1,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(1,:))

plotting = C_1000(D_1);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(3,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(3,:))

plotting = C_100(D_2);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(2,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(2,:))

plotting = C_1000(D_2);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(4,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(4,:))

set(gca,'Xtick',1:4,'XTicklabel',{'Round','Ceil','Floor','Pre_{round}'});
ylabel('Dimensionality')

legend('20.02, C = 100','20.02, C = 1000','21.02, C = 100','21.02, C = 1000')

axis square
xlim([0.5,4.5])

%% Direct comparison of ceil and pre_round
clf;
cmap1 = colormap(redblue(7));
cmap = cmap1([1,3,7,5],:);
for i = 1:4
plot(1,0,'.','markersize',20,'color',cmap(i,:));
hold all
end


plotting = C_100(D_1);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
sorted_meth([1,3]) = [];
meth_order([1,3]) = [];
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(1,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(1,:));

plotting = C_1000(D_1);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
sorted_meth([1,3]) = [];
meth_order([1,3]) = [];
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(3,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(3,:))

plotting = C_100(D_2);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
sorted_meth([1,3]) = [];
meth_order([1,3]) = [];
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(2,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(2,:))

plotting = C_1000(D_2);
[sorted_meth,meth_order] = sort(Network_Rejection_Table.conv_meth(plotting));
sorted_meth([1,3]) = [];
meth_order([1,3]) = [];
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'linewidth',2,'color',cmap(4,:));
plot(sorted_meth,Network_Rejection_Table.WCM_RejectionDn(plotting(meth_order)),'.','markersize',20,'color',cmap(4,:))

set(gca,'Xtick',[2,4],'XTicklabel',{'Ceil','Pre_{round}'});
ylabel('Dimensionality')

legend('20.02, C = 100','20.02, C = 1000','21.02, C = 100','21.02, C = 1000')

axis square
xlim([1.5,4.5])
