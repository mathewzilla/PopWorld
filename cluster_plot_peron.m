%% Code to load and plot peron crcns clustering results

%% Look at cluster ID for data vs events versions
% f = dir('example_clustering/*.mat');
f = dir('Results_batch2/Clustered_an197522*');
f(end) = [];
for i = linspace(2,numel(f),numel(f)/2)
    clf
%     load(['example_clustering/',f(i-1).name]);
    load(['Results_batch2/',f(i-1).name]); 
    subplot(1,2,1)
    plot(Full.QmaxCluster);
    title('Calcium')

    
%     load(['example_clustering/',f(i).name]);
    load(['Results_batch2/',f(i).name]); 
    subplot(1,2,2)
    plot(Full.QmaxCluster);
    title('Peron events')
    suptitle(f(i).name)
    pause
end