%% Code to load and plot peron crcns clustering results

%% Look at cluster ID for data vs events versions
f = dir('example_clustering/*.mat');
for i = linspace(2,numel(f),numel(f)/2)
    load(['example_clustering/',f(i-1).name]); 
    subplot(1,2,1)
    plot(Full.QmaxCluster);
    title('Calcium')

    subplot(1,2,2)
    load(['example_clustering/',f(i).name]); 
    plot(Full.QmaxCluster);
    title('Peron events')
    pause
end