% consensus_explorer.m
%
% Script to load two datasets and attempt clustering - trying to understand
% why one succeeds and the other fails

clear all

blnLabels = 0;      % write node labels? Omit for large networks
fontsize = 6;

clusterpars.nreps = 1000;
clusterpars.nLouvain = 5;

fname1 = 'an197522_2013_02_20_data_s_sv_1';
fname2 = 'an197522_2013_02_21_data_s_sv_1';

D1 = load(['Results_reject_preround/Rejected_', fname1,'.mat']);
    
D2 = load(['Results_reject_preround/Rejected_', fname2,'.mat']);

load(['Clustering_Results_preround/Network_Clustering_Table']);
%% Checking LowDSpace code works as expected
B1 = D1.Data.A - D1.Data.ExpA;
% [Data1.Dspace,Data1.ixpos,Data1.Dn,Data1.EigEst,Data1.Nspace,Data1.ixneg,Data1.Dneg,Data1.NEigEst] = LowDSpace(B1,D1.Data.E,D1.pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors
            
[V1,egs1] = eig(B1);  
egs1 = diag(egs1);
[egs1,ix] = sort(egs1,'descend');  
V1 = V1(:,ix);  

B2 = D2.Data.A - D2.Data.ExpA;
% [Data2.Dspace,Data2.ixpos,Data2.Dn,Data2.EigEst,Data2.Nspace,Data2.ixneg,Data2.Dneg,Data2.NEigEst] = LowDSpace(B2,D2.Data.E,D2.pars.I); % to just obtain low-dimensional projection; Data.Dn = number of retained eigenvectors
[V2,egs2] = eig(B2);  
egs2 = diag(egs2);
[egs2,ix] = sort(egs2,'descend');  
V2 = V2(:,ix);  

clf;
plot(egs1);
hold all
plot(egs2);
legend('egs1','egs2')


%% cluster - with noise rejection
% construct new null model
W1 = D1.Data.Asignal_final;
P1 = D1.Data.ExpA(D1.Data.ixSignal_Final,D1.Data.ixSignal_Final); % extract relevant part of null model
L1 = 1 + D1.Data.Dn;
M1 = 1 + D1.Data.Dn;

% W2 = D2.Data.Asignal_final;
% P2 = D2.Data.ExpA(D2.Data.ixSignal_Final,D2.Data.ixSignal_Final); % extract relevant part of null model
% L2 = 1 + D2.Data.Dn;
% M2 = 1 + D2.Data.Dn;

W = W1; P = P1; L = L1; M = M1;
% W = W2; P = P2; L = L2; M = M2;

%%
nreps = 1000;     % of each distance metric
nreps_cc = 50;  % CC gets its own smaller nreps
dims = 'all';   % use all embedding dimensions for each k-means clustering
blnExplore = 1;

% internal parameters

nIDs = size(W,1);     % number of nodes of the weight matrix
m = sum(sum(W))/2;    % number of unique links (or total unique weights)

blnConverged = 0;       % stopping flag
ctr = 1;                % iterations of consensus

% cluster signal network
B = W - P;          % initial modularity matrix, given data matrix W and specified null model P
[V,egs] = eig(B,'vector');
[~,ix] = sort(egs,'descend');    % sort into descending order
V = V(:,ix);                       % ditto the eigenvectors 

C = kmeansSweep(V(:,1:M-1),L,M,nreps,dims);  % find groups in embedding dimensions: sweep from L to M
% C = kmeansSweep(V(:,1:M),L,M,nreps,dims);  % find groups in embedding dimensions: sweep from L to M

for iQ = 1:size(C,2)
    Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering
end

%% Consensus clustering loop
while ~blnConverged
    % make consensus
    Allowed = (Q > 0);       % only take consensus using Q-positive clusterings...
    CCons = makeConsensusMatrix(C(:,Allowed));
    
    %% check convergence
    [blnConverged,grpscon] = CheckConvergenceConsensus(CCons);
    title(['Run ',num2str(ctr)]); drawnow;
%     keyboard
    %% sort out what to do next
    if blnConverged
        Qcon = computeQ(grpscon,B,m);  % compute Q, and exit
    else
        ctr = ctr+1;  % increment consensus iteration counter
        if ctr > 50
            % do escape if not converging
            warning('Did not converge in 50 iterations - exiting without consensus answer')
            grpscon = [];
            Qcon = 0;
            blnConverged = 0;
            return
        else
            
            % if a single size of groups was requested, and not exploring
            % just use the same number of groups as the original requested
            % set: find consensus in that space.
            if L == M && ~blnExplore
                disp('I am not exploring')
                
                T = sum(reshape(Allowed,nreps_cc,1+M-L));  % count how many at each K were retained
                [D,~,Mcons] = EmbedConsensusNull(CCons,'sweep',L:M,T);  % option 'expected' available as well as 'sweep'
                % do k-means sweep using D, restricted to original M groups
                % (thus M-1 dimensions)
                if Mcons >= M
                    C = kmeansSweep(D(:,1:M-1),L,M,nreps_cc,dims);  % find groups in embedding dimensions
                elseif isempty(D)
                    C = [];
                else
                    % use all of D if less than M returned
                    C = kmeansSweep(D,L,M,nreps_cc,dims);  % find groups in embedding dimensions 
                end
            end
            
%             % find upper limit of groups - replicate this code when using
%             null model for consensus 
%             [D,~,Mcons] = EmbedConsensusWishart(CCons);
%              % do k-means sweep using found M
%             C = kmeansSweep(D,L,Mcons,nreps,dims);  % find groups in embedding dimensions
           

            % keyboard
            if L~=M || blnExplore
                disp('I am exploring')
                reps = numel(Allowed)/(1+M-L);           % To catch difference in nreps on first run
                T = sum(reshape(Allowed,reps,1+M-L));  % count how many at each K were retained
                [D,~,Mcons] = EmbedConsensusNull(CCons,'sweep',L:M,T);  % option 'expected' available as well as 'sweep'
                M = Mcons;
                 % do k-means sweep using found M
                C = kmeansSweep(D,L,M,nreps_cc,dims);  % find groups in embedding dimensions
            end
            

%             % do Laplacian on consensus matrix, using original M
%             D = ProjectLaplacian(CCons,M);
%             keyboard
%             C = kmeansSweep(D,M,M,nreps,dims);  % find groups in embedding dimensions

            if isempty(C)
                warning('Consensus matrix projection is empty - exiting without consensus answer')
                grpscon = [];
                Qcon = 0;
                blnConverged = 0;               
                return
            else
             
                % compute Q
                Q = zeros(size(C,2),1);
                for iQ = 1:size(C,2)
                    Q(iQ) = computeQ(C(:,iQ),B,m); % compute modularity Q for each clustering using original modularity matrix
                end
            end
        end
    end
end






