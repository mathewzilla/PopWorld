function Q = mymodularity(W,C)
% Copied modularity computation code from LouvainCommunityUDnondeterm.m
    nG = numel(unique(C));
    N = numel(C);
    % null model
    inks = sum(W);
    outks = sum(W');
    m = sum(inks);  % number of edges
    P = (outks' * inks) ./ m;   % null model
    
    % community matrix
    S = zeros(N,nG);

    for loop = 1:nG
        S(:,loop) = (C == loop);
    end
    
    Q = trace(S' * (W-P) * S) ./ m;   % scaled by 1/2m...
    
%     % for debugging...
%     rQ = 0;
%     k = sum(W);
%     for i = 1:N
%         % compute node i's modularity for current assignment
%         Gcurr = find(C == C(i));    % all nodes in current group
%         nullmodel = k(i) .* k(Gcurr) ./ m; % expected number of links in current assignment
%         q_i = sum(W(i,Gcurr) - nullmodel);     
%         rQ = rQ + q_i;
%     end
%     rQ = rQ ./ m;   % running Q
    % keyboard
end