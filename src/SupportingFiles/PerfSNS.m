function [avgSP,maxSP]= PerfSNS(G,sNodes)
% [out1,out2]= PerfSNS(G,sNodes)
% OUTPUTS:
%   avgSP -   average shortest path length from each node to its closest
%             server node (returns -1 for invalid input data)
%   maxSP -   maximum shortest path length between any pair of server nodes
%             (returns -1 for invalid input data)
% INPUTS:
%   G -       graph of the network
%   sNodes -  a row array with server nodes
    
    nNodes= numnodes(G);
    if length(sNodes)<1
        avgSP= -1;
        maxSP= -1;
        return
    end
    if (max(sNodes)>nNodes || min(sNodes)<1 || length(unique(sNodes))<length(sNodes))
        avgSP= -1;
        maxSP= -1;
        return
    end
    clients= setdiff(1:nNodes,sNodes);
    dist= distances(G,sNodes,clients);
    if length(sNodes)>1
        avgSP= sum(min(dist))/nNodes;
        maxSP= max(max(distances(G,sNodes,sNodes)));
    else
        avgSP= sum(dist)/nNodes;
        maxSP= 0;
    end
end