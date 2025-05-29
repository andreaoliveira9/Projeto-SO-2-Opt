function [bestScore, bestNodes, totalIterations, bestFoundTime] = GRASP_SNS(G, time, n, r, Cmax, seed)
    % GRASP_SNS GRASP algorithm for Server Node Selection with Cmax constraint
    % Outputs:
    %   bestScore - best average shortest path length
    %   bestNodes - best solution found
    %   totalIterations - number of local improvement loops performed
    %   bestFoundTime - time when best solution was found

    if nargin >= 6 && ~isempty(seed)
        rng(seed); % Definir seed do gerador de números aleatórios
    end

    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;

    globalStartTime = tic;
    D = distances(G);

    while true
        elapsed = toc(globalStartTime);
        if elapsed >= time
            break;
        end
    
        timeLeft = time - elapsed;
        if timeLeft < 0.01  % margem mínima para executar algo útil
            break;
        end
    
        currentNodes = GreedyRandomizedConstruction(G, D, n, r, Cmax);
        if toc(globalStartTime) >= time
            break;
        end

        if isempty(currentNodes) || length(unique(currentNodes)) < n
            continue;
        end

        [currentScore, maxSP] = PerfSNS(G, currentNodes);
        if maxSP > Cmax
            continue;
        end
        
        % Checa tempo restante antes da busca local
        if toc(globalStartTime) >= time
            break;
        end

        % --- Call local search function ---
        [currentNodes, currentScore, localIterations, ~] = ...
            LocalSearch_SA_HC(G, currentNodes, currentScore, Cmax, globalStartTime, time);

        totalIterations = totalIterations + localIterations;

        if currentScore < bestScore
            bestScore = currentScore;
            bestNodes = currentNodes;
            bestFoundTime = toc(globalStartTime);
        end
    end
end


function nodes = GreedyRandomizedConstruction(G, D, n, r, Cmax)
    % Generate a greedy randomized solution for the Server Node Selection problem
    % with Cmax constraint
    % Input:
    %   G - graph of the network
    %   D - precomputed distances matrix
    %   n - number of nodes to include in the solution
    %   r - size of the Restricted Candidate List (RCL)
    %   Cmax - maximum allowed shortest path length between any pair of server nodes
    % Output:
    %   nodes - array of selected nodes, empty if no feasible solution found
    
    numNodes = numnodes(G);
    nodes = zeros(1, n);
    remaining = 1:numNodes;
    
    % Start with a random node with good centrality
    centrality = zeros(1, numNodes);
    for i = 1:numNodes
        centrality(i) = sum(D(i, :));  % Lower is better (total distance to all other nodes)
    end
    
    % Select first node from RCL based on centrality
    [~, sortedIndices] = sort(centrality);
    rclSize = min(r, numNodes);
    rcl = sortedIndices(1:rclSize);
    selectedIdx = randi(rclSize);
    nodes(1) = rcl(selectedIdx);
    
    % Remove selected node from remaining
    remaining(remaining == nodes(1)) = [];
    
    % Flag to indicate if we can't find a valid solution
    invalidSolution = false;
    
    % For each subsequent node
    for k = 2:n
        % Calculate benefit of adding each remaining node
        benefit = zeros(1, length(remaining));
        validCandidates = false(1, length(remaining));
        
        for i = 1:length(remaining)
            candidate = remaining(i);
            tempSolution = [nodes(1:k-1), candidate];
            
            % Check if adding this node would violate the Cmax constraint
            % Only needed for k > 1
            if k > 1
                maxDistBetweenServers = 0;
                for a = 1:k-1
                    serverA = nodes(a);
                    % Distance to the candidate node
                    dist = D(serverA, candidate);
                    maxDistBetweenServers = max(maxDistBetweenServers, dist);
                end
                
                % Skip this candidate if it violates Cmax
                if maxDistBetweenServers > Cmax
                    continue;
                end
            end
            
            % Mark as valid candidate
            validCandidates(i) = true;
            
            % Calculate a proxy measure for solution quality
            % (Sum of distances from all nodes to their closest server)
            totalDistance = 0;
            for node = 1:numNodes
                if ~ismember(node, tempSolution)
                    dist = min(D(node, tempSolution));
                    totalDistance = totalDistance + dist;
                end
            end
            
            % We want to minimize, so negative benefit
            benefit(i) = -totalDistance;
        end
        
        % Check if we have any valid candidates
        if ~any(validCandidates)
            invalidSolution = true;
            break;
        end
        
        % Get indices of valid candidates
        validIndices = find(validCandidates);
        
        % Get benefits of valid candidates
        validBenefits = benefit(validCandidates);
        
        % Build RCL with top r valid candidates
        [~, sortedIndices] = sort(validBenefits, 'descend');
        
        % Determine size of RCL (minimum of r and valid remaining nodes)
        rclSize = min(r, length(validIndices));
        rcl = remaining(validIndices(sortedIndices(1:rclSize)));
        
        % Randomly select a node from the RCL
        selectedIdx = randi(rclSize);
        selectedNode = rcl(selectedIdx);
        
        % Add the selected node to the solution
        nodes(k) = selectedNode;
        
        % Remove selected node from remaining
        remaining = remaining(remaining ~= selectedNode);
    end
    
    % Return empty array if no valid solution found
    if invalidSolution
        nodes = [];
    end
end

function [bestNodes, bestScore, localIterations, improved] = ...
    LocalSearch_SA_HC(G, currentNodes, currentScore, Cmax, globalStartTime, maxTime)

    numNodes = numnodes(G);
    localIterations = 0;
    improved = true;
    bestNodes = currentNodes;
    bestScore = currentScore;

    while improved && toc(globalStartTime) < maxTime
        localIterations = localIterations + 1;
        improved = false;
        notSelected = setdiff(1:numNodes, bestNodes);
        bestNeighborScore = bestScore;
        bestSwap = [0, 0];

        for i = 1:length(bestNodes)
            for j = 1:length(notSelected)
                if toc(globalStartTime) >= maxTime
                    return;
                end

                neighborNodes = bestNodes;
                neighborNodes(i) = notSelected(j);

                [neighborScore, neighborMaxSP] = PerfSNS(G, neighborNodes);

                if neighborMaxSP > Cmax
                    continue;
                end

                if neighborScore < bestNeighborScore
                    bestNeighborScore = neighborScore;
                    bestSwap = [i, notSelected(j)];
                    improved = true;
                end
            end
        end

        if improved
            bestNodes(bestSwap(1)) = bestSwap(2);
            bestScore = bestNeighborScore;
        end
    end
end