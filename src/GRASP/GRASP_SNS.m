function [bestScore, bestNodes, totalIterations, bestFoundTime] = GRASP_SNS_with_Cmax(G, time, n, r, Cmax)
    % GRASP_SNS_with_Cmax GRASP algorithm for Server Node Selection with Cmax constraint
    % Input:
    %   G - graph of the network
    %   time - time to run the method (in seconds)
    %   n - number of nodes to include in a solution
    %   r - size of Restricted Candidate List (RCL)
    %   Cmax - maximum allowed shortest path length between any pair of server nodes
    % Output:
    %   bestScore - best average shortest path length
    %   bestNodes - list of selected nodes
    %   totalIterations - total number of iterations across all runs
    %   bestFoundTime - time when the best solution was found
    
    % Get the number of nodes in the graph
    numNodes = numnodes(G);
    
    % Initialize best solution
    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;
    
    % Record start time
    globalStartTime = tic;
    
    % Precompute distances matrix
    D = distances(G);
    
    % GRASP loop
    while toc(globalStartTime) < time
        % Greedy randomized construction phase
        currentNodes = GreedyRandomizedConstruction(G, D, n, r, Cmax);
        
        % Skip this iteration if no feasible solution was found during construction
        if isempty(currentNodes)
            continue;
        end
        
        % Evaluate the current solution
        [currentScore, maxSP] = PerfSNS(G, currentNodes);
        
        % Skip if not satisfying Cmax constraint
        if maxSP > Cmax
            continue;
        end
        
        % Local improvement loop (SA-HC)
        improved = true;
        localIterations = 0;
        
        while improved && (toc(globalStartTime) < time)
            localIterations = localIterations + 1;
            improved = false;
            
            % Get all nodes not in the current solution
            notSelected = setdiff(1:numNodes, currentNodes);
            
            bestNeighborScore = currentScore;
            bestSwap = [0, 0]; % [index in currentNodes, node from notSelected]
            validNeighborFound = false;
            
            % Evaluate all possible swaps (exhaustive neighborhood search)
            for i = 1:n
                for j = 1:length(notSelected)
                    % Create neighbor by swapping one node
                    neighborNodes = currentNodes;
                    neighborNodes(i) = notSelected(j);
                    
                    % Check if the neighbor satisfies Cmax constraint
                    [neighborScore, neighborMaxSP] = PerfSNS(G, neighborNodes);
                    
                    % Only consider neighbors that satisfy Cmax constraint
                    if neighborMaxSP <= Cmax && neighborScore < bestNeighborScore
                        bestNeighborScore = neighborScore;
                        bestSwap = [i, j];
                        improved = true;
                        validNeighborFound = true;
                    end
                end
            end
            
            % Apply the best swap if it improves the solution and respects constraints
            if validNeighborFound
                currentNodes(bestSwap(1)) = notSelected(bestSwap(2));
                currentScore = bestNeighborScore;
            end
        end
        
        % Update total iterations
        totalIterations = totalIterations + localIterations;
        
        % Update best solution if this run was better
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