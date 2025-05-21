
function [bestScore, bestNodes, totalIterations, bestFoundTime, totalEvaluated, totalValid, discardedCount] = GRASP_SNS(G, time, n, r, Cmax)
    % GRASP_SNS GRASP algorithm for Server Node Selection with Cmax constraint
    % Outputs:
    %   bestScore - best average shortest path length
    %   bestNodes - best solution found
    %   totalIterations - number of local improvement loops performed
    %   bestFoundTime - time when best solution was found
    %   totalEvaluated - total number of unique solutions evaluated
    %   totalValid - total number of valid solutions satisfying Cmax
    %   discardedCount - number of solutions discarded due to Cmax violation

    numNodes = numnodes(G);
    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;
    discardedCount = 0;
    evaluatedSolutions = containers.Map();
    totalEvaluated = 0;
    totalValid = 0;

    globalStartTime = tic;
    D = distances(G);

    while toc(globalStartTime) < time
        [currentNodes, discardedCountGreedy] = GreedyRandomizedConstruction(G, D, n, r, Cmax);
        discardedCount = discardedCount + discardedCountGreedy;

        if isempty(currentNodes) || length(unique(currentNodes)) < n
            continue;
        end

        key = mat2str(sort(currentNodes));
        if isKey(evaluatedSolutions, key)
            continue;
        end

        totalEvaluated = totalEvaluated + 1;

        [currentScore, maxSP] = PerfSNS(G, currentNodes);
        if maxSP > Cmax
            discardedCount = discardedCount + 1;
            continue;
        end

        totalValid = totalValid + 1;
        evaluatedSolutions(key) = currentScore;

        improved = true;
        localIterations = 0;

        while improved && (toc(globalStartTime) < time)
            localIterations = localIterations + 1;
            improved = false;
            notSelected = setdiff(1:numNodes, currentNodes);
            bestNeighborScore = currentScore;
            bestSwap = [0, 0];

            for i = 1:n
                for j = 1:length(notSelected)
                    neighborNodes = currentNodes;
                    neighborNodes(i) = notSelected(j);
                    neighborKey = mat2str(sort(neighborNodes));
                    if isKey(evaluatedSolutions, neighborKey)
                        continue;
                    end

                    totalEvaluated = totalEvaluated + 1;
                    [neighborScore, neighborMaxSP] = PerfSNS(G, neighborNodes);
                    if neighborMaxSP > Cmax
                        discardedCount = discardedCount + 1;
                        continue;
                    end

                    totalValid = totalValid + 1;
                    evaluatedSolutions(neighborKey) = neighborScore;

                    if neighborScore < bestNeighborScore
                        bestNeighborScore = neighborScore;
                        bestSwap = [i, notSelected(j)];
                        improved = true;
                    end
                end
            end

            if improved
                currentNodes(bestSwap(1)) = bestSwap(2);
                currentScore = bestNeighborScore;
                evaluatedSolutions(mat2str(sort(currentNodes))) = currentScore;
            end
        end

        totalIterations = totalIterations + localIterations;

        if currentScore < bestScore
            bestScore = currentScore;
            bestNodes = currentNodes;
            bestFoundTime = toc(globalStartTime);
        end
    end
end


function [nodes, discardedCount] = GreedyRandomizedConstruction(G, D, n, r, Cmax)
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
    discardedCount = 0;
    
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
                    discardedCount = discardedCount + 1;
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