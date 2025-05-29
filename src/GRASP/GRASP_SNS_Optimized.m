function [bestScore, bestNodes, totalIterations, bestFoundTime] = GRASP_SNS_Optimized(G, time, n, r, Cmax)
    % GRASP_SNS_Otimizado - Versão otimizada do GRASP_SNS
    % 
    % Principais otimizações:
    % 1. Pré-filtragem de candidatos baseada em Cmax
    % 2. Cache de distâncias e estruturas auxiliares
    % 3. Avaliação incremental de soluções
    % 4. Busca local mais eficiente com early stopping
    
    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;
    
    globalStartTime = tic;
    
    % Pré-computação de todas as distâncias (uma única vez)
    D = distances(G);
    numNodes = numnodes(G);
    
    % PRÉ-FILTRAGEM: Criar grafo de candidatos válidos baseado em Cmax
    % Dois nós podem estar na mesma solução se a distância entre eles <= Cmax
    validPairs = D <= Cmax;
    
    % Criar lista de adjacências para candidatos válidos (mais eficiente que matriz)
    validNeighbors = cell(numNodes, 1);
    for i = 1:numNodes
        validNeighbors{i} = find(validPairs(i, :));
    end
    
    % Pré-computar centralidades (inverso da soma das distâncias)
    centrality = 1 ./ sum(D, 2)';
    [~, centralityOrder] = sort(centrality, 'descend');
    
    % Cache para evitar recálculos na busca local
    scoreCache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    while true
        elapsed = toc(globalStartTime);
        if elapsed >= time
            break;
        end
        
        timeLeft = time - elapsed;
        if timeLeft < 0.01
            break;
        end
        
        % Construção gulosa randomizada otimizada
        currentNodes = GreedyRandomizedConstruction_Otimizada(D, validNeighbors, ...
            centralityOrder, n, r, numNodes);
            
        if toc(globalStartTime) >= time
            break;
        end
        
        if isempty(currentNodes) || length(unique(currentNodes)) < n
            continue;
        end
        
        % Avaliação rápida usando cache
        nodeKey = mat2str(sort(currentNodes));
        if isKey(scoreCache, nodeKey)
            currentScore = scoreCache(nodeKey);
        else
            [currentScore, maxSP] = PerfSNS(G, currentNodes);
            if maxSP > Cmax
                continue;
            end
            scoreCache(nodeKey) = currentScore;
        end
        
        if toc(globalStartTime) >= time
            break;
        end
        
        % Busca local otimizada
        [currentNodes, currentScore, localIterations] = ...
            LocalSearch_Otimizada(G, D, currentNodes, currentScore, Cmax, ...
            validNeighbors, scoreCache, globalStartTime, time);
        
        totalIterations = totalIterations + localIterations;
        
        if currentScore < bestScore
            bestScore = currentScore;
            bestNodes = currentNodes;
            bestFoundTime = toc(globalStartTime);
        end
    end
end

function nodes = GreedyRandomizedConstruction_Otimizada(D, validNeighbors, ...
    centralityOrder, n, r, numNodes)
    % Construção gulosa randomizada otimizada com pré-filtragem
    
    nodes = zeros(1, n);
    usedNodes = false(1, numNodes);
    
    % Selecionar primeiro nó da lista de centralidade
    rclSize = min(r, numNodes);
    firstNodeIdx = randi(rclSize);
    nodes(1) = centralityOrder(firstNodeIdx);
    usedNodes(nodes(1)) = true;
    
    % Manter conjunto de candidatos válidos (interseção das vizinhanças válidas)
    validCandidates = validNeighbors{nodes(1)};
    validCandidates = validCandidates(~usedNodes(validCandidates));
    
    for k = 2:n
        if isempty(validCandidates)
            nodes = []; % Solução inválida
            return;
        end
        
        % Calcular benefícios apenas para candidatos válidos
        benefits = zeros(1, length(validCandidates));
        
        for i = 1:length(validCandidates)
            candidate = validCandidates(i);
            tempSolution = [nodes(1:k-1), candidate];
            
            % Avaliação incremental mais eficiente
            totalDistance = 0;
            for node = 1:numNodes
                if ~usedNodes(node) && node ~= candidate
                    minDist = min(D(node, tempSolution));
                    totalDistance = totalDistance + minDist;
                end
            end
            
            benefits(i) = -totalDistance; % Negativo porque queremos minimizar
        end
        
        % Construir RCL com top-r candidatos
        [~, sortedIdx] = sort(benefits, 'descend');
        rclSize = min(r, length(validCandidates));
        rcl = validCandidates(sortedIdx(1:rclSize));
        
        % Seleção aleatória do RCL
        selectedNode = rcl(randi(rclSize));
        nodes(k) = selectedNode;
        usedNodes(selectedNode) = true;
        
        % Atualizar candidatos válidos (interseção)
        newValidCandidates = validNeighbors{selectedNode};
        validCandidates = intersect(validCandidates, newValidCandidates);
        validCandidates = validCandidates(~usedNodes(validCandidates));
    end
end

function [bestNodes, bestScore, localIterations] = ...
    LocalSearch_Otimizada(G, D, currentNodes, currentScore, Cmax, ...
    validNeighbors, scoreCache, globalStartTime, maxTime)
    % Busca local otimizada com avaliação incremental e cache
    
    localIterations = 0;
    improved = true;
    bestNodes = currentNodes;
    bestScore = currentScore;
    numNodes = size(D, 1);
    
    % Pré-computar conjunto de nós não selecionados
    usedNodes = false(1, numNodes);
    usedNodes(bestNodes) = true;
    notSelected = find(~usedNodes);
    
    while improved && toc(globalStartTime) < maxTime
        localIterations = localIterations + 1;
        improved = false;
        bestNeighborScore = bestScore;
        bestSwap = [0, 0];
        
        % Ordenar movimentos por potencial de melhoria
        moves = [];
        for i = 1:length(bestNodes)
            currentNode = bestNodes(i);
            % Apenas considerar nós válidos (que respeitam Cmax com os outros)
            validCandidatesForSwap = validNeighbors{currentNode};
            validCandidatesForSwap = intersect(validCandidatesForSwap, notSelected);
            
            % Verificar compatibilidade com outros nós da solução
            otherNodes = bestNodes(bestNodes ~= currentNode);
            for candidate = validCandidatesForSwap
                if all(D(candidate, otherNodes) <= Cmax)
                    moves = [moves; i, candidate];
                end
            end
        end
        
        % Embaralhar movimentos para diversificação
        moves = moves(randperm(size(moves, 1)), :);
        
        for moveIdx = 1:size(moves, 1)
            if toc(globalStartTime) >= maxTime
                return;
            end
            
            i = moves(moveIdx, 1);
            candidate = moves(moveIdx, 2);
            
            neighborNodes = bestNodes;
            neighborNodes(i) = candidate;
            
            % Usar cache quando possível
            nodeKey = mat2str(sort(neighborNodes));
            if isKey(scoreCache, nodeKey)
                neighborScore = scoreCache(nodeKey);
            else
                [neighborScore, neighborMaxSP] = PerfSNS(G, neighborNodes);
                if neighborMaxSP > Cmax
                    continue;
                end
                scoreCache(nodeKey) = neighborScore;
            end
            
            if neighborScore < bestNeighborScore
                bestNeighborScore = neighborScore;
                bestSwap = [i, candidate];
                improved = true;
                
                % Early stopping na primeira melhoria (first improvement)
                break;
            end
        end
        
        if improved
            % Atualizar solução e estruturas auxiliares
            oldNode = bestNodes(bestSwap(1));
            newNode = bestSwap(2);
            
            bestNodes(bestSwap(1)) = newNode;
            bestScore = bestNeighborScore;
            
            % Atualizar notSelected
            notSelected = notSelected(notSelected ~= newNode);
            notSelected = [notSelected, oldNode];
        end
    end
end