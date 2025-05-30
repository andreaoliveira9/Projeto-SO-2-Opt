function [bestScore, bestNodes, totalIterations, bestFoundTime] = GRASP_SNS_Optimized(G, time, n, r, Cmax, seed)
% GRASP_SNS_Optimized - Versão otimizada do algoritmo GRASP para seleção de nós servidor
%
% INPUTS:
%   G     - grafo da rede (objeto graph do MATLAB)
%   time  - tempo limite de execução em segundos
%   n     - número de nós servidor a selecionar
%   r     - tamanho da Lista Restrita de Candidatos (RCL)
%   Cmax  - distância máxima permitida entre qualquer par de servidores
%   seed  - semente para o gerador de números aleatórios (opcional)
%
% OUTPUTS:
%   bestScore       - melhor valor da função objetivo encontrado
%   bestNodes       - melhor solução encontrada (array com os nós selecionados)
%   totalIterations - número total de iterações de busca local realizadas
%   bestFoundTime   - tempo em que a melhor solução foi encontrada
%
% OTIMIZAÇÕES IMPLEMENTADAS:
%   1. Pré-filtragem de candidatos baseada em Cmax
%   2. Cache de distâncias e estruturas auxiliares
%   3. Avaliação incremental de soluções
%   4. Busca local mais eficiente com early stopping

    % Configuração inicial do gerador de números aleatórios
    if nargin >= 6 && ~isempty(seed)
        rng(seed);
    end

    % Inicialização das variáveis de retorno
    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;
    
    globalStartTime = tic;
    
    % Pré-computação única da matriz de distâncias
    D = distances(G);
    numNodes = numnodes(G);
    
    % OTIMIZAÇÃO 1: Pré-filtragem baseada na restrição Cmax
    % Identifica pares de nós que podem coexistir numa solução válida
    validPairs = D <= Cmax;
    
    % Conversão para lista de adjacências (mais eficiente que matriz)
    validNeighbors = cell(numNodes, 1);
    for i = 1:numNodes
        validNeighbors{i} = find(validPairs(i, :));
    end
    
    % Pré-computação das centralidades ordenadas por importância
    centrality = 1 ./ sum(D, 2)';
    [~, centralityOrder] = sort(centrality, 'descend');
    
    % OTIMIZAÇÃO 2: Cache para evitar recálculos de soluções já avaliadas
    scoreCache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    % Loop principal do GRASP otimizado
    while true
        elapsed = toc(globalStartTime);
        if elapsed >= time
            break;
        end
        
        timeLeft = time - elapsed;
        if timeLeft < 0.01
            break;
        end
        
        % FASE 1: Construção gulosa randomizada com pré-filtragem
        currentNodes = GreedyRandomizedConstruction_Optimized(D, validNeighbors, ...
            centralityOrder, n, r, numNodes);
            
        if toc(globalStartTime) >= time
            break;
        end
        
        % Validação da solução construída
        if isempty(currentNodes) || length(unique(currentNodes)) < n
            continue;
        end
        
        % OTIMIZAÇÃO 3: Avaliação com cache para evitar recálculos
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
        
        % FASE 2: Busca local otimizada com early stopping
        [currentNodes, currentScore, localIterations] = ...
            LocalSearch_SA_HC_Optimized(G, D, currentNodes, currentScore, Cmax, ...
            validNeighbors, scoreCache, globalStartTime, time);
        
        totalIterations = totalIterations + localIterations;
        
        % Atualização da melhor solução global
        if currentScore < bestScore
            bestScore = currentScore;
            bestNodes = currentNodes;
            bestFoundTime = toc(globalStartTime);
        end
    end
end

function nodes = GreedyRandomizedConstruction_Optimized(D, validNeighbors, ...
    centralityOrder, n, r, numNodes)
% GreedyRandomizedConstruction_Optimized- Construção gulosa com pré-filtragem otimizada
%
% INPUTS:
%   D               - matriz de distâncias pré-calculada
%   validNeighbors  - lista de adjacências de candidatos válidos por Cmax
%   centralityOrder - nós ordenados por centralidade decrescente
%   n               - número de nós a selecionar
%   r               - tamanho da Lista Restrita de Candidatos (RCL)
%   numNodes        - número total de nós no grafo
%
% OUTPUTS:
%   nodes - array com os nós selecionados (vazio se solução infeasível)
    
    nodes = zeros(1, n);
    usedNodes = false(1, numNodes); % Máscara booleana para nós já selecionados
    
    % Seleção do primeiro nó: escolha aleatória entre os mais centrais
    rclSize = min(r, numNodes);
    firstNodeIdx = randi(rclSize);
    nodes(1) = centralityOrder(firstNodeIdx);
    usedNodes(nodes(1)) = true;
    
    % Inicialização do conjunto de candidatos válidos
    validCandidates = validNeighbors{nodes(1)};
    validCandidates = validCandidates(~usedNodes(validCandidates));
    
    % Construção iterativa mantendo compatibilidade com Cmax
    for k = 2:n
        if isempty(validCandidates)
            nodes = []; % Retorna solução vazia se não há candidatos válidos
            return;
        end
        
        % Avaliação do benefício apenas dos candidatos pré-filtrados
        benefits = zeros(1, length(validCandidates));
        
        for i = 1:length(validCandidates)
            candidate = validCandidates(i);
            tempSolution = [nodes(1:k-1), candidate];
            
            % Cálculo incremental da função objetivo
            totalDistance = 0;
            for node = 1:numNodes
                if ~usedNodes(node) && node ~= candidate
                    minDist = min(D(node, tempSolution));
                    totalDistance = totalDistance + minDist;
                end
            end
            
            benefits(i) = -totalDistance; % Negativo para minimização
        end
        
        % Construção da RCL com os melhores candidatos
        [~, sortedIdx] = sort(benefits, 'descend');
        rclSize = min(r, length(validCandidates));
        rcl = validCandidates(sortedIdx(1:rclSize));
        
        % Seleção aleatória dentro da RCL
        selectedNode = rcl(randi(rclSize));
        nodes(k) = selectedNode;
        usedNodes(selectedNode) = true;
        
        % Atualização dos candidatos válidos por interseção
        newValidCandidates = validNeighbors{selectedNode};
        validCandidates = intersect(validCandidates, newValidCandidates);
        validCandidates = validCandidates(~usedNodes(validCandidates));
    end
end

function [bestNodes, bestScore, localIterations] = ...
    LocalSearch_SA_HC_Optimized(G, D, currentNodes, currentScore, Cmax, ...
    validNeighbors, scoreCache, globalStartTime, maxTime)
% LocalSearch_SA_HC_Optimized - Busca local otimizada com cache e early stopping
%
% INPUTS:
%   G               - grafo da rede
%   D               - matriz de distâncias pré-calculada
%   currentNodes    - solução inicial (array com nós selecionados)
%   currentScore    - valor da função objetivo da solução inicial
%   Cmax            - distância máxima permitida entre servidores
%   validNeighbors  - lista de adjacências de candidatos válidos
%   scoreCache      - cache de avaliações já computadas
%   globalStartTime - tempo de início da execução global
%   maxTime         - tempo máximo total de execução
%
% OUTPUTS:
%   bestNodes       - melhor solução encontrada após busca local
%   bestScore       - valor da função objetivo da melhor solução
%   localIterations - número de iterações de busca local realizadas
    
    localIterations = 0;
    improved = true;
    bestNodes = currentNodes;
    bestScore = currentScore;
    numNodes = size(D, 1);
    
    % Pré-computação do conjunto de nós não selecionados
    usedNodes = false(1, numNodes);
    usedNodes(bestNodes) = true;
    notSelected = find(~usedNodes);
    
    % Loop de melhoria com early stopping
    while improved && toc(globalStartTime) < maxTime
        localIterations = localIterations + 1;
        improved = false;
        bestNeighborScore = bestScore;
        bestSwap = [0, 0];
        
        % Pré-filtragem de movimentos válidos baseada em Cmax
        moves = [];
        for i = 1:length(bestNodes)
            currentNode = bestNodes(i);
            % Candidatos que respeitam Cmax com o nó atual
            validCandidatesForSwap = validNeighbors{currentNode};
            validCandidatesForSwap = intersect(validCandidatesForSwap, notSelected);
            
            % Verificação de compatibilidade com demais nós da solução
            otherNodes = bestNodes(bestNodes ~= currentNode);
            for candidate = validCandidatesForSwap
                if all(D(candidate, otherNodes) <= Cmax)
                    moves = [moves; i, candidate];
                end
            end
        end
        
        % Randomização da ordem dos movimentos para diversificação
        moves = moves(randperm(size(moves, 1)), :);
        
        % Exploração da vizinhança com avaliação em cache
        for moveIdx = 1:size(moves, 1)
            if toc(globalStartTime) >= maxTime
                return;
            end
            
            i = moves(moveIdx, 1);
            candidate = moves(moveIdx, 2);
            
            % Geração da solução vizinha
            neighborNodes = bestNodes;
            neighborNodes(i) = candidate;
            
            % Avaliação com cache para evitar recálculos
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
            
            % Critério de melhoria com early stopping
            if neighborScore < bestNeighborScore
                bestNeighborScore = neighborScore;
                bestSwap = [i, candidate];
                improved = true;
                
                % OTIMIZAÇÃO: Para na primeira melhoria encontrada
                break;
            end
        end
        
        % Aplicação da melhor melhoria e atualização das estruturas
        if improved
            oldNode = bestNodes(bestSwap(1));
            newNode = bestSwap(2);
            
            bestNodes(bestSwap(1)) = newNode;
            bestScore = bestNeighborScore;
            
            % Atualização eficiente do conjunto de nós não selecionados
            notSelected = notSelected(notSelected ~= newNode);
            notSelected = [notSelected, oldNode];
        end
    end
end