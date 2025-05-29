function [bestScore, bestNodes, totalIterations, bestFoundTime] = GRASP_SNS(G, time, n, r, Cmax, seed)
% GRASP_SNS - Algoritmo GRASP para seleção de nós servidor com restrição de distância máxima
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
%   bestScore      - melhor valor da função objetivo encontrado
%   bestNodes      - melhor solução encontrada (array com os nós selecionados)
%   totalIterations - número total de iterações de busca local realizadas
%   bestFoundTime  - tempo em que a melhor solução foi encontrada
    
    % Configuração inicial do gerador de números aleatórios
    if nargin >= 6 && ~isempty(seed)
        rng(seed);
    end

    % Inicialização das variáveis de retorno
    bestScore = Inf;
    bestNodes = [];
    totalIterations = 0;
    bestFoundTime = 0;

    % Controle de tempo e cálculo prévio das distâncias
    globalStartTime = tic;
    D = distances(G); % Matriz de distâncias entre todos os pares de nós

    % Loop principal do GRASP: executa até esgotar o tempo limite
    while true
        elapsed = toc(globalStartTime);
        if elapsed >= time
            break;
        end
    
        % Verificação de tempo restante suficiente
        timeLeft = time - elapsed;
        if timeLeft < 0.01
            break;
        end
    
        % FASE 1: Construção gulosa randomizada
        currentNodes = GreedyRandomizedConstruction(G, D, n, r, Cmax);
        if toc(globalStartTime) >= time
            break;
        end

        % Validação da solução construída
        if isempty(currentNodes) || length(unique(currentNodes)) < n
            continue;
        end

        % Avaliação da solução e verificação da restrição Cmax
        [currentScore, maxSP] = PerfSNS(G, currentNodes);
        if maxSP > Cmax
            continue;
        end
        
        % Verificação de tempo antes da busca local
        if toc(globalStartTime) >= time
            break;
        end

        % FASE 2: Busca local com Hill Climbing
        [currentNodes, currentScore, localIterations, ~] = ...
            LocalSearch_SA_HC(G, currentNodes, currentScore, Cmax, globalStartTime, time);

        totalIterations = totalIterations + localIterations;

        % Atualização da melhor solução encontrada
        if currentScore < bestScore
            bestScore = currentScore;
            bestNodes = currentNodes;
            bestFoundTime = toc(globalStartTime);
        end
    end
end


function nodes = GreedyRandomizedConstruction(G, D, n, r, Cmax)
% GreedyRandomizedConstruction - Construção gulosa randomizada para seleção de nós servidor
%
% INPUTS:
%   G     - grafo da rede
%   D     - matriz de distâncias pré-calculada entre todos os pares de nós
%   n     - número de nós a selecionar
%   r     - tamanho da Lista Restrita de Candidatos (RCL)
%   Cmax  - distância máxima permitida entre servidores
%
% OUTPUTS:
%   nodes - array com os nós selecionados (vazio se solução inviável)
    
    numNodes = numnodes(G);
    nodes = zeros(1, n);
    remaining = 1:numNodes; % Nós ainda não selecionados
    
    % Cálculo da centralidade de cada nó (soma das distâncias a todos os outros)
    centrality = zeros(1, numNodes);
    for i = 1:numNodes
        centrality(i) = sum(D(i, :));
    end
    
    % Seleção do primeiro nó: escolha randomizada entre os mais centrais
    [~, sortedIndices] = sort(centrality);
    rclSize = min(r, numNodes);
    rcl = sortedIndices(1:rclSize); % Lista restrita de candidatos
    selectedIdx = randi(rclSize);
    nodes(1) = rcl(selectedIdx);
    
    % Remoção do nó selecionado da lista de candidatos
    remaining(remaining == nodes(1)) = [];
    
    invalidSolution = false;
    
    % Seleção iterativa dos nós restantes
    for k = 2:n
        % Avaliação do benefício de cada nó candidato
        benefit = zeros(1, length(remaining));
        validCandidates = false(1, length(remaining));
        
        for i = 1:length(remaining)
            candidate = remaining(i);
            tempSolution = [nodes(1:k-1), candidate];
            
            % Verificação da restrição Cmax: distância máxima entre servidores
            maxDistBetweenServers = 0;
            for a = 1:k-1
                serverA = nodes(a);
                dist = D(serverA, candidate);
                maxDistBetweenServers = max(maxDistBetweenServers, dist);
            end
            
            % Exclusão de candidatos que violam a restrição
            if maxDistBetweenServers > Cmax
                continue;
            end
            
            validCandidates(i) = true;
            
            % Cálculo da função objetivo: soma das distâncias mínimas
            totalDistance = 0;
            for node = 1:numNodes
                if ~ismember(node, tempSolution)
                    dist = min(D(node, tempSolution));
                    totalDistance = totalDistance + dist;
                end
            end
            
            % Benefício negativo para minimização
            benefit(i) = -totalDistance;
        end
        
        % Verificação de existência de candidatos válidos
        if ~any(validCandidates)
            invalidSolution = true;
            break;
        end
        
        % Construção da lista restrita de candidatos (RCL)
        validIndices = find(validCandidates);
        validBenefits = benefit(validCandidates);
        [~, sortedIndices] = sort(validBenefits, 'descend');
        
        rclSize = min(r, length(validIndices));
        rcl = remaining(validIndices(sortedIndices(1:rclSize)));
        
        % Seleção aleatória dentro da RCL
        selectedIdx = randi(rclSize);
        selectedNode = rcl(selectedIdx);
        
        nodes(k) = selectedNode;
        remaining = remaining(remaining ~= selectedNode);
    end
    
    % Retorno de solução vazia se infeasível
    if invalidSolution
        nodes = [];
    end
end

function [bestNodes, bestScore, localIterations, improved] = ...
    LocalSearch_SA_HC(G, currentNodes, currentScore, Cmax, globalStartTime, maxTime)
% LocalSearch_SA_HC - Busca local usando Hill Climbing com movimentos de swap
%
% INPUTS:
%   G               - grafo da rede
%   currentNodes    - solução inicial (array com nós selecionados)
%   currentScore    - valor da função objetivo da solução inicial
%   Cmax           - distância máxima permitida entre servidores
%   globalStartTime - tempo de início da execução global
%   maxTime        - tempo máximo total de execução
%
% OUTPUTS:
%   bestNodes       - melhor solução encontrada após busca local
%   bestScore       - valor da função objetivo da melhor solução
%   localIterations - número de iterações de busca local realizadas
%   improved        - flag indicando se houve melhoria na última iteração

    numNodes = numnodes(G);
    localIterations = 0;
    improved = true;
    bestNodes = currentNodes;
    bestScore = currentScore;

    % Loop de melhoria: continua enquanto houver melhorias
    while improved && toc(globalStartTime) < maxTime
        localIterations = localIterations + 1;
        improved = false;
        notSelected = setdiff(1:numNodes, bestNodes); % Nós não selecionados
        bestNeighborScore = bestScore;
        bestSwap = [0, 0];

        % Exploração da vizinhança: swap de cada nó selecionado
        for i = 1:length(bestNodes)
            for j = 1:length(notSelected)
                if toc(globalStartTime) >= maxTime
                    return;
                end

                % Geração de solução vizinha por swap
                neighborNodes = bestNodes;
                neighborNodes(i) = notSelected(j);

                % Avaliação da solução vizinha
                [neighborScore, neighborMaxSP] = PerfSNS(G, neighborNodes);

                % Verificação da restrição Cmax
                if neighborMaxSP > Cmax
                    continue;
                end

                % Identificação de melhoria (critério Hill Climbing)
                if neighborScore < bestNeighborScore
                    bestNeighborScore = neighborScore;
                    bestSwap = [i, notSelected(j)];
                    improved = true;
                end
            end
        end

        % Aplicação da melhor melhoria encontrada
        if improved
            bestNodes(bestSwap(1)) = bestSwap(2);
            bestScore = bestNeighborScore;
        end
    end
end