function [bestScore, bestNodes, generations, bestFoundTime] = GA_SNS(G, time, n, populationSize, mutationProb, elitistParam, Cmax, seed)
    % GA_SNS Genetic Algorithm for Server Node Selection problem
    % Input:
    %   G              - graph of the network
    %   time           - time to run the method (in seconds)
    %   n              - number of nodes to include in a solution
    %   populationSize - size of the population
    %   mutationProb   - probability of mutation (0 to 1)
    %   elitistParam   - maximum number of parents that can survive to next generation
    %   Cmax           - maximum allowed shortest path length between any pair of selected nodes
    %   seed           - semente para o gerador de números aleatórios (opcional)
    % Output:
    %   bestScore     - best average shortest path length
    %   bestNodes     - list of selected nodes
    %   generations   - number of generations created
    %   bestFoundTime - time when the best solution was found
    
    % Configuração inicial do gerador de números aleatórios
    if nargin >= 8 && ~isempty(seed)
        rng(seed);
    end

    numNodes = numnodes(G);
    bestScore = Inf;
    bestNodes = [];
    bestFoundTime = 0;
    startTime = tic;
    
    % Initialize population and evaluate initial solutions
    [population, fitness, feasible, bestScore, bestNodes, bestFoundTime] = ...
        initializePopulation(G, populationSize, numNodes, n, Cmax, bestScore, bestNodes, bestFoundTime, startTime);
    
    generations = 1;
    
    % Main evolutionary loop
    while true
        elapsed = toc(startTime);
        if elapsed >= time
            break;
        end
    
        timeLeft = time - elapsed;
        if timeLeft < 0.01  % margem de segurança mínima para uma geração
            break;
        end

        % Generate offspring population
        [offspringPopulation, offspringFitness, offspringFeasible, bestScore, bestNodes, bestFoundTime] = ...
            generateOffspring(G, population, fitness, feasible, populationSize, mutationProb, ...
            numNodes, n, Cmax, bestScore, bestNodes, bestFoundTime, startTime, time);

        if toc(startTime) >= time
            break;
        end
    
        % Apply elitism and survival selection
        [population, fitness, feasible] = applySurvivalSelection(population, fitness, feasible, ...
            offspringPopulation, offspringFitness, offspringFeasible, populationSize, elitistParam);
        
        generations = generations + 1;
    end
end

function [population, fitness, feasible, bestScore, bestNodes, bestFoundTime] = ...
    initializePopulation(G, populationSize, numNodes, n, Cmax, bestScore, bestNodes, bestFoundTime, startTime)
    % initializePopulation - Initialize population and evaluate initial solutions
    % Input:
    %   G              - graph of the network
    %   populationSize - size of the population
    %   numNodes       - total number of nodes in the graph
    %   n              - number of nodes to include in a solution
    %   Cmax           - maximum allowed shortest path length between any pair of selected nodes
    %   bestScore      - current best score (will be updated)
    %   bestNodes      - current best solution (will be updated)
    %   bestFoundTime  - time when best solution was found (will be updated)
    %   startTime      - algorithm start time reference
    % Output:
    %   population     - cell array with initial population
    %   fitness        - fitness values for each individual
    %   feasible       - feasibility status for each individual
    %   bestScore      - updated best score
    %   bestNodes      - updated best solution
    %   bestFoundTime  - updated time when best solution was found
    
    population = cell(populationSize, 1);
    fitness = zeros(populationSize, 1);
    feasible = false(populationSize, 1);
    
    for i = 1:populationSize
        population{i} = randperm(numNodes, n);
        [avgSP, maxSP] = PerfSNS(G, population{i});
        
        % Check if solution is feasible (maxSP <= Cmax)
        if maxSP <= Cmax
            feasible(i) = true;
            fitness(i) = avgSP;
            
            if avgSP < bestScore
                bestScore = avgSP;
                bestNodes = population{i};
                bestFoundTime = toc(startTime);
            end
        else
            % If infeasible, penalize the fitness
            fitness(i) = avgSP + (maxSP - Cmax); % Penalty proportional to constraint violation
        end
    end
end

function [offspringPopulation, offspringFitness, offspringFeasible, bestScore, bestNodes, bestFoundTime] = ...
    generateOffspring(G, population, fitness, feasible, populationSize, mutationProb, ...
    numNodes, n, Cmax, bestScore, bestNodes, bestFoundTime, startTime, timeLimit)
    % generateOffspring - Generate offspring population through selection, crossover, and mutation
    % Input:
    %   G                  - graph of the network
    %   population         - current population
    %   fitness            - fitness values of current population
    %   feasible           - feasibility status of current population
    %   populationSize     - size of the population
    %   mutationProb       - probability of mutation (0 to 1)
    %   numNodes           - total number of nodes in the graph
    %   n                  - number of nodes to include in a solution
    %   Cmax               - maximum allowed shortest path length between any pair of selected nodes
    %   bestScore          - current best score (will be updated)
    %   bestNodes          - current best solution (will be updated)
    %   bestFoundTime      - time when best solution was found (will be updated)
    %   startTime          - algorithm start time reference
    %   timeLimit          - maximum time for algorithm execution
    % Output:
    %   offspringPopulation - cell array with offspring population
    %   offspringFitness    - fitness values for each offspring
    %   offspringFeasible   - feasibility status for each offspring
    %   bestScore           - updated best score
    %   bestNodes           - updated best solution
    %   bestFoundTime       - updated time when best solution was found
    
    offspringPopulation = cell(populationSize, 1);
    offspringFitness = zeros(populationSize, 1);
    offspringFeasible = false(populationSize, 1);

    for i = 1:populationSize
        if toc(startTime) >= timeLimit
            break;
        end
        
        % Parent selection
        [parent1Idx, parent2Idx] = selectParents(fitness, feasible);
        
        % Crossover
        child = crossover(population{parent1Idx}, population{parent2Idx}, n);

        % Mutation
        if rand() < mutationProb
            child = mutate(child, numNodes);
        end

        offspringPopulation{i} = child;
        
        if toc(startTime) >= timeLimit
            break;
        end
        
        % Evaluate offspring
        [offspringFitness(i), offspringFeasible(i), bestScore, bestNodes, bestFoundTime] = ...
            evaluateAndUpdateBest(G, child, Cmax, bestScore, bestNodes, bestFoundTime, startTime);
    end
end

function [parent1Idx, parent2Idx] = selectParents(fitness, feasible)
    % selectParents - Select two different parents using tournament selection
    % Input:
    %   fitness    - fitness values of current population
    %   feasible   - feasibility status of current population
    % Output:
    %   parent1Idx - index of first selected parent
    %   parent2Idx - index of second selected parent
    
    parent1Idx = tournamentSelection(fitness, feasible);
    parent2Idx = tournamentSelection(fitness, feasible);
    
    while parent2Idx == parent1Idx
        parent2Idx = tournamentSelection(fitness, feasible);
    end
end

function [solutionFitness, solutionFeasible, bestScore, bestNodes, bestFoundTime] = ...
    evaluateAndUpdateBest(G, solution, Cmax, bestScore, bestNodes, bestFoundTime, startTime)
    % evaluateAndUpdateBest - Evaluate a solution and update best solution if necessary
    % Input:
    %   G             - graph of the network
    %   solution      - solution to be evaluated
    %   Cmax          - maximum allowed shortest path length between any pair of selected nodes
    %   bestScore     - current best score (will be updated)
    %   bestNodes     - current best solution (will be updated)
    %   bestFoundTime - time when best solution was found (will be updated)
    %   startTime     - algorithm start time reference
    % Output:
    %   solutionFitness  - fitness value of the evaluated solution
    %   solutionFeasible - feasibility status of the evaluated solution
    %   bestScore        - updated best score
    %   bestNodes        - updated best solution
    %   bestFoundTime    - updated time when best solution was found
    
    [avgSP, maxSP] = PerfSNS(G, solution);
    
    % Check if solution is feasible (maxSP <= Cmax)
    if maxSP <= Cmax
        solutionFeasible = true;
        solutionFitness = avgSP;
        
        if avgSP < bestScore
            bestScore = avgSP;
            bestNodes = solution;
            bestFoundTime = toc(startTime);
        end
    else
        % If infeasible, penalize the fitness
        solutionFeasible = false;
        solutionFitness = avgSP + (maxSP - Cmax); % Penalty proportional to constraint violation
    end
end

function [newPopulation, newFitness, newFeasible] = applySurvivalSelection(population, fitness, feasible, ...
    offspringPopulation, offspringFitness, offspringFeasible, populationSize, elitistParam)
    % applySurvivalSelection - Apply elitism and survival selection to form next generation
    % Input:
    %   population          - current population
    %   fitness             - fitness values of current population
    %   feasible            - feasibility status of current population
    %   offspringPopulation - offspring population
    %   offspringFitness    - fitness values of offspring population
    %   offspringFeasible   - feasibility status of offspring population
    %   populationSize      - size of the population
    %   elitistParam        - maximum number of parents that can survive to next generation
    % Output:
    %   newPopulation - selected population for next generation
    %   newFitness    - fitness values for next generation
    %   newFeasible   - feasibility status for next generation
    
    % Apply elitism: Keep top elitistParam solutions from previous generation
    [filteredParents, filteredParentFitness, filteredParentFeasible] = ...
        applyElitism(population, fitness, feasible, elitistParam, populationSize);
    
    % Combine parents and offspring
    selectionPool = [filteredParents; offspringPopulation];
    selectionPoolFitness = [filteredParentFitness; offspringFitness];
    selectionPoolFeasible = [filteredParentFeasible; offspringFeasible];
    
    % Sort by feasibility first, then by fitness
    [selectionPool, selectionPoolFitness, selectionPoolFeasible] = ...
        sortSolutionsByFeasibilityAndFitness(selectionPool, selectionPoolFitness, selectionPoolFeasible);
    
    % Select best solutions for next generation
    [newPopulation, newFitness, newFeasible] = ...
        selectSurvivors(selectionPool, selectionPoolFitness, selectionPoolFeasible, populationSize);
end

function [filteredParents, filteredParentFitness, filteredParentFeasible] = ...
    applyElitism(population, fitness, feasible, elitistParam, populationSize)
    % applyElitism - Apply elitism by keeping the best solutions from current generation
    % Input:
    %   population     - current population
    %   fitness        - fitness values of current population
    %   feasible       - feasibility status of current population
    %   elitistParam   - maximum number of parents that can survive to next generation
    %   populationSize - size of the population
    % Output:
    %   filteredParents        - elite solutions from current generation
    %   filteredParentFitness  - fitness values of elite solutions
    %   filteredParentFeasible - feasibility status of elite solutions
    
    [parentFitness, parentIndices] = sort(fitness);
    filteredParentCount = min(elitistParam, populationSize);
    filteredParents = cell(filteredParentCount, 1);
    filteredParentFitness = zeros(filteredParentCount, 1);
    filteredParentFeasible = false(filteredParentCount, 1);

    for i = 1:filteredParentCount
        filteredParents{i} = population{parentIndices(i)};
        filteredParentFitness(i) = parentFitness(i);
        filteredParentFeasible(i) = feasible(parentIndices(i));
    end
end

function [sortedPool, sortedFitness, sortedFeasible] = ...
    sortSolutionsByFeasibilityAndFitness(selectionPool, selectionPoolFitness, selectionPoolFeasible)
    % sortSolutionsByFeasibilityAndFitness - Custom sorting: First by feasibility, then by fitness
    % Input:
    %   selectionPool         - combined pool of solutions
    %   selectionPoolFitness  - fitness values of solutions in pool
    %   selectionPoolFeasible - feasibility status of solutions in pool
    % Output:
    %   sortedPool     - solutions sorted by feasibility then fitness
    %   sortedFitness  - fitness values in sorted order
    %   sortedFeasible - feasibility status in sorted order
    
    numSolutions = length(selectionPoolFitness);
    sortIndices = 1:numSolutions;
    
    for i = 1:numSolutions
        for j = i+1:numSolutions
            if selectionPoolFeasible(i) < selectionPoolFeasible(j) || ...
               (selectionPoolFeasible(i) == selectionPoolFeasible(j) && selectionPoolFitness(i) > selectionPoolFitness(j))
                % Swap
                tempFitness = selectionPoolFitness(i);
                tempFeasible = selectionPoolFeasible(i);
                tempSol = selectionPool{i};
                tempIdx = sortIndices(i);
                
                selectionPoolFitness(i) = selectionPoolFitness(j);
                selectionPoolFeasible(i) = selectionPoolFeasible(j);
                selectionPool{i} = selectionPool{j};
                sortIndices(i) = sortIndices(j);
                
                selectionPoolFitness(j) = tempFitness;
                selectionPoolFeasible(j) = tempFeasible;
                selectionPool{j} = tempSol;
                sortIndices(j) = tempIdx;
            end
        end
    end
    
    sortedPool = selectionPool;
    sortedFitness = selectionPoolFitness;
    sortedFeasible = selectionPoolFeasible;
end

function [newPopulation, newFitness, newFeasible] = ...
    selectSurvivors(selectionPool, selectionPoolFitness, selectionPoolFeasible, populationSize)
    % selectSurvivors - Select the best solutions for the next generation
    % Input:
    %   selectionPool         - pool of candidate solutions
    %   selectionPoolFitness  - fitness values of candidate solutions
    %   selectionPoolFeasible - feasibility status of candidate solutions
    %   populationSize        - size of the population
    % Output:
    %   newPopulation - selected solutions for next generation
    %   newFitness    - fitness values of selected solutions
    %   newFeasible   - feasibility status of selected solutions
    
    newPopulation = cell(populationSize, 1);
    newFitness = zeros(populationSize, 1);
    newFeasible = false(populationSize, 1);

    for i = 1:populationSize
        newPopulation{i} = selectionPool{i};
        newFitness(i) = selectionPoolFitness(i);
        newFeasible(i) = selectionPoolFeasible(i);
    end
end

function idx = tournamentSelection(fitness, feasible)
    % tournamentSelection - Tournament selection favoring feasible solutions
    % Input:
    %   fitness  - fitness values of current population
    %   feasible - feasibility status of current population
    % Output:
    %   idx - index of selected individual
    tournamentSize = 2;
    candidates = randi(length(fitness), 1, tournamentSize);
    
    % Check if any candidates are feasible
    feasibleCandidates = candidates(feasible(candidates));
    
    if ~isempty(feasibleCandidates)
        % If there are feasible candidates, select the best among them
        [~, bestIdx] = min(fitness(feasibleCandidates));
        idx = feasibleCandidates(bestIdx);
    else
        % If no feasible candidates, select based on lowest fitness value
        [~, bestIdx] = min(fitness(candidates));
        idx = candidates(bestIdx);
    end
end

function child = crossover(parent1, parent2, n)
    % crossover - Crossover operation between two parents
    % Input:
    %   parent1 - first parent solution
    %   parent2 - second parent solution
    %   n       - number of nodes in a solution
    % Output:
    %   child - offspring solution
    selected = parent1;
    uniqueToParent2 = setdiff(parent2, parent1);

    if ~isempty(uniqueToParent2)
        numToSwap = randi([1, floor(n/2)]);
        numToSwap = min(numToSwap, length(uniqueToParent2));
        posToReplace = randperm(n, numToSwap);
        elemsToInsert = uniqueToParent2(randperm(length(uniqueToParent2), numToSwap));
        selected(posToReplace) = elemsToInsert;
    end

    child = selected;
end

function mutated = mutate(solution, numNodes)
    % mutate - Mutation operation on a solution
    % Input:
    %   solution - solution to be mutated
    %   numNodes - total number of nodes in the graph
    % Output:
    %   mutated - mutated solution
    mutated = solution;
    notSelected = setdiff(1:numNodes, solution);

    if ~isempty(notSelected)
        pos = randi(length(solution));
        newNode = notSelected(randi(length(notSelected)));
        mutated(pos) = newNode;
    end
end