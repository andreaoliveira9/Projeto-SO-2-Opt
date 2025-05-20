function [bestScore, bestNodes, generations, bestFoundTime] = GA_SNS(G, time, n, populationSize, mutationProb, elitistParam, Cmax)
    % GA_SNS Genetic Algorithm for Server Node Selection problem
    % Input:
    %   G - graph of the network
    %   time - time to run the method (in seconds)
    %   n - number of nodes to include in a solution
    %   populationSize - size of the population
    %   mutationProb - probability of mutation (0 to 1)
    %   elitistParam - maximum number of parents that can survive to next generation
    %   Cmax - maximum allowed shortest path length between any pair of selected nodes
    % Output:
    %   bestScore - best average shortest path length
    %   bestNodes - list of selected nodes
    %   generations - number of generations created
    %   bestFoundTime - time when the best solution was found
    
    numNodes = numnodes(G);
    bestScore = Inf;
    bestNodes = [];
    bestFoundTime = 0;
    startTime = tic;
    
    population = cell(populationSize, 1);
    fitness = zeros(populationSize, 1);
    feasible = false(populationSize, 1);
    
    % Initialize population
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
    
    generations = 1;
    
    while toc(startTime) < time
        offspringPopulation = cell(populationSize, 1);
        offspringFitness = zeros(populationSize, 1);
        offspringFeasible = false(populationSize, 1);
    
        for i = 1:populationSize
            % Tournament selection favors feasible solutions
            parent1Idx = tournamentSelection(fitness, feasible);
            parent2Idx = tournamentSelection(fitness, feasible);
            while parent2Idx == parent1Idx
                parent2Idx = tournamentSelection(fitness, feasible);
            end
    
            child = crossover(population{parent1Idx}, population{parent2Idx}, n);
    
            if rand() < mutationProb
                child = mutate(child, numNodes);
            end
    
            offspringPopulation{i} = child;
            [avgSP, maxSP] = PerfSNS(G, child);
            
            % Check if solution is feasible (maxSP <= Cmax)
            if maxSP <= Cmax
                offspringFeasible(i) = true;
                offspringFitness(i) = avgSP;
                
                if avgSP < bestScore
                    bestScore = avgSP;
                    bestNodes = child;
                    bestFoundTime = toc(startTime);
                end
            else
                % If infeasible, penalize the fitness
                offspringFitness(i) = avgSP + (maxSP - Cmax); % Penalty proportional to constraint violation
            end
        end
    
        % Elitism: Keep top elitistParam solutions from previous generation
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
    
        % Combine parents and offspring
        selectionPool = [filteredParents; offspringPopulation];
        selectionPoolFitness = [filteredParentFitness; offspringFitness];
        selectionPoolFeasible = [filteredParentFeasible; offspringFeasible];
        
        % Prioritize feasible solutions in sorting
        numSolutions = length(selectionPoolFitness);
        sortIndices = 1:numSolutions;
        
        % Custom sorting: First by feasibility, then by fitness
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
        
        % Take the best solutions for the next generation
        newPopulation = cell(populationSize, 1);
        newFitness = zeros(populationSize, 1);
        newFeasible = false(populationSize, 1);
    
        for i = 1:populationSize
            newPopulation{i} = selectionPool{i};
            newFitness(i) = selectionPoolFitness(i);
            newFeasible(i) = selectionPoolFeasible(i);
        end
    
        population = newPopulation;
        fitness = newFitness;
        feasible = newFeasible;
        generations = generations + 1;
    end
end

% Tournament selection helper function - modified to favor feasible solutions
function idx = tournamentSelection(fitness, feasible)
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

% Crossover helper function
function child = crossover(parent1, parent2, n)
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

% Mutation helper function
function mutated = mutate(solution, numNodes)
    mutated = solution;
    notSelected = setdiff(1:numNodes, solution);

    if ~isempty(notSelected)
        pos = randi(length(solution));
        newNode = notSelected(randi(length(notSelected)));
        mutated(pos) = newNode;
    end
end