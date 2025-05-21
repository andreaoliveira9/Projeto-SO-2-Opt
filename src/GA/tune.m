
% tune_ga.m
% Testa diferentes combinações de parâmetros para o algoritmo GA e identifica os melhores

addpath('SupportingFiles');
clear; clc;

% --- Parâmetros do problema ---
n = 12;
time = 30;
Cmax = 1000;
numRuns = 10;
populationSizes = [20, 50, 100, 150];
mutationProbs = [0.05, 0.1, 0.2];
elitistParams = [1, 5, 10];

% --- Carregamento dos dados ---
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
G = graph(L);

% --- Resultados ---
results = [];

fprintf('Testando diferentes configurações de parâmetros GA...\n');

for p = populationSizes
    for m = mutationProbs
        for e = elitistParams
            scores = zeros(1, numRuns);
            times = zeros(1, numRuns);
            evaluations = zeros(1, numRuns);
            validCounts = zeros(1, numRuns);
            discards = zeros(1, numRuns);
            fprintf('\n--- Pop=%d Mut=%.2f Elite=%d ---\n', p, m, e);

            for i = 1:numRuns
                [score, nodes, ~, foundTime] = GA_SNS(G, time, n, p, m, e, Cmax);
                [avgSP, ~] = PerfSNS(G, nodes);
                scores(i) = avgSP;
                times(i) = foundTime;
                fprintf('Run %d: SP = %.4f | Tempo = %.2fs\n', i, avgSP, foundTime);
            end

            results = [results; p, m, e, ...
                min(scores), mean(scores), max(scores), ...
                mean(times)];

            fprintf('Resumo para Pop=%d Mut=%.2f Elite=%d -> Min: %.4f | Média: %.4f | Max: %.4f | Tempo médio: %.4f\n', ...
                p, m, e, min(scores), mean(scores), max(scores), mean(times));
        end
    end
end

% --- Tabela final ---
fprintf('\n==================== Resultados Finais ====================\n');
fprintf(' Pop | Mut | Elite |   Min SP   |   Média SP   |   Max SP   | Tempo Médio (s)\n');
fprintf('------------------------------------------------------------\n');
for i = 1:size(results,1)
    fprintf('%3d | %.2f | %5d | %10.4f | %12.4f | %9.4f | %14.2f\n', ...
        results(i,1), results(i,2), results(i,3), ...
        results(i,4), results(i,5), results(i,6), ...
        results(i,7));
end
