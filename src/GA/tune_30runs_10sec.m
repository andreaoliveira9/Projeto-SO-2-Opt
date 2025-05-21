
% tune_ga.m
% Testa diferentes combinações de parâmetros para o algoritmo GA e identifica os melhores

addpath('SupportingFiles');
clear; clc;

% --- Parâmetros do problema ---
n = 12;
time = 10;
Cmax = 1000;
numRuns = 30;
populationSizes = [20, 50, 100, 150, 200];
mutationProbs = [0.05, 0.1, 0.2, 0.5, 0.7];
elitistParams = [1, 5, 10];

% --- Carregamento dos dados ---
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
G = graph(L);

% --- Resultados ---
results = [];

configId = 1;
fprintf('Testando diferentes configurações de parâmetros GA...\n');

for p = populationSizes
    for m = mutationProbs
        for e = elitistParams
            scores = zeros(1, numRuns);
            times = zeros(1, numRuns);
            evaluations = zeros(1, numRuns);
            validCounts = zeros(1, numRuns);
            discards = zeros(1, numRuns);

            parfor i = 1:numRuns
                [score, nodes, ~, foundTime, totalEval, totalVal, discarded] = GA_SNS(G, time, n, p, m, e, Cmax);
                [avgSP, ~] = PerfSNS(G, nodes);
                scores(i) = avgSP;
                times(i) = foundTime;
                evaluations(i) = totalEval;
                validCounts(i) = totalVal;
                discards(i) = discarded;
            end

            results = [results; p, m, e, ...
                min(scores), mean(scores), max(scores), ...
                mean(times), mean(evaluations), mean(validCounts), mean(discards)];

            fprintf('Pop=%d Mut=%.2f Elite=%d -> Min: %.4f | Média: %.4f | Max: %.4f | Tempo médio: %.4f\n', ...
                p, m, e, min(scores), mean(scores), max(scores), mean(times));
        end
    end
end

% --- Tabela final ---
fprintf('\n================== Resultados Finais GA ==================\n');
fprintf('Pop | Mut | Elite | MinSP | MedSP | MaxSP | Tempo | Aval | Válidas | Descartadas\n');
fprintf('-------------------------------------------------------------------------------\n');
for i = 1:size(results,1)
    fprintf('%3d | %.2f | %5d | %6.2f | %6.2f | %6.2f | %5.2f | %5.0f | %8.0f | %11.0f\n', ...
        results(i,1), results(i,2), results(i,3), ...
        results(i,4), results(i,5), results(i,6), ...
        results(i,7), results(i,8), results(i,9), results(i,10));
end
