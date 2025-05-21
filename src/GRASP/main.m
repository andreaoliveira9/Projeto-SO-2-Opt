% main.m
% Execução de 10 runs do algoritmo GRASP com verificação de Cmax

addpath('SupportingFiles');
clear; clc;

% --- Parâmetros do problema ---
n = 12;             % número de nós a selecionar
r = 3;              % tamanho da RCL
time = 30;          % tempo de execução por execução
Cmax = 1000;        % limite máximo permitido entre pares de servidores
numRuns = 10;

% --- Carregamento dos dados ---
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
G = graph(L);

% --- Resultados ---
allScores = zeros(1, numRuns);
allMaxSP = zeros(1, numRuns);
allTimes = zeros(1, numRuns);
allIterations = zeros(1, numRuns);
allSolutions = zeros(numRuns, n);
allEvaluated = zeros(1, numRuns);
allValid = zeros(1, numRuns);
allDiscarded = zeros(1, numRuns);

fprintf('Executando GRASP %d vezes...\n\n', numRuns);

for i = 1:numRuns
    fprintf('--- Execução %d ---\n', i);
    [score, nodes, iterations, foundTime, totalEval, totalVal, discarded] = ...
        GRASP_SNS(G, time, n, r, Cmax);
    [avgSP, maxSP] = PerfSNS(G, nodes);

    allScores(i) = avgSP;
    allMaxSP(i) = maxSP;
    allTimes(i) = foundTime;
    allIterations(i) = iterations;
    allSolutions(i, :) = nodes;
    allEvaluated(i) = totalEval;
    allValid(i) = totalVal;
    allDiscarded(i) = discarded;

    fprintf('Média SP: %.4f | Max SP: %.4f | Tempo: %.2fs | Iterações: %d\n', ...
            avgSP, maxSP, foundTime, iterations);
    fprintf('Total Avaliados: %d | Válidos: %d | Descartados: %d\n\n', ...
            totalEval, totalVal, discarded);
end

% --- Estatísticas ---
fprintf('====================\n');
fprintf('Resultados Finais (%d execuções)\n', numRuns);
fprintf('====================\n');
fprintf('Média função objetivo: %.4f\n', mean(allScores));
fprintf('Valor mínimo: %.4f\n', min(allScores));
fprintf('Valor máximo: %.4f\n', max(allScores));
fprintf('Tempo médio até melhor solução: %.2fs\n', mean(allTimes));
fprintf('Iterações médias: %.2f\n', mean(allIterations));
fprintf('Soluções avaliadas (média): %.1f\n', mean(allEvaluated));
fprintf('Soluções válidas (média): %.1f\n', mean(allValid));
fprintf('Descartadas por Cmax (média): %.1f\n', mean(allDiscarded));

% --- Melhor solução encontrada ---
[~, bestIdx] = min(allScores);
bestNodes = allSolutions(bestIdx, :);

fprintf('\nMelhor solução encontrada (run %d):\n', bestIdx);
fprintf('Nodes: [%s]\n', num2str(sort(bestNodes)));

% --- Plot ---
figure;
plotTopology(Nodes, Links, bestNodes);
title('Melhor solução GRASP (entre 10 runs)');