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
allSolutions = zeros(numRuns, n);
allIterations = zeros(1, numRuns);
allLocalSearchIterations = zeros(1, numRuns);

fprintf('Executando GRASP %d vezes...\n', numRuns);

for i = 1:numRuns
    fprintf('\n--- Execução %d ---\n', i);
    [score, nodes, totalIterations, localSeachIterations, foundTime] = GRASP_SNS_Optimized(G, time, n, r, Cmax);
    [avgSP, maxSP] = PerfSNS(G, nodes);

    allScores(i) = avgSP;
    allMaxSP(i) = maxSP;
    allTimes(i) = foundTime;
    allSolutions(i, :) = nodes;
    allIterations(i) = totalIterations;
    allLocalSearchIterations(i) = localSeachIterations;

    fprintf('Média SP: %.4f | Max SP: %.4f | Tempo: %.2fs | Iterações Totais: %d | Iterações Local Search: %d\n', ...
            avgSP, maxSP, foundTime, totalIterations, localSeachIterations);
end

% --- Estatísticas ---
fprintf('\n====================\n');
fprintf('Resultados Finais (%d execuções)\n', numRuns);
fprintf('====================\n');
fprintf('Média função objetivo: %.4f\n', mean(allScores));
fprintf('Valor mínimo: %.4f\n', min(allScores));
fprintf('Valor máximo: %.4f\n', max(allScores));
fprintf('Tempo médio até melhor solução: %.2fs\n', mean(allTimes));
fprintf('Iterações totais: %.4f\n', mean(allIterations));
fprintf('Iterações Local Search: %.4f\n', mean(allLocalSearchIterations));

% --- Melhor solução encontrada ---
[~, bestIdx] = min(allScores);
bestNodes = allSolutions(bestIdx, :);

fprintf('\nMelhor solução encontrada (run %d):\n', bestIdx);
fprintf('Nodes: [%s]\n', num2str(sort(bestNodes)));

% --- Plot ---
figure;
plotTopology(Nodes, Links, bestNodes);
title('Melhor solução GRASP (entre 10 runs)');