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

fprintf('Executando GRASP %d vezes...\n\n', numRuns);

for i = 1:numRuns
    fprintf('--- Execução %d ---\n', i);
    [score, nodes, iterations, foundTime] = GRASP_SNS(G, time, n, r, Cmax);
    [avgSP, maxSP] = PerfSNS(G, nodes);

    allScores(i) = avgSP;
    allMaxSP(i) = maxSP;
    allTimes(i) = foundTime;
    allIterations(i) = iterations;
    allSolutions(i, :) = nodes;

    fprintf('Média SP: %.4f | Max SP: %.4f | Tempo: %.2fs | Iterações: %d\n\n', ...
            avgSP, maxSP, foundTime, iterations);
end

% --- Estatísticas ---
fprintf('====================\n');
fprintf('Resultados Finais (%d execuções)\n', numRuns);
fprintf('====================\n');
fprintf('Média dos valores da função objetivo: %.4f\n', mean(allScores));
fprintf('Valor mínimo da função objetivo: %.4f\n', min(allScores));
fprintf('Valor máximo da função objetivo: %.4f\n', max(allScores));
fprintf('Tempo médio até melhor solução: %.2fs\n', mean(allTimes));
fprintf('Iterações médias: %.2f\n', mean(allIterations));

% --- Melhor solução encontrada ---
[~, bestIdx] = min(allScores);
bestNodes = allSolutions(bestIdx, :);

fprintf('\nMelhor solução encontrada (run %d):\n', bestIdx);
fprintf('Nodes: [%s]\n', num2str(sort(bestNodes)));

% --- Plot ---
figure;
plotTopology(Nodes, Links, bestNodes);
title('Melhor solução GRASP (entre 10 runs)');