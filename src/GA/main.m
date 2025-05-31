% main.m
% Execução de 10 runs do algoritmo Genético com verificação de Cmax
addpath('SupportingFiles');
clear; clc;

% --- Parâmetros do problema ---
n = 12;                 % número de nós a selecionar
populationSize = 100;   % tamanho da população
mutationProb = 0.2;     % probabilidade de mutação
elitistParam = 10;      % parâmetro elitista (qtd de pais que sobrevivem)
time = 30;              % tempo de execução por execução
Cmax = 1000;            % limite máximo permitido entre pares de servidores
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

fprintf('Executando GA %d vezes...\n', numRuns);

for i = 1:numRuns
    fprintf('\n--- Execução %d ---\n', i);
    [score, nodes, ~, foundTime] = GA_SNS(G, time, n, populationSize, mutationProb, elitistParam, Cmax);
    [avgSP, maxSP] = PerfSNS(G, nodes);

    allScores(i) = avgSP;
    allMaxSP(i) = maxSP;
    allTimes(i) = foundTime;
    allSolutions(i, :) = nodes;

    fprintf('Média SP: %.4f | Max SP: %.4f | Tempo: %.2fs\n', ...
            avgSP, maxSP, foundTime);
end

% --- Estatísticas ---
fprintf('\n====================\n');
fprintf('Resultados Finais (%d execuções)\n', numRuns);
fprintf('====================\n');
fprintf('Média função objetivo: %.4f\n', mean(allScores));
fprintf('Valor mínimo: %.4f\n', min(allScores));
fprintf('Valor máximo: %.4f\n', max(allScores));
fprintf('Tempo médio até melhor solução: %.2fs\n', mean(allTimes));

% --- Melhor solução encontrada ---
[~, bestIdx] = min(allScores);
bestNodes = allSolutions(bestIdx, :);

fprintf('\nMelhor solução encontrada (run %d):\n', bestIdx);
fprintf('Nodes: [%s]\n', num2str(sort(bestNodes)));

% --- Plot ---
figure;
plotTopology(Nodes, Links, bestNodes);
title('Melhor solução GA (entre 10 runs)');