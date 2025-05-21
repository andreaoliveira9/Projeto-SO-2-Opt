
% tune_grasp.m
% Testa diferentes valores de r para o algoritmo GRASP e identifica os melhores parâmetros

addpath('SupportingFiles');
clear; clc;

% --- Parâmetros do problema ---
n = 12;
time = 30;
Cmax = 1000;
numRuns = 10;
rValues = [1, 2, 3, 5, 7, 10, 20, 50];

% --- Carregamento dos dados ---
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
G = graph(L);

% --- Resultados por valor de r ---
results = [];

fprintf('Testando diferentes valores de r...\n');

for idx = 1:length(rValues)
    r = rValues(idx);
    scores = zeros(1, numRuns);
    times = zeros(1, numRuns);
    fprintf('\n--- r = %d ---\n', r);

    for i = 1:numRuns
        [score, nodes, ~, foundTime] = GRASP_SNS(G, time, n, r, Cmax);
        [avgSP, ~] = PerfSNS(G, nodes);
        scores(i) = avgSP;
        times(i) = foundTime;
        fprintf('Run %d: SP = %.4f | Tempo = %.2fs\n', i, avgSP, foundTime);
    end

    results(idx, :) = [r, min(scores), mean(scores), max(scores), mean(times)];

    fprintf('Resumo para r = %d -> Min: %.4f | Média: %.4f | Max: %.4f | Tempo médio: %.4f\n', ...
        r, min(scores), mean(scores), max(scores), mean(times));
end

% --- Tabela final ---
fprintf('\n==================== Resultados Finais ====================\n');
fprintf(' r |   Min SP   |   Média SP   |   Max SP   | Tempo Médio (s)\n');
fprintf('------------------------------------------------------------\n');
for i = 1:size(results,1)
    fprintf('%2d | %10.4f | %12.4f | %9.4f | %14.2f\n', ...
        results(i,1), results(i,2), results(i,3), results(i,4), results(i,5));
end
