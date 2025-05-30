% Carregar dados
Nodes = load('SupportingFiles/Nodes200.txt');
Links = load('SupportingFiles/Links200.txt');
L = load('SupportingFiles/L200.txt');
n_servers = 12; % Número de servidores a selecionar
Cmax = 1000; % Distância máxima entre controladores

G = graph(L);
D = distances(G); % Matriz de distâncias mínimas

N = size(Nodes, 1);

fid = fopen('ILP/opt_problem.lp', 'w');

% Função objetivo
fprintf(fid, 'min: ');
terms = [];
for s = 1:N
    for i = 1:N
        if s ~= i
            terms{end+1} = sprintf('%.6f*g_%d_%d', D(s,i), s, i);
        end
    end
end
fprintf(fid, '%s;\n', strjoin(terms, ' + '));
fprintf(fid, '\n');

% Restrição: exatamente n servidores
fprintf(fid, '%s = %d;\n', strjoin(arrayfun(@(i) sprintf('z_%d', i), 1:N, 'UniformOutput', false), ' + '), n_servers);
fprintf(fid, '\n');

% Restrição: cada nó é servido por exatamente um servidor
for s = 1:N
    fprintf(fid, '%s = 1;\n', strjoin(arrayfun(@(i) sprintf('g_%d_%d', s, i), 1:N, 'UniformOutput', false), ' + '));
end
fprintf(fid, '\n');

% Restrição: só pode ser servido por um nó que seja servidor
for s = 1:N
    for i = 1:N
        fprintf(fid, 'g_%d_%d - z_%d <= 0;\n', s, i, i);
    end
end
fprintf(fid, '\n');

% Restrição: distância máxima entre controladores não pode exceder Cmax
for i = 1:N
    for j = i+1:N
        if D(i,j) > Cmax
            fprintf(fid, 'z_%d +  z_%d <= 1;\n', i, j);
        end
    end
end
fprintf(fid, '\n');

% Variáveis binárias
fprintf(fid, '\n');
for i = 1:N
    fprintf(fid, 'bin z_%d;\n', i);
end
for s = 1:N
    for i = 1:N
        fprintf(fid, 'bin g_%d_%d;\n', s, i);
    end
end

fclose(fid);
disp('Ficheiro LP gerado: opt_problem.lp');