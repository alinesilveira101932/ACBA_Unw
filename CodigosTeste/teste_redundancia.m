%% Teste de redundancia

function t = teste_redundancia(S,b,c,d)

 % % Mostra as dimensões antes de qualquer operação
 %    fprintf('\n=== DEBUG teste_redundancia ===\n');
 %    fprintf('size(S): %d x %d\n', size(S,1), size(S,2));
 %    fprintf('size(b): %d x %d\n', size(b,1), size(b,2));
 %    fprintf('size(c): %d x %d\n', size(c,1), size(c,2));
 %    fprintf('size(d): %d x %d\n', size(d,1), size(d,2));
% A restrição c'*x <= d será redundante se e somente se t <= 0

Sa = [S;c'];
ba = [b;d+1]; % Para limitar superiormente o custo

x = linprog(-c,Sa,ba);
t = (c'*x - d);
