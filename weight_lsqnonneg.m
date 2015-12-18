function x = weight_lsqnonneg(C,d,w)
% C - n * n
% d - n * 1
% w - n * 1
% Minimize the || (C * x - d) .* w ||
% <=>>.333
sw=sqrt(w);
%sw=w;
lambda = 2;
CC=diag(sw) * C; dd = d.*sw;
[m n] = size(CC);
CC = [CC; lambda * eye(n)]; dd = [dd; zeros(n, 1)];
x = lsqnonneg(CC, dd);

% options.print_level = 2;    % Increase print level.
% 
% % mm=5; nn= 5;
% % C_mat = rand(mm, nn);
% % d_vec = rand(nn, 1);
% 
% C_mat = diag(sw) * C; d_vec = d.*sw;
% [m n] = size(C_mat);
% 
% %C_mat - 32 * 8
% %d_vec - 32 * 1
% 
% % [x, g, info] = bcls( A, b, bl, bu, x, c, damp, options )
% 
% options.usolve = ; 
% [x, g, info] = bcls(C_mat, d_vec, zeros(n, 1), inf*ones(n, 1), rand(n, 1), [], 0, options);