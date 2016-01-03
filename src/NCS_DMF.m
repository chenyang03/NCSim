% Y. Liao, P. Geurts and G. Leduc. Network Distance Prediction Based on Decentralized Matrix Factorization. Proc. of IFIP Networking 2010. 

function [out_host, in_host, fperr] = NCS_DMF(B_change, dim, K, iteration, lambda, show_log_on)

fperr = [];
C=15;
%out_host_cont = out_host;
%in_host_cont = in_host;
%delta = 0.0000002;
%K=32;
N=length(B_change);

if (N < K)
    % NMF Directly %
    [out_host, in_host] = mysvd(B_change, dim);
   return;
end

neighbor = zeros(N,K);
error_in = zeros(1,N);
error_out = zeros(1,N);
error = zeros(1, N);

for i=1:N
    out_host(i,:)=rand(1,dim);
    in_host(:,i)=rand(dim,1);
    tmp = randperm(N); 

    current_point = 1;
    neighbor(i, :) = tmp(1:K);
end
if (show_log_on == 1)
    fprintf('\n');
end
for round=1:iteration
    for i=1:N
        % d_i_to / d_i_from
        %neighbor(i, :)
        d_i_to = B_change(i, neighbor(i, :));
        d_i_from = B_change(neighbor(i, :), i)';
        % X_{i}, Y_{i}
        Xi = out_host(neighbor(i, :), :);
        Yi= in_host(:, neighbor(i, :))';
        % C = MRDIVIDE(A,B)
        % A*INV(B)


        lambda = 50;
        out_host(i, :) = mrdivide(d_i_to * Yi, Yi' * Yi + lambda * eye(dim));
        in_host(:, i) = mrdivide(d_i_from * Xi, Xi' * Xi + lambda * eye(dim))';
        
    end
    if (show_log_on == 1)
        tmp_npre = NPRE(relative_error(out_host * in_host, B_change));    
        fprintf('%.3f ', tmp_npre);    
        fperr = [fperr tmp_npre];
    end
end
%fprintf('\nDMF        (Dimension # = %d, lambda = %3d): <50th, 90th, Avg.>', dim, lambda);
%fprintf('%.3f, %.3f, %.3f', median(relative_error(out_host * in_host, B_change)), NPRE(relative_error(out_host * in_host, B_change)), mean(relative_error(out_host * in_host, B_change)));
