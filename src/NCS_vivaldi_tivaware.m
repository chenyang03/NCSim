% G. Wang, B. Zhang, T.S.E. Ng. Towards Network Triangle Inequality Violation Aware Distributed Systems. In Proc. of ACM IMC, 2007.

function [coord, fperror_newhost] = NCS_vivaldi_tivaware(D_change, dim, N, K, converge_on)
% parameters
delta = 0.25;
ce = 0.25;
%iteration = 2000;
%fprintf('Vivaldi\n');


% in TIV aware Vivaldi, each host has 2K neighbors...
neighbor = zeros(N,2*K);
predicted_matrix = zeros(N, N);
host_age = zeros(1, N);

neighbor_count = zeros(N);
fperror_flashcrowd = [];
new_host_scale = 20;

% evaluate the convergence? %
if (converge_on == 0)
    % no: only run for 30 rounds
    iteration = 2000;
else
    % yes: run more 30 rounds and see the convergence procedure of 10 newly
    % joinging hosts
    
    iteration = 4000;
end

fperror_newhost = [];
% N
for i=1:N
      coord(i,:)=zeros(1,dim);
     % Randomly selected Neighbors
      tmp = randperm(N);
      point = 1;
      for j=1:N
          if (D_change(i, j) >= 0 && i ~= j)
              neighbor_count(i) = neighbor_count(i) + 1;
          end
      end
      if (neighbor_count(i) > 2*K)
          neighbor_count(i) = 2*K;
      end
      j = 1;
      for ii = 1:N
          if (D_change(i, tmp(ii)) >= 0 && i ~= tmp(ii))
              neighbor(i, j) = tmp(ii);
              j=j+1;
              if (j > neighbor_count(i))
                  break;
              end              
          end
      end
    error(i) = 1;
end
% for i=1:N
%     coord(i,:)=zeros(1,dim);
%     error(i)=1;
% end

delta_temp = delta;
predicted_matrix = zeros(new_host_scale,N);

for j = 1:iteration
    if (converge_on == 1 & j == floor(iteration/2));
        rand_seq = randperm(N);
        new_hosts =  rand_seq(1:new_host_scale);
        existing_hosts = rand_seq(new_host_scale+1:N);
        % Remove the NC values of the 10 new hosts        
        for ii = 1:new_host_scale
            coord(new_hosts(ii),:)=zeros(1,dim);
            error(new_hosts(ii)) = 1;
            host_age(new_hosts(ii)) = 0;
        end
        % new
    end
    
    for i = 1:N
        host_age(i) = host_age(i) + 1;
        ID_neigh = neighbor(i, ceil(K*rand));
        while (ID_neigh <= 0)
            ID_neigh = neighbor(i, ceil(K*rand));
        end
        
        temp_delta = delta;
        if not(D_change(ID_neigh, i) <= 0 | D_change(i, ID_neigh) <= 0)
            w = error(i) / (error(i) + error(ID_neigh));
            e_s = abs(Dist(coord(i,:), coord(ID_neigh,:)) - D_change(i, ID_neigh)) / D_change(i, ID_neigh);
            error(i) = e_s * ce * w + error(i) * (1 - ce * w);
            delta = delta_temp * w;
            Force = delta * ( D_change(i,ID_neigh) - Dist(coord(i,:),coord(ID_neigh,:)) ) * Direct(coord(i,:),coord(ID_neigh,:));
            coord(i,:) = coord(i,:) + Force;
        end
        
        if (mod(host_age(i), 100) == 0)
            % neighbor - renew %
            % score each neighbor
            shurnk_score = zeros(1, 2*K);
            for kk=1:2*K
                if (neighbor(i, kk) == 0)
                    shurnk_score(kk) = -1;
                    continue;
                end
                measured_dist = D_change(i, neighbor(i, kk));
                predicted_dist = Dist(coord(i, :), coord(neighbor(i, kk)));
                score(kk) = predicted_dist / (measured_dist+eps);                
                if (measured_dist < 0)
                    shurnk_score(kk) = -1;
                end
            end
            % remove bad neighbor
            % If the prediction ratio of an edge is very small, it means this edge is shrunk a lot and it is more likely to cause severe TIV. So we remove the 32
            % nodes with smallest prediction ratios among the 64 neighbor candidates, and the remaining 32 nodes are used as the neighbors in the next iteration.
            [yy, ii] = sort(shurnk_score);
            % ii(1...2*K); the first 1...K will be removed
            smaller_neighbor_set = neighbor(i, ii(K+1:2*K));            
            % involve more new neighbor
            neighbor(i, 1:K) = smaller_neighbor_set;
            for ii=K+1:2*K
                % select neighbor(i, ii)
                while(1)
                    tmp = ceil(N*rand);
                    ok = 1;
                    for jj=1:ii-1
                        if (tmp == neighbor(i, jj))
                            ok = 0;
                            break;
                        end
                    end
                    if (ok == 1)                        
                        break;
                    end
                end
                neighbor(i, ii) = tmp;
            end
        end
        
    end;
    
    
    if (j >= floor(iteration/2) & mod((j-floor(iteration/2)), 64) == 0 & converge_on == 1)
        predict_matrix_newhosts = zeros(new_host_scale, N);
        for x = 1:new_host_scale
            for y = 1:N
                predicted_matrix_newhosts(x,y) = Dist(coord(new_hosts(x),:), coord(y,:)); 
            end;
        end;        
        abs_err = absolute_error(predicted_matrix_newhosts, D_change(new_hosts, :));
        fperror_newhost = [fperror_newhost; median(abs_err)];
    end

    
end;

if (converge_on == 1)
    fprintf('New Host: ');
    for i=1:length(fperror_newhost)
        fprintf('%.3f ', fperror_newhost(i));
    end
    fprintf('\n');
end