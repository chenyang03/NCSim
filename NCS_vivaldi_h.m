% F. Dabek, R. Cox, and F. Kaashoek. Vivaldi: A Decentralized Network Coordinate System. In Proc. of ACM SIGCOMM, 2004.

function [coord] = NCS_vivaldi_h(D_change, dim, N, K)

% parameters

% N - Number of hosts
% K - Number of neighbors

delta = 0.25;
ce = 0.25;
iteration = 2000;
neighbor = zeros(N,K);
for i=1:N
    coord(i,:)=zeros(1,dim);
    coord(i, dim) = 0.01;

    % Randomly selected Neighbors
      tmp = randperm(N);
      point = 1;
      for j=1:K
          if (i ~= tmp(point))
              neighbor(i, j) = tmp(point);
              point = point + 1;
          else
             while(D_change(i, tmp(point))== 0 || D_change(i, tmp(point)) < 0)
                point = point + 1;
             end
              neighbor(i, j) = tmp(point);
              point = point+1;
          end
      end    
%      tmp = randperm(N); 
%      neighbor(i,:)= tmp(1:K);
    error(i) = 1;
end
delta_temp = delta;
predicted_matrix = zeros(N,N);
for j = 1:iteration
    for i = 1:N
        % ID_neigh = ceil(N*rand);   %choose a neighbor of node i
        ID_neigh = neighbor(i, ceil(K*rand));
        temp_delta = delta;
        if (D_change(ID_neigh, i) <= 0 | D_change(i, ID_neigh) <= 0)
            continue;
        end;
%        if (Select_delta)
            w = error(i) / (error(i) + error(ID_neigh));
            e_s = abs(Dist_h(coord(i,:), coord(ID_neigh,:)) - D_change(i, ID_neigh)) / D_change(i, ID_neigh);
            error(i) = e_s * ce * w + error(i) * (1 - ce * w);
            if error(i) > 1
                error(i) = 1;
            end;
            delta = delta_temp * w;
%        end;
        Force = delta * ( D_change(i,ID_neigh) - Dist_h(coord(i,:),coord(ID_neigh,:)) ) * Direct_h(coord(i,:),coord(ID_neigh,:));
        coord(i,:) = coord(i,:) + Force;
%          coord(i, 1:(dim-1)) = coord(i, 1:(dim-1)) + Force(1:(dim-1));
%          coord(i, dim) = coord(i, dim) + Force(dim);
        if (coord(i, dim) < 0.01)
            coord(i, dim) = 0.01;
        end
    end;
    
end;
