function [coord] = vivaldi(D_change, dim, N, K)
% parameters
delta = 0.25;
ce = 0.25;
iteration = 2000;
%fprintf('Vivaldi\n');
coord = zeros(N, dim);
neighbor = zeros(N,K);

neighbor_count = zeros(N);

%N

fperror = [];
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
      if (neighbor_count(i) > K)
          neighbor_count(i) = K;
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
%           
%           
%       for j=1:neighbor_count(i)
%           
%           for ii=1:N
%               
%           end
%       % Actual Neighbor maybe less than K %
%           if (i ~= tmp(point))
%               neighbor(i, j) = tmp(point);
%               point = point + 1;
%           else
%              while(D_change(i, tmp(point))== 0 || D_change(i, tmp(point)) < 0)
%                 point = point + 1;
%              end
%               neighbor(i, j) = tmp(point);
%               point = point+1;
%           end
%       end    
% %      tmp = randperm(N); 
%      neighbor(i,:)= tmp(1:K);
    error(i) = 1;
end
% for i=1:N
%     coord(i,:)=zeros(1,dim);
%     error(i)=1;
% end

delta_temp = delta;
predicted_matrix = single(zeros(N,N));
%predicted_matrix = single(predicted_matrix);

for j = 1:iteration
    for i = 1:N
        % ID_neigh = ceil(N*rand);   %choose a neighbor of node i
        if (neighbor_count(i) == 0)
            continue;
        end
        ID_neigh = neighbor(i, ceil(neighbor_count(i)*rand));
        temp_delta = delta;
%        ID_neigh = ceil(N*rand);
%    D_change(ID_neigh, i)
        if (D_change(ID_neigh, i) <= 0 | D_change(i, ID_neigh) <= 0)
            continue;
        end;

%        if (Select_delta)
            w = error(i) / (error(i) + error(ID_neigh));
            e_s = abs(Dist(coord(i,:), coord(ID_neigh,:)) - D_change(i, ID_neigh)) / D_change(i, ID_neigh);
            error(i) = e_s * ce * w + error(i) * (1 - ce * w);
%            if error(i) > 1
%                error(i) = 1;
%            end;
            delta = delta_temp * w;
%        end;
        Force = delta * ( D_change(i,ID_neigh) - Dist(coord(i,:),coord(ID_neigh,:)) ) * Direct(coord(i,:),coord(ID_neigh,:));
        coord(i,:) = coord(i,:) + Force;
    end;
    
    %%%%%%%%%    
end;

%coord
%result_matrix

%fprintf('\n');