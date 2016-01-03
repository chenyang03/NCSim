% a simplified Phoenix NC system without considering node churn, distance variation, etc
function [out_host, in_host, fpre_newhost, fpre_flashcrowd] = phoenix(D, dim, N, K, C, converge_on, new_hosts)

% N: Number of nodes
% K: Number of Neighbors

D_change = D;

fpre_newhost=0;
fpre_flashcrowd=0;

% Parameters
C = 10;
%new_host_scale = 20;
if (converge_on == 1)
    new_host_scale = length(new_hosts);
end

% Per Round -> Per Second ...

% evaluate the convergence? %
if (converge_on == 0)
    % no: only run for 30 rounds
    round_bound = 30;
else
    % yes: run more 30 rounds and see the convergence procedure of 10 newly
    % joinging hosts
    
    round_bound = 60;
end

% 30*64 = 
closest_num =0;
result_matrix_seq = [];

if (N == 1)
%    [out_host, in_host] = NMF(D, dim);
    out_host = zeros(1, dim);
    in_host = zeros(dim, 1);
    result_matrix_seq = [];
    if (converge_on == 1)
        predicted_matrix = out_host*in_host;
        for round=1:round_bound
            result_matrix_seq = [result_matrix_seq; D(1, 1)];
        end
    end
    return
end

if (N < K)
    % NMF Directly %
    length(D);
    [out_host, in_host] = NMF(D, dim);
    result_matrix_seq = [];
    if (converge_on == 1)
        predicted_matrix = out_host*in_host;
        for round=1:round_bound
            result_matrix_seq = [result_matrix_seq; predicted_matrix];
        end
    end
    
   return;
end

neighbor = zeros(N,K);
error_in = zeros(1,N);
error_out = zeros(1,N);
for i=1:N
    out_host(i,:)=rand(1,dim);
    in_host(:,i)=rand(dim,1);

    % Randomly selected Neighbors
      tmp = randperm(N);
      
      point = 1;
      for j=1:K
          neighbor(i, j) = -1;
      end
      for j=1:K
          % fill in neighbor(i, j) %
          if (i ~= tmp(point) && D_change(i, tmp(point)) > 0)
              neighbor(i, j) = tmp(point);
              point = point + 1;
          else
             while(i == tmp(point) || D_change(i, tmp(point)) <= 0) % || D_change(i, tmp(point)) < 0)
                point = point + 1;
                if (point > N)
                    break;
                end
             end
             
%              printf('%d ', j);
%              neighbor(i, j) = tmp(point);
%              point = point+1;
          end
              if (point > N)
                  break;
              end              
      end
      

end

% Out : X in the paper, In : Y in the paprt
%delta_out = delta;
%delta_in = delta;

%for j=1:15

w=zeros(N, K, dim);
h=zeros(N, dim, K);
D_host2landmark=zeros(N, K);
D_host2landmark_out=zeros(N, K);
D_host2landmark_in=zeros(N, K);

%for j = 1:100
%for j = 1:(iteration/2)
fpre_newhost = [];
fpre_flashcrowd = [];

for round = 1:round_bound
%    fprintf('22');
    if (round == round_bound)        
        % Stat %
        average_node_weight = zeros(1, N);
        node_weight_count = zeros(1, N);
        %%%%%%%%%%%%
    end
    if (converge_on == 1 & round == floor(round_bound/2));
%        rand_seq = randperm(N);
%        new_hosts =  rand_seq(1:new_host_scale);
        % Remove the NC values of the 10 new hosts        
        for ii = 1:new_host_scale
            out_host(new_hosts(ii),:)=rand(1,dim);
            in_host(:,new_hosts(ii))=rand(dim,1);
        end
        % new
    end

%     
    if (round == 1)
        % First round : Joining %

        for i=1:N
            tmp = randperm(N);
            neighbor(i, :) = tmp(1:K);
        end
    end

    
    
    for i = 1:N
        
        % K: number of landmarks
        % w: Kxd, h: dxK. The position vectors of all landmarks

%         if (round == 1 & i < K)
%             continue; % Already has NC
%         end

        new_hit = 0; % if a new joining host in convergence test
        
        if (converge_on == 1 & round == floor(round_bound/2))
            for kk = 1:new_host_scale
                if (new_hosts(kk) == i)
                    new_hit = 1;
                    break;
                end
            end
        end
        
        if (round == 1 | new_hit == 1) % Initial NC Calculation
            target_host = neighbor(i, :);
            %target_host = reference_hosts(i, :);
            target_host = target_host(find(target_host>0));
            target_host = target_host(find(D(i, target_host)>=0));
            %target_host
            
            weight_out_vec = zeros(1, length(target_host)); % K -> actual_K
            weight_in_vec = zeros(1, length(target_host)); % K -> actual_K

            temp_w = out_host(target_host,:);
            temp_h = in_host(:, target_host);
            temp_D_host2landmark = D_change(i, target_host);
             temp_D_host2landmark_out = temp_D_host2landmark(1:length(target_host));
             temp_D_host2landmark_in = temp_D_host2landmark(1:length(target_host));

            for ii=1:1:length(target_host)
                if (D(i,  target_host(ii)) < 0)
                    weight_out_vec(ii) = eps;             
                    weight_in_vec(ii) = eps;           
                else
                    weight_out_vec(ii) = 1;             
                    weight_in_vec(ii) = 1;           
                end
            end
%             size(temp_h')
%             size(temp_D_host2landmark_in')
%             size(sqrt(weight_in_vec)')
            t = weight_lsqnonneg(temp_h', temp_D_host2landmark_in', sqrt(weight_in_vec)');
            out_host(i, :) = t';        
            in_host(:, i) = weight_lsqnonneg(temp_w, temp_D_host2landmark_out', sqrt(weight_out_vec)');
        end
             
        % target 
        
        target_host = neighbor(i, :);
        %target_host = reference_hosts(i, :);
        %target_host
        target_host = target_host(find(target_host>0));
        %target_host
        
        actual_K = length(target_host);
        
        temp_w = out_host(target_host,:);
        temp_h = in_host(:, target_host);
        temp_D_host2landmark = D_change(i, target_host);
        


        % get the score of all hosts in neighbor(i, :) %
        score_out_vec = zeros(1, actual_K); % K -> actual_K
        score_in_vec = zeros(1, actual_K); % K -> actual_K
        score_aver_vec = zeros(1, actual_K); % K -> actual_K
        
        weight_out_vec = zeros(1, actual_K); % K -> actual_K
        weight_in_vec = zeros(1, actual_K); % K -> actual_K

        
        for index_nb=1:actual_K % K -> actual_K

            predict_ii_in = temp_w(index_nb, :) * in_host(:, i);
            predict_ii_out = out_host(i, :) * temp_h(:, index_nb);

 
           s1 = abs(predict_ii_out - D_change(i, target_host(index_nb)));% / (D_change(i, neighbor(i, index_nb))+eps);
           s2 = abs(predict_ii_in - D_change(i, target_host(index_nb)));% / (D_change(i, neighbor(i, index_nb))+eps);

           score_out_vec(index_nb) = s1;
           score_in_vec(index_nb) = s2;

        end
        
        out_threshold = median(score_out_vec);
        in_threshold = median(score_in_vec);


        for ii=1:actual_K % K -> actual_K
            if (score_out_vec(ii) < out_threshold)
                weight_out_vec(ii) = 1;
            else
                if (score_out_vec(ii) < out_threshold * C)% || score_out_vec(ii)  < out_upbound_threshold)
                    weight_out_vec(ii) = (out_threshold/score_out_vec(ii))^2;           
                else
                    weight_out_vec(ii) = eps;
                end            
            end
            if (score_in_vec(ii) < in_threshold)
                weight_in_vec(ii) = 1;
            else
                if (score_in_vec(ii) < in_threshold * C)% || score_in_vec(ii)  < in_upbound_threshold)
                    weight_in_vec(ii) = (in_threshold/score_in_vec(ii))^2;           
%                    weight_in_vec(ii) = 1 - score_in_vec(ii);
                else
                    weight_in_vec(ii) = eps;%1/C;
                end 
            end
        end

%         if (round == round_bound)                
%             
%             for ii=1:actual_K
%                 tmp_seq = target_host(ii);
%                 %tmp_seq
%                 average_node_weight(tmp_seq) = average_node_weight(tmp_seq) + weight_in_vec(ii) + weight_out_vec(ii);
%                 node_weight_count(tmp_seq) = node_weight_count(tmp_seq) + 1;
%             end
%         end
        
        
        
%         temp_w = temp_w(index_nblist_out, :);
%         temp_h = temp_h(:, index_nblist_in);
%         temp_D_host2landmark = temp_D_host2landmark(index_nblist);
         
         temp_D_host2landmark_out = temp_D_host2landmark(1:actual_K); % K -> actual_K
         temp_D_host2landmark_in = temp_D_host2landmark(1:actual_K); % K -> actual_K
         
      
        
        % x = lsqnonneg(C,d) 
        % returns the vector x that minimizes norm(C*x-d) subject to x >= 0. C and d must be real.

        %  out_host(i, :) [1*d] * In_NCs[d*K] = D_host2landmark[1*K];        
        %  => out_host(i, :) * h = D_hosts2landmark[1*K];
        %  => h' * out_host(i, :)' = D_host2landmark';
%        size(temp_h)
%        size(temp_D_host2landmark)



%        t = lsqnonneg(temp_h', temp_D_host2landmark_in');
        t = weight_lsqnonneg(temp_h', temp_D_host2landmark_in', sqrt(weight_in_vec)');
        out_host(i, :) = t';
        
        %  Out_NCs[K*d] * in_host(:, i) [d*1] = D_host2landmark'[K*1];
        %  => w * in_host(:, i) = D_hosts2landmark'[K*1];
        
%        in_host(:, i) = lsqnonneg(temp_w, temp_D_host2landmark_out');
        in_host(:, i) = weight_lsqnonneg(temp_w, temp_D_host2landmark_out', sqrt(weight_out_vec)');

%         w = backup_w;
%         h = backup_h;
%         D_host2landmark = backup_D_host2landmark;
    
    end;
    
    % Stat %
%     if (converge_on == 1 & round <= floor(round_bound/2))
%         predicted_matrix = out_host*in_host;
%         rerr = absolute_error(predicted_matrix, D);        
%         fpre_flashcrowd = [fpre_flashcrowd; median(rerr)];
%     end
    if (converge_on == 1 & round >= floor(round_bound/2))
        %predicted_matrix = out_host*in_host;
        predicted_matrix = out_host(new_hosts, :) * in_host(:, :);
        rerr = absolute_error(predicted_matrix, D(new_hosts, :));
        %rerr=absolute_error(predicted_matrix, D);
        fpre_newhost = [fpre_newhost; median(rerr)];
    end
end;

% normalize_weight = average_node_weight ./ (node_weight_count+eps)
% all_err = Detailed_Error_SVD(D, dim);
% 
% for i=1:N
%     fprintf('(%.3f, %.3f) ', all_err(i), normalize_weight(i));
% end

% Error - Weight - Relationship
%average_node_weight./node_weight_count

%re=(abs(out_host * in_host-D)./(D+0.1));
%mean(re)

%figure;
%plot(mean(re), average_node_weight./node_weight_count)

if (converge_on == 1)
    fprintf('New Host: ');
    for i=1:floor(round_bound/2)
        fprintf('%.3f ', fpre_newhost(i));
    end
    fprintf('\n');
%     fprintf('FlashCrowd: ');
%     for i=1:floor(round_bound/2)
%         fprintf('%.3f ', fpre_flashcrowd(i));
%     end
%     fprintf('\n');

end

% rerr=relative_error(out_host * in_host, D);
% fprintf('Lambda = %.3f, Relative Error: %.3f\n', lambda, NPRE(rerr));
fpre_flashcrowd = [];