function [final_rank_accuracy_result] = rank_accuracy(estimate, real)
    
    percentage_seq = [0.01:0.01:0.1 0.2:0.1:1];
       estimate = max(eps, estimate);
    raw_rtt = real(:);
    tmp = find(raw_rtt<=0);
    raw_rtt(tmp) = [];
    raw_rtt = ceil(raw_rtt);
    sorted_raw_rtt = sort(raw_rtt);
    
    %   1    2   3     4   5    6    7     8   9   10  11  12  13  14  15
    % 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6
    
    %  16  17  18  19
    %  0.7 0.8 0.9 1
       
    
    threshold_distance_seq = ceil(sorted_raw_rtt(ceil(percentage_seq*length(raw_rtt))));
    

    max_distance = ceil(max(raw_rtt));
    
    range_accuracy_result = zeros(1, max_distance); % from 1ms to 300ms
%    reverse_range_accuracy_result = zeros(1, max_distance); % from 1ms to 300ms
    number_measured_links = zeros(1, max_distance);
    number_predicted_links = zeros(1, max_distance);
    number_measured_predicted_links = zeros(1, max_distance);
    
    [m n] = size(estimate);
    for i=1:m
        for j=1:n
            if (real(i, j) <=0 | i == j)
                continue;
            end
            tmp = ceil(real(i, j));
            for ii = tmp:max_distance
                number_measured_links(ii) = number_measured_links(ii) + 1;
            end
            tmp2 = max(ceil(real(i, j)), ceil(estimate(i, j)));
            for ii = tmp2:max_distance
                number_measured_predicted_links(ii) = number_measured_predicted_links(ii) + 1;
            end
            tmp3 = ceil(estimate(i, j)+0.01);
            if (tmp3 <= 0)
                fprintf('Warning: %d\n', tmp3);
            end
            %fprintf('%.2f ', estimate(i, j));
            for ii=tmp3:max_distance
                number_predicted_links(ii) = number_predicted_links(ii) + 1;
            end
        end
    end
    for i=1:max_distance
        range_accuracy_result(i) = number_measured_predicted_links(i) / number_measured_links(i);
        if (number_predicted_links(i) > 0)
            reverse_range_accuracy_result(i) = number_measured_predicted_links(i) / number_predicted_links(i);
        end
    end
    fprintf('(%.2f %.2f %.2f %.2f)\n', range_accuracy_result(threshold_distance_seq(1)), range_accuracy_result(threshold_distance_seq(5)), range_accuracy_result(threshold_distance_seq(10)), range_accuracy_result(threshold_distance_seq(14)));
    %fprintf('(%.2f %.2f %.2f)\n', reverse_range_accuracy_result(20), reverse_range_accuracy_result(40), reverse_range_accuracy_result(80));
    final_rank_accuracy_result = [];
     for i=1:length(threshold_distance_seq)
         final_rank_accuracy_result = [final_rank_accuracy_result range_accuracy_result(threshold_distance_seq(i))];    
     end
    