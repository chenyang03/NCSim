function [point] = seek_percentage(input_list, percentage)
    m = length(input_list);
    for i = 1:m
        if (input_list(i) > percentage)
            break;
        end
    end
    %return
    point = i*1.0/m;
    