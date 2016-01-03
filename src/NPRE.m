function NPRE_Value = NPRE(input_vector)
    temp = sort(input_vector);
%    length(temp)
    NPRE_Value = temp(ceil(length(temp)*0.9));