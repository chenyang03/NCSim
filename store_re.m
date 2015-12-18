function output = store_re(input, range, precision)

m = ceil(input / range * precision);
for i=1:length(m)
    if (m(i) < 0)
        fprintf('%d warning!!', m(i));
    end
end
%m
count = zeros(1, precision);
total = zeros(1, precision);
for i = 1:length(m)
    if (m(i) < precision)
        count(m(i)+1) = count(m(i)+1) + 1;
    end
end
total_number = length(m);

for i = 1:1000
    if (i == 1)
        total(i) = count(1);
    else
        total(i) = total(i-1) + count(i);
    end
end

output = total/total_number;
%h3=plot(sorted(pos_array), m_array, str);
%h3 = plot(0:1/1000:1-1/1000, total/total_number, str);
