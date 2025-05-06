function s = find_zero_divide(func, min, max, minimal_interval, abs_thr)
    while true
        j_max = func(max);
        j_min = func(min);
        if j_min*j_max > 0
            max = max * 2;
        else
            break;
        end
    end
    while max-min > minimal_interval
        s = (min+max)/2;
        j = func(s);
        if abs(j_min)<abs_thr && abs(j_max)<abs_thr
            if j_min < 0
                s = min;
            else
                s = max;
            end
            break;
        elseif j*j_min > 0
            min = s;
            j_min = j;
        else
            max = s;
            j_max = j;
        end
    end
end