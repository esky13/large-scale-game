function to = find_ti(xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu)
    delta_x = [xA0-xD0;yA0-yD0];
    delta_v = [vAx0-vDx0;vAy0-vDy0];
    func = @(t)cal_delta_o(t,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
    d_func = @(t)cal_d_delta_o(t,delta_x,delta_v,uA,uD,mu);
    dd_func = @(t)cal_dd_delta_o(t,delta_x,delta_v,uA,uD,mu);
    if delta_v'*delta_v+mu*delta_x'*delta_v < 0
        to = find_zero_divide(func,0,mu*norm(delta_x)/(uA+uD)+1/mu);
    elseif delta_x'*delta_v > 0
        to = find_zero_divide(func,0,log(2)/mu);
    else
        t_d_peak = find_zero_divide(dd_func,0,log(2)/mu);
        j = d_func(t_d_peak);
        if j < 0
            to = find_zero_divide(func,0,1);
        else
            t_d_zero1 = find_zero_divide(d_func,0,t_d_peak);
            t_d_zero2 = find_zero_divide(d_func,t_d_peak,2*t_d_peak);
            
            local_min = func(t_d_zero1);
            local_max = func(t_d_zero2);
            if local_min > 0
                to = find_zero_divide(func,t_d_zero2,2*t_d_zero2);
            elseif local_max < 0
                to = find_zero_divide(func,0,t_d_zero1);
            else
                t1 = find_zero_divide(func,0,t_d_zero1);
                t2 = find_zero_divide(func,t_d_zero1,t_d_zero2);
                t3 = find_zero_divide(func,t_d_zero2,2*t_d_zero2);
                to = [t1,t2,t3];
            end
        end
    end
end

function delta_o = cal_delta_o(t,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu)
    [rAc, xAc, yAc] = getIsochrones(t,xA0,yA0,vAx0,vAy0,uA,mu);
    [rDc, xDc, yDc] = getIsochrones(t,xD0,yD0,vDx0,vDy0,uD,mu);
    delta_o = (xAc-xDc)^2+(yAc-yDc)^2-(rAc-rDc)^2;
end

function d_delta_o = cal_d_delta_o(t,delta_x,delta_v,uA,uD,mu)
    e = exp(-mu*t);
    d_delta_o = e*(delta_v'*delta_v/mu*(1-e)+delta_x'*delta_v)... 
        -(uA-uD)^2/mu^2*(t-(1-e)/mu)*(1-e);
end

function dd_delta_o = cal_dd_delta_o(t,delta_x,delta_v,uA,uD,mu)
    e = exp(-mu*t);
    dd_delta_o = e*(delta_v'*delta_v*(2*e-1)-mu*delta_x'*delta_v)... 
        -(uA-uD)^2/mu^2*((1-e)^2+mu*e*(t-(1-e)/mu));
end

function s = find_zero_divide(func, min, max)
    while true
        j_max = func(max);
        j_min = func(min);
        if j_min*j_max > 0
            max = max * 2;
        else
            break;
        end
    end
    while true
        s = (min+max)/2;
        j = func(s);
        if abs(j_min)<1e-9 && abs(j_max)<1e-9
            if j_min > 0
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
