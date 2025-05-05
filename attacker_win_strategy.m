function [t_reach,theta_reach,thetaD,valid] = attacker_win_strategy(tmin,tmax,x0,y0,vx0,vy0,u,xD0,yD0,vDx0,vDy0,uD,mu)
tf = tmax;
if x0*vx0+y0*vy0<0 && (vx0^2+vy0^2)/mu+x0*vx0+y0*vy0 > 0
    [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
    [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
    [dd_xc,dd_yc,dd_rc] = get_ddIso(tf,vx0,vy0,u,mu);
    dd_j = d_xc^2+d_yc^2-d_rc^2+xc*dd_xc+yc*dd_yc-rc*dd_rc;
    t1 = 0; t2 = tf;
    while dd_j > 0
        t2 = 2*t2;
        tf = t2;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
        [dd_xc,dd_yc,dd_rc] = get_ddIso(tf,vx0,vy0,u,mu);
        dd_j = d_xc^2+d_yc^2-d_rc^2+xc*dd_xc+yc*dd_yc-rc*dd_rc;
    end
    while abs(dd_j) > 1e-3
        if dd_j > 0
            t1 = tf;
        else
            t2 = tf;
        end
        tf = (t1+t2)/2;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
        [dd_xc,dd_yc,dd_rc] = get_ddIso(tf,vx0,vy0,u,mu);
        dd_j = d_xc^2+d_yc^2-d_rc^2+xc*dd_xc+yc*dd_yc-rc*dd_rc;
    end

    t_climb = tf;
    [rc, xc, yc] = getIsochrones(t_climb,x0,y0,vx0,vy0,u,mu);
    [d_xc,d_yc,d_rc] = get_dIso(t_climb,vx0,vy0,u,mu);
    d_j = xc*d_xc+yc*d_yc-rc*d_rc;
    if d_j < 0
        tf = find_zerospoints(0,Inf,x0,y0,vx0,vy0,u,mu);% tmax
    else
        [t_minimal,dis_min] = find_centerdis_localmin(0,t_climb,x0,y0,vx0,vy0,u,mu);
        [t_maximum,dis_max] = find_centerdis_localmax(t_climb,Inf,x0,y0,vx0,vy0,u,mu);
        if dis_min > 0
            t_maximum_s = t_maximum;
            while true
                [rc, xc, yc] = getIsochrones(t_maximum_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_maximum_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j > 0
                    t_maximum_s = t_maximum_s+1e-5;
                else
                    break;
                end
            end
            tf = find_zerospoints(t_maximum_s,Inf,x0,y0,vx0,vy0,u,mu);
        elseif dis_max < 0
            t_minimal_s = t_minimal;
            while true
                [rc, xc, yc] = getIsochrones(t_minimal_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_minimal_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j > 0
                    t_minimal_s = t_minimal_s-1e-5;
                else
                    break;
                end
            end
            tf = find_zerospoints(0,t_minimal_s,x0,y0,vx0,vy0,u,mu);
        else
            t_minimal_s = t_minimal;
            while true
                [rc, xc, yc] = getIsochrones(t_minimal_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_minimal_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j > 0
                    t_minimal_s = t_minimal_s-1e-5;
                else
                    break;
                end
            end
            tf1 = find_zerospoints(0,t_minimal_s,x0,y0,vx0,vy0,u,mu);
            t_minimal_s = t_minimal;
            while true
                [rc, xc, yc] = getIsochrones(t_minimal_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_minimal_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j < 0
                    t_minimal_s = t_minimal_s+1e-5;
                else
                    break;
                end
            end
            t_maximum_s = t_maximum;
            while true
                [rc, xc, yc] = getIsochrones(t_maximum_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_maximum_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j < 0
                    t_maximum_s = t_maximum_s-1e-5;
                else
                    break;
                end
            end
            tf2 = find_zerospoints(t_minimal_s,t_maximum_s,x0,y0,vx0,vy0,u,mu);
            t_maximum_s = t_maximum;
            while true
                [rc, xc, yc] = getIsochrones(t_maximum_s,x0,y0,vx0,vy0,u,mu);
                [d_xc,d_yc,d_rc] = get_dIso(t_maximum_s,vx0,vy0,u,mu);
                d_j = xc*d_xc+yc*d_yc-rc*d_rc;
                if d_j > 0
                    t_maximum_s = t_maximum_s+1e-5;
                else
                    break;
                end
            end
            tf3 = find_zerospoints(t_maximum_s,Inf,x0,y0,vx0,vy0,u,mu);
            tf = [tf1,tf2,tf3];
        end
    end
elseif (vx0^2+vy0^2)/mu+x0*vx0+y0*vy0<0
    tf = find_zerospoints(0,Inf,x0,y0,vx0,vy0,u,mu);
else
    tf = tmax;
    [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
    j = xc^2+yc^2-rc^2;
    t1 = 0; t2 = tf;
    while j > 0
        t2 = 2*t2;
        tf = t2;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        j = xc^2+yc^2-rc^2;
    end
    while abs(j) > 1e-3
        if j > 0
            t1 = tf;
        else
            t2 = tf;
        end
        tf = (t1+t2)/2;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        j = xc^2+yc^2-rc^2;
    end
end
t_reach = 0;
theta_reach = 0;
thetaD = 0;
valid = false;
max_dis = -Inf;
for t = tf
    t_reach = t;%
    if t > tmax
        continue
    elseif t < tmin
        valid = true;
    else
        [~, xc, yc] = getIsochrones(t,x0,y0,vx0,vy0,u,mu);
        thetaA = atan2(-yc,-xc);
        step = (t-tmin)/100; t_cur = t;
        while true
            [xAt,yAt] = get_XY(t_cur,thetaA,u,x0,y0,vx0,vy0,mu);
            [rDc, xDc, yDc] = getIsochrones(t_cur,xD0,yD0,vDx0,vDy0,uD,mu);
            dis_cur = (xDc-xAt)^2+(yDc-yAt)^2-rDc^2;
            if dis_cur < 0
                break
            end
            if t_cur < tmin
                valid = true;
                break
            end
            t_cur = t_cur - step;
        end
    end
    if valid
        [rDc, xDc, yDc] = getIsochrones(t,xD0,yD0,vDx0,vDy0,uD,mu);
        dis = sqrt(xDc^2+yDc^2)-rDc;
        if dis > max_dis
            max_dis = dis;
            thetaD = atan2(-yDc,-xDc);
            t_reach = t;
            [~, xAc, yAc] = getIsochrones(t,x0,y0,vx0,vy0,u,mu);
            theta_reach = atan2(-yAc,-xAc);
        end
    end
end
end

function [d_x,d_y,d_r] = get_dIso(t,vx0,vy0,u,mu)
d_x = vx0*exp(-mu*t);
d_y = vy0*exp(-mu*t);
d_r = u/mu*(1-exp(-mu*t));
end

function [dd_x,dd_y,dd_r] = get_ddIso(t,vx0,vy0,u,mu)
dd_x = -mu*vx0*exp(-mu*t);
dd_y = -mu*vy0*exp(-mu*t);
dd_r = u*exp(-mu*t);
end

function tf = find_zerospoints(t_start,t_end,x0,y0,vx0,vy0,u,mu)
tf = t_start;
beta = 0.8; eps = 1e-4;

[rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
[d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
j = xc^2+yc^2-rc^2;
d_j = xc*d_xc+yc*d_yc-rc*d_rc;

while abs(j)>eps
    t = 1;
    delta = j/d_j;
    while true
        temp_tf = tf - t*delta;
        [rc, xc, yc] = getIsochrones(temp_tf,x0,y0,vx0,vy0,u,mu);
        temp_j = xc^2+yc^2-rc^2;
        if temp_tf < t_start || temp_tf > t_end || abs(temp_j) > abs(j)
            t = t*beta;
        else
            break;
        end
    end
    tf = tf - t*delta;
    [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
    [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
    j = xc^2+yc^2-rc^2;
    d_j = xc*d_xc+yc*d_yc-rc*d_rc;
end
end

function [tf,j] = find_centerdis_localmin(t_start,t_end,x0,y0,vx0,vy0,u,mu)
    tf = t_start;
    [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
    [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
    j = xc^2+yc^2-rc^2;
    d_j = xc*d_xc+yc*d_yc-rc*d_rc;

    eps = 5*1e-4;
    alpha = 0.5; beta = 0.8;
    while abs(d_j)>eps
        t = 1;
        while true
            temp_tf = tf-alpha*t*d_j;
            [rc, xc, yc] = getIsochrones(temp_tf,x0,y0,vx0,vy0,u,mu);
            temp_j = xc^2+yc^2-rc^2;
            if temp_tf < t_start || temp_tf > t_end
                t = beta*t;
                continue
            end
            if temp_j < j
                break;
            else
                t = beta*t;
            end
        end
        tf = tf-alpha*t*d_j;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
        j = xc^2+yc^2-rc^2;
        d_j = xc*d_xc+yc*d_yc-rc*d_rc;
    end
end

function [tf,j] = find_centerdis_localmax(t_start,t_end,x0,y0,vx0,vy0,u,mu)
    tf = t_start;
    [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
    [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
    j = xc^2+yc^2-rc^2;
    d_j = xc*d_xc+yc*d_yc-rc*d_rc;

    eps = 5*1e-4;
    alpha = 0.5; beta = 0.8;
    while abs(d_j)>eps
        t = 1;
        while true
            temp_tf = tf+alpha*t*d_j;
            [rc, xc, yc] = getIsochrones(temp_tf,x0,y0,vx0,vy0,u,mu);
            temp_j = xc^2+yc^2-rc^2;
            if temp_tf < t_start || temp_tf > t_end
                t = beta*t;
                continue
            end
            if temp_j > j
                break;
            else
                t = beta*t;
            end
        end
        tf = tf+alpha*t*d_j;
        [rc, xc, yc] = getIsochrones(tf,x0,y0,vx0,vy0,u,mu);
        [d_xc,d_yc,d_rc] = get_dIso(tf,vx0,vy0,u,mu);
        j = xc^2+yc^2-rc^2;
        d_j = xc*d_xc+yc*d_yc-rc*d_rc;
    end
end
