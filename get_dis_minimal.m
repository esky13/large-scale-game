function [t_minimal,thetaA,thetaD,dis_min,sx1,sy1] = get_dis_minimal(tmin,tmax,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu,ti1)
delta = (tmax-tmin)/100;
minimal_interval = 1e-10;
prevs1 = 0; prevs2 = 0;

while true
    [rDc, xDc, yDc] = getIsochrones(tmin,xD0,yD0,vDx0,vDy0,uD,mu);
    [rAc, xAc, yAc] = getIsochrones(tmin,xA0,yA0,vAx0,vAy0,uA,mu);
    [hasroot,x1,y1,x2,y2]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
    if hasroot
        break
    end
    tmin = tmin+delta;
end
if x1^2+y1^2<x2^2+y2^2
    dis_min = x1^2+y1^2;
    thetaA = atan2(y1-yAc,x1-xAc);
    thetaD = atan2(y1-yDc,x1-xDc);
    sx1=x1;sy1=y1;
else
    dis_min = x2^2+y2^2;
    thetaA = atan2(y2-yAc,x2-xAc);
    thetaD = atan2(y2-yDc,x2-xDc);
    sx1=x2;sy1=y2;
end
t_minimal = tmin;

while true
    [rDc, xDc, yDc] = getIsochrones(tmax,xD0,yD0,vDx0,vDy0,uD,mu);
    [rAc, xAc, yAc] = getIsochrones(tmax,xA0,yA0,vAx0,vAy0,uA,mu);
    [hasroot,x1,y1,x2,y2]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
    if hasroot
        break
    end
    tmax = tmax-delta;
end
if x1^2+y1^2<dis_min
    dis_min = x1^2+y1^2;
    thetaA = atan2(y1-yAc,x1-xAc);
    thetaD = atan2(y1-yDc,x1-xDc);
    t_minimal = tmax;
    sx1=x1;sy1=y1;
elseif x2^2+y2^2<dis_min
    dis_min = x2^2+y2^2;
    thetaA = atan2(y2-yAc,x2-xAc);
    thetaD = atan2(y2-yDc,x2-xDc);
    t_minimal = tmax;
    sx1=x2;sy1=y2;
end
if tmax < tmin
    return
end
t_step = tmin:delta:tmax;
for t = t_step
    [rDc, xDc, yDc] = getIsochrones(t,xD0,yD0,vDx0,vDy0,uD,mu);
    [rAc, xAc, yAc] = getIsochrones(t,xA0,yA0,vAx0,vAy0,uA,mu);
    [hasroot,x1,y1,x2,y2]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
    if ~hasroot
        continue
    end
    thetaA1 = atan2(y1-yAc,x1-xAc);
    thetaD1 = atan2(y1-yDc,x1-xDc);
    thetaA2 = atan2(y2-yAc,x2-xAc);
    thetaD2 = atan2(y2-yDc,x2-xDc);

    s1 = cal_dis_derivative(thetaA1,thetaD1,t,x1,y1,uA,uD,vAx0,vAy0,vDx0,vDy0,mu);
    if prevs1*s1 < 0
        left = t-delta; right = t; t_cur = left;
        s_left = prevs1;
        while right-left > minimal_interval
            t_cur = (left+right)/2;
            [rDc, xDc, yDc] = getIsochrones(t_cur,xD0,yD0,vDx0,vDy0,uD,mu);
            [rAc, xAc, yAc] = getIsochrones(t_cur,xA0,yA0,vAx0,vAy0,uA,mu);
            [~,x1,y1,~,~]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
            thetaA1 = atan2(y1-yAc,x1-xAc);
            thetaD1 = atan2(y1-yDc,x1-xDc);
            s_cur = cal_dis_derivative(thetaA1,thetaD1,t_cur,x1,y1,uA,uD,vAx0,vAy0,vDx0,vDy0,mu);
            if abs(s_cur)<1e-5
                break
            end
            if s_left*s_cur > 0
                left = t_cur;
                s_left = s_cur;
            else
                right = t_cur;
            end
        end
        if x1^2+y1^2 < dis_min
            % scatter(x1,y1,'filled')
            [F,nu,critria] = get_condition(thetaA1,thetaD1,t_cur,xA0,yA0,vAx0,vAy0,vDx0,vDy0,uA,uD,mu);
            if nu > 0 && critria > 0
                % b_intersect = judge_intersect(thetaA1,ti1,t_cur,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
                % if ~b_intersect
                    t_minimal = t_cur;
                    dis_min = x1^2+y1^2;
                    thetaA = thetaA1;
                    thetaD = thetaD1;
                    sx1=x1;sy1=y1;
                % end
            end
        end
    end
        
    s2 = cal_dis_derivative(thetaA2,thetaD2,t,x2,y2,uA,uD,vAx0,vAy0,vDx0,vDy0,mu);
    if prevs2*s2 < 0
        left = t-delta; right = t; t_cur = left;
        s_left = prevs2;
        while right-left > minimal_interval
            t_cur = (left+right)/2;
            [rDc, xDc, yDc] = getIsochrones(t_cur,xD0,yD0,vDx0,vDy0,uD,mu);
            [rAc, xAc, yAc] = getIsochrones(t_cur,xA0,yA0,vAx0,vAy0,uA,mu);
            [~,~,~,x2,y2]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
            thetaA2 = atan2(y2-yAc,x2-xAc);
            thetaD2 = atan2(y2-yDc,x2-xDc);
            s_cur = cal_dis_derivative(thetaA2,thetaD2,t_cur,x2,y2,uA,uD,vAx0,vAy0,vDx0,vDy0,mu);
            if abs(s_cur)<1e-5
                break
            end
            if s_left*s_cur > 0
                left = t_cur;
                s_left = s_cur;
            else
                right = t_cur;
            end
        end

        if x2^2+y2^2 < dis_min
            % scatter(x2,y2,'filled')
            [F,nu,critria] = get_condition(thetaA2,thetaD2,t_cur,xA0,yA0,vAx0,vAy0,vDx0,vDy0,uA,uD,mu);
            if nu > 0 && critria > 0
                % b_intersect = judge_intersect(thetaA2,ti1,t_cur,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
                % if ~b_intersect
                    t_minimal = t_cur;
                    dis_min = x2^2+y2^2;
                    thetaA = thetaA2;
                    thetaD = thetaD2;
                    sx1=x2;sy1=y2;
                % end
            end
        end
    end

    prevs1 = s1;prevs2 = s2;
end
end

function s = cal_dis_derivative(thetaA,thetaD,t,xf,yf,uA,uD,vAx0,vAy0,vDx0,vDy0,mu)
    [vAx,vAy] = get_vXY(t,thetaA,uA,vAx0,vAy0,mu);
    [vDx,vDy] = get_vXY(t,thetaD,uD,vDx0,vDy0,mu);
    phi = atan2(yf,xf);
    s = (sin(thetaA-phi)*(vDx*cos(thetaD)+vDy*sin(thetaD)) ...
        -sin(thetaD-phi)*(vAx*cos(thetaA)+vAy*sin(thetaA)))/sin(thetaA-thetaD);
end

function [F,nu,critria] = get_condition(thetaA,thetaD,tf,xA0,yA0,vAx0,vAy0,vDx0,vDy0,uA,uD,mu)
[xAf,yAf] = get_XY(tf,thetaA,uA,xA0,yA0,vAx0,vAy0,mu);
phi = atan2(yAf,xAf);
nu = sin(thetaA-phi)/sin(thetaA-thetaD);
gamma = sqrt(1+nu^2-2*nu*cos(thetaD-phi));
critria = sin(thetaA-phi)/sin(thetaD-phi);
F = nu*(vDx0*cos(thetaD)+vDy0*sin(thetaD))-gamma*(vAx0*cos(thetaA)+vAy0*sin(thetaA))+1/mu*(exp(mu*tf)-1)*(nu*uD-gamma*uA);
end