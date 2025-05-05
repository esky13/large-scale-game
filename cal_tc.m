function tc = cal_tc(ts,vx0,vy0,u,mu)
v02 = vx0^2+vy0^2;
tc = ts/2;
eps = 1e-5;
l = 0; r = ts;
while true
    f = u^2/mu*tc-exp(-mu*tc)*(u^2/mu^2-exp(-mu*tc)*(u^2/mu^2-v02));
    if abs(f) < eps
        break
    end
    if f > 0
        r = tc;
        tc = (tc+l)/2;
    else
        l = tc;
        tc = (tc+r)/2;
    end
end
end
