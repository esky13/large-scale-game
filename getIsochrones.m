function [rc,xc,yc] = getIsochrones(t,x0,y0,vx0,vy0,u,mu)
    rc = u/mu*(t+exp(-mu*t)/mu-1/mu);
    xc = x0+vx0/mu*(1-exp(-mu*t));
    yc = y0+vy0/mu*(1-exp(-mu*t));
end