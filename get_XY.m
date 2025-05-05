function [x,y] = get_XY(t,theta,u,x0,y0,vx0,vy0,mu)
x = 1/mu*vx0*(1-exp(-mu*t))+u*cos(theta)/mu.*(t-(1-exp(-mu*t))/mu)+x0;
y = 1/mu*vy0*(1-exp(-mu*t))+u*sin(theta)/mu.*(t-(1-exp(-mu*t))/mu)+y0;
end