function [vx,vy] = get_vXY(t,theta,u,vx0,vy0,mu)
vx = vx0*exp(-mu*t)+u*cos(theta)/mu*(1-exp(-mu*t));
vy = vy0*exp(-mu*t)+u*sin(theta)/mu*(1-exp(-mu*t));
end