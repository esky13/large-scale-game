function draw_trajectory(thetaA,thetaD,tf,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax)
    t = 0:tf/100:tf;
    sA = zeros([length(t),2]);
    sD = zeros([length(t),2]);
    [sA(:,1),sA(:,2)] = get_XY(t,xA0,yA0,uA,mu,thetaA,vAx0,vAy0);
    [sD(:,1),sD(:,2)] = get_XY(t,xD0,yD0,uD,mu,thetaD,vDx0,vDy0);
    figure(ax);
    plot(sA(:,1),sA(:,2),'r','LineWidth',1);
    plot(sD(:,1),sD(:,2),'b','LineWidth',1);
end

function [x,y] = get_XY(t,x0,y0,u,mu,theta,vx0,vy0)
    a = 1-exp(-mu*t);
    x = x0 + (vx0/mu-u*cos(theta)/mu^2)*a + u*cos(theta)/mu*t;
    y = y0 + (vy0/mu-u*sin(theta)/mu^2)*a + u*sin(theta)/mu*t;
end