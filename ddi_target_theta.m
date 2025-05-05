function [theta,t] = ddi_target_theta(x,y,vx,vy,u,tx,ty,mu)
    t = find_to(x,y,vx,vy,u,tx,ty,0,0,0,mu);
    t = t(1);
    [~, xc, yc] = getIsochrones(t,x,y,vx,vy,u,mu);
    theta = atan2(ty-yc,tx-xc);
end

