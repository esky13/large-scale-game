function [thetaP,thetaE,tx,ty] = find_optimal_RA(xP,yP,xE,yE,alpha)
    xc = 1/(1-alpha^2)*(xE-alpha^2*xP);
    yc = 1/(1-alpha^2)*(yE-alpha^2*yP);
    rc = alpha/(1-alpha^2)*(sqrt((xE-xP)^2+(yE-yP)^2));
    tx = xc-cos(atan2(yc,xc))*rc;
    ty = yc-sin(atan2(yc,xc))*rc;
    thetaP = atan2(ty-yP,tx-xP);
    thetaE = atan2(ty-yE,tx-xE);
end