function [tx,ty,thetaA2,thetaD2] = approximate_opti_control(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu)
i = 1;
I = 6;
[~,~,tx0,ty0]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
tx = tx0;
ty = ty0;
t_fA = sqrt((tx-xA)^2+(ty-yA)^2)/(uA/mu);
t_fD = sqrt((tx-xD)^2+(ty-yD)^2)/(uD/mu);
t_f = min(t_fA,t_fD);
[thetaA1,~] = ddi_target_theta(xA,yA,vAx,vAy,uA,tx,ty,mu);
[thetaD1,~] = ddi_target_theta(xD,yD,vDx,vDy,uD,tx,ty,mu);
% simulate to Te
[TeA,~,~] = simulate_to_Te(tx,ty,xA,yA,vAx,vAy,uA,mu);
[TeD,~,~] = simulate_to_Te(tx,ty,xD,yD,vDx,vDy,uD,mu);
Te = min(TeA,TeD);
if Te>t_f
    thetaA2 = thetaA1;
    thetaD2 = thetaD1;
    return 
end
[xA1,yA1] = get_XY(Te,thetaA1,uA,xA,yA,vAx,vAy,mu);
[xD1,yD1] = get_XY(Te,thetaD1,uD,xD,yD,vDx,vDy,mu);

[~,~,tx1,ty1]=find_optimal_RA(xD1,yD1,xA1,yA1,uA/uD);
[thetaA2,~] = ddi_target_theta(xA,yA,vAx,vAy,uA,tx1,ty1,mu);
[thetaD2,~] = ddi_target_theta(xD,yD,vDx,vDy,uD,tx1,ty1,mu);
[~,bxA,byA] = simulate_2_to_Te(tx,ty,tx1,ty1,xA,yA,vAx,vAy,uA,mu);
[~,bxD,byD] = simulate_2_to_Te(tx,ty,tx1,ty1,xD,yD,vDx,vDy,uD,mu);
[xA,yA] = get_XY(Te,thetaA2,uA,xA,yA,vAx,vAy,mu);
[xD,yD] = get_XY(Te,thetaD2,uD,xD,yD,vDx,vDy,mu);

[vAx1,vAy1] = get_vXY(Te,thetaA1,uA,vAx,vAy,mu);
[vAx2,vAy2] = get_vXY(Te,thetaA2,uA,vAx,vAy,mu);
delta_theta_trueA = atan2(vAy2,vAx2)-atan2(vAy1,vAx1);

[vDx1,vDy1] = get_vXY(Te,thetaD1,uD,vDx,vDy,mu);
[vDx2,vDy2] = get_vXY(Te,thetaD2,uD,vDx,vDy,mu);
delta_theta_trueD = atan2(vDy2,vDx2)-atan2(vDy1,vDx1);

alpha_xA = ((xA-xA1)/delta_theta_trueA-bxA+yA);
alpha_yA = (xA+byA-(yA-yA1)/delta_theta_trueA);
alpha_xD = ((xD-xD1)/delta_theta_trueD-bxD+yD);
alpha_yD = (xD+byD-(yD-yD1)/delta_theta_trueD);

xA2 = xA;
yA2 = yA;
xD2 = xD;
yD2 = yD;

while true

    if i>=I
        break
    end
    [~,~,tx,ty]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
    [thetaA2,tA] = ddi_target_theta(xA,yA,vAx,vAy,uA,tx,ty,mu);
    [thetaD2,tD] = ddi_target_theta(xD,yD,vDx,vDy,uD,tx,ty,mu);

    delta_theta_xA = (tx1-tx)/(-ty1+alpha_xA+bxA);
    delta_theta_yA = (ty1-ty)/(tx1-alpha_yA+byA); 
    % avg_delta_theta_A = (delta_theta_xA+delta_theta_yA)/2;
    delta_Ax = delta_theta_xA*(-yA2+alpha_xA+bxA);
    delta_Ay = delta_theta_yA*(xA2-alpha_yA+byA);
    if abs(-ty1+alpha_xA+bxA)<1e-2
        delta_Ax = 0;
    end
    if abs(tx1-alpha_yA+byA)<1e-2
        delta_Ay = 0;
    end

    est_Ax = xA2-delta_Ax;
    est_Ay = yA2-delta_Ay;

    delta_theta_xD = (tx1-tx)/(-ty1+alpha_xD+bxD);
    delta_theta_yD = (ty1-ty)/(tx1-alpha_yD+byD); 
    % avg_delta_theta_D = (delta_theta_xD+delta_theta_yD)/2;
    delta_Dx = delta_theta_xD*(-yD2+alpha_xD+bxD);
    delta_Dy = delta_theta_yD*(xD2-alpha_yD+byD);
    if abs(-ty1+alpha_xD+bxD)<1e-2
        delta_Dx = 0;
    end
    if abs(tx1-alpha_yD+byD)<1e-2
        delta_Dy = 0;
    end

    est_Dx = xD2-delta_Dx;
    est_Dy = yD2-delta_Dy;

    xA = est_Ax;
    yA = est_Ay;
    xD = est_Dx;
    yD = est_Dy;

    i = i+1;
end

end

