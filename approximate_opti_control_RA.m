function [tx,ty,thetaA2,thetaD2] = approximate_opti_control_RA(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu)
i = 1;
I = 6;
RA_ddi_problem = ddi_problem(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
sm_solver = simple_motion_solver();

% 第一次规划
[~,~,tx0,ty0] = sm_solver.find_optimal_RA(xD,yD,xA,yA,uA/uD);
tx = tx0;
ty = ty0;
t_fA = sqrt((tx-xA)^2+(ty-yA)^2)/(uA/mu);
t_fD = sqrt((tx-xD)^2+(ty-yD)^2)/(uD/mu);
t_f = min(t_fA,t_fD);
[thetaA1,~] = RA_ddi_problem.ddi_model_A.ddi_target_theta(tx,ty);
[thetaD1,~] = RA_ddi_problem.ddi_model_D.ddi_target_theta(tx,ty);
% simulate to Te
TeA = RA_ddi_problem.ddi_model_A.simulate_to_Te(tx,ty);
TeD = RA_ddi_problem.ddi_model_D.simulate_to_Te(tx,ty);
Te = min(TeA,TeD);
if Te>t_f
    thetaA2 = thetaA1;
    thetaD2 = thetaD1;
    return 
end
[xA1,yA1] = RA_ddi_problem.ddi_model_A.get_XY_from0(Te,thetaA1);
[xD1,yD1] = RA_ddi_problem.ddi_model_D.get_XY_from0(Te,thetaD1);

% 第二次规划
[~,~,tx1,ty1] = sm_solver.find_optimal_RA(xD1,yD1,xA1,yA1,uA/uD);
[thetaA2,~] = RA_ddi_problem.ddi_model_A.ddi_target_theta(tx1,ty1);
[thetaD2,~] = RA_ddi_problem.ddi_model_D.ddi_target_theta(tx1,ty1);
% Te 时刻状态
[xA,yA] = RA_ddi_problem.ddi_model_A.get_XY_from0(Te,thetaA2);
[xD,yD] = RA_ddi_problem.ddi_model_D.get_XY_from0(Te,thetaD2);
[vAx1,vAy1] = RA_ddi_problem.ddi_model_A.get_vXY_from0(Te,thetaA1);
[vAx2,vAy2] = RA_ddi_problem.ddi_model_A.get_vXY_from0(Te,thetaA2);
delta_theta_trueA = atan2(vAy2,vAx2)-atan2(vAy1,vAx1);

[vDx1,vDy1] = RA_ddi_problem.ddi_model_D.get_vXY_from0(Te,thetaD1);
[vDx2,vDy2] = RA_ddi_problem.ddi_model_D.get_vXY_from0(Te,thetaD2);
delta_theta_trueD = atan2(vDy2,vDx2)-atan2(vDy1,vDx1);

% alpha
alpha_xA = ((xA-xA1)/delta_theta_trueA+yA);
alpha_yA = (xA-(yA-yA1)/delta_theta_trueA);
alpha_xD = ((xD-xD1)/delta_theta_trueD+yD);
alpha_yD = (xD-(yD-yD1)/delta_theta_trueD);

xA2 = xA;
yA2 = yA;
xD2 = xD;
yD2 = yD;

while true
    if i>=I
        break
    end
    [~,~,tx,ty] = sm_solver.find_optimal_RA(xD,yD,xA,yA,uA/uD);
    [thetaA2,~] = RA_ddi_problem.ddi_model_A.ddi_target_theta(tx,ty);
    [thetaD2,~] = RA_ddi_problem.ddi_model_D.ddi_target_theta(tx,ty);

    delta_theta_xA = (tx1-tx)/(-ty1+alpha_xA);
    delta_theta_yA = (ty1-ty)/(tx1-alpha_yA); 
    avg_delta_theta_A = (delta_theta_xA+delta_theta_yA)/2;
    delta_Ax = avg_delta_theta_A*(-yA2+alpha_xA);
    delta_Ay = avg_delta_theta_A*(xA2-alpha_yA);
    if abs(-ty1+alpha_xA)<1e-2
        delta_Ax = 0;
    end
    if abs(tx1-alpha_yA)<1e-2
        delta_Ay = 0;
    end

    est_Ax = xA2-delta_Ax;
    est_Ay = yA2-delta_Ay;

    delta_theta_xD = (tx1-tx)/(-ty1+alpha_xD);
    delta_theta_yD = (ty1-ty)/(tx1-alpha_yD); 
    avg_delta_theta_D = (delta_theta_xD+delta_theta_yD)/2;
    delta_Dx = avg_delta_theta_D*(-yD2+alpha_xD);
    delta_Dy = avg_delta_theta_D*(xD2-alpha_yD);
    if abs(-ty1+alpha_xD)<1e-2
        delta_Dx = 0;
    end
    if abs(tx1-alpha_yD)<1e-2
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

