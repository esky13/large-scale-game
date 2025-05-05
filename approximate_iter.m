% close all
mu = 1;
xA0 = -14.8229; yA0 = 10.3939;
xD0 = -0.9002; yD0 = 2.6429;
vAx0 = 0.4186; vAy0 = 0.1537; 
vDx0 = -0.538; vDy0 = 0.7667;

% xA0 = -20; yA0 = 10;
% xD0 = -1; yD0 = 2;
% vAx0 = 1; vAy0 = 0; 
% vDx0 = 0; vDy0 = 2;

uA = 1; uD = 2;

ax1 = figure;hold on 
axis equal
xlabel('x');
ylabel('y');
set(gca,'FontName','Times New Roman','FontSize',20)
scatter(0,0,'green','filled')
scatter(xA0,yA0,'red','filled')
scatter(xD0,yD0,'blue','filled')

t = 0;
dt = 0.001;
xA = xA0;
yA = yA0;
vAx = vAx0;
vAy = vAy0;
xD = xD0;
yD = yD0;
vDx = vDx0;
vDy = vDy0;

[op_tx,op_ty,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA_double_integrator(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
[op_vx,op_vy] = get_vXY(t_minimal,op_thetaA,uA,vAx0,vAy0,mu);
draw_ADR(ti,to,xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu,ax1,[0 0 0],true);
draw_trajectory(op_thetaA,op_thetaD,t_minimal,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
scatter(op_tx,op_ty,'magenta','filled','pentagram')
tx_lst = [];
ty_lst = [];
xA_lst = [];
yA_lst = [];
xD_lst = [];
yD_lst = [];
thetaA_lst = [];
thetaD_lst = [];
tx2_lst = [];
ty2_lst = [];
i = 1;
I = 10;
[~,~,tx0,ty0]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
tx = tx0;
ty = ty0;
[thetaA1,~] = ddi_target_theta(xA0,yA0,vAx0,vAy0,uA,tx,ty,mu);
[thetaD1,~] = ddi_target_theta(xD0,yD0,vDx0,vDy0,uD,tx,ty,mu);
% simulate to Te
[TeA,~,~] = simulate_to_Te(tx,ty,xA0,yA0,vAx0,vAy0,uA,mu);
[TeD,~,~] = simulate_to_Te(tx,ty,xD0,yD0,vDx0,vDy0,uD,mu);
Te = min(TeA,TeD);
[xA1,yA1] = get_XY(Te,thetaA1,uA,xA0,yA0,vAx0,vAy0,mu);
scatter(xA1,yA1,20,[0,0.8,0.3],'filled');
[xD1,yD1] = get_XY(Te,thetaD1,uD,xD0,yD0,vDx0,vDy0,mu);
scatter(xD1,yD1,20,[0.5,0.7,0],'filled');

[~,~,tx1,ty1]=find_optimal_RA(xD1,yD1,xA1,yA1,uA/uD);
[thetaA2,~] = ddi_target_theta(xA0,yA0,vAx0,vAy0,uA,tx1,ty1,mu);
[thetaD2,~] = ddi_target_theta(xD0,yD0,vDx0,vDy0,uD,tx1,ty1,mu);
[xA,yA] = get_XY(Te,thetaA2,uA,xA0,yA0,vAx0,vAy0,mu);
scatter(xA,yA,20,[0,0.8,0.3],'filled');
[xD,yD] = get_XY(Te,thetaD2,uD,xD0,yD0,vDx0,vDy0,mu);
scatter(xD,yD,20,[0.5,0.7,0],'filled');

[vAx1,vAy1] = get_vXY(Te,thetaA1,uA,vAx0,vAy0,mu);
[vAx2,vAy2] = get_vXY(Te,thetaA2,uA,vAx0,vAy0,mu);
delta_theta_trueA = atan2(vAy2,vAx2)-atan2(vAy1,vAx1);

[vDx1,vDy1] = get_vXY(Te,thetaD1,uD,vDx0,vDy0,mu);
[vDx2,vDy2] = get_vXY(Te,thetaD2,uD,vDx0,vDy0,mu);
delta_theta_trueD = atan2(vDy2,vDx2)-atan2(vDy1,vDx1);

alpha_xA = (xA-xA1)/delta_theta_trueA+yA;
alpha_yA = xA-(yA-yA1)/delta_theta_trueA;
alpha_xD = (xD-xD1)/delta_theta_trueD+yD;
alpha_yD = xD-(yD-yD1)/delta_theta_trueD;

t_fA = sqrt((tx-xA0)^2+(ty-yA0)^2)/(uA/mu);
t_fD = sqrt((tx-xD0)^2+(ty-yD0)^2)/(uD/mu);
t_f = min(t_fA,t_fD);
simulate_2_to_tf(tx0,ty0,tx1,ty1,xA0,yA0,vAx0,vAy0,uA,mu,t_f,Te);

xA2 = xA;
yA2 = yA;
xD2 = xD;
yD2 = yD;
% find terminal
while true
    if true
        tx_lst(end+1) = tx;
        ty_lst(end+1) = ty;
        thetaA_lst(end+1) = thetaA2;
        thetaD_lst(end+1) = thetaD2;
    end

    if i>=I
        break
    end
    [~,~,tx,ty]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
    scatter(tx,ty,50,'yellow','filled','pentagram')
    [thetaA2,tA] = ddi_target_theta(xA0,yA0,vAx0,vAy0,uA,tx,ty,mu);
    [thetaD2,tD] = ddi_target_theta(xD0,yD0,vDx0,vDy0,uD,tx,ty,mu);
    tt = min(tA,tD);
    draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)

    [vAx,vAy] = get_vXY(Te,thetaA2,uA,vAx0,vAy0,mu);
    delta_theta_trueA = atan2(vAy1,vAx1)-atan2(vAy,vAx);

    delta_theta_xA = (tx1-tx)/(-ty1+alpha_xA);
    delta_theta_yA = (ty1-ty)/(tx1-alpha_yA); 
    % avg_delta_theta_A = (delta_theta_xA+delta_theta_yA)/2;
    delta_Ax = delta_theta_xA*(-yA2+alpha_xA);
    delta_Ay = delta_theta_yA*(xA2-alpha_yA);
    if abs(-ty1+alpha_xA)<1e-1
        delta_Ax = 0;
    end
    if abs(tx1-alpha_yA)<1e-1
        delta_Ay = 0;
    end

    est_Ax = xA2-delta_Ax;
    est_Ay = yA2-delta_Ay;
    % scatter(est_Ax,est_Ay,50,[0,0.8,0.3],"square",'filled','+');
    %

    [vDx,vDy] = get_vXY(Te,thetaD2,uD,vDx0,vDy0,mu);
    delta_theta_trueD = atan2(vDy1,vDx1)-atan2(vDy,vDx);

    delta_theta_xD = (tx1-tx)/(-ty1+alpha_xD);
    delta_theta_yD = (ty1-ty)/(tx1-alpha_yD); 
    % avg_delta_theta_D = (delta_theta_xD+delta_theta_yD)/2;
    delta_Dx = delta_theta_xD*(-yD2+alpha_xD);
    delta_Dy = delta_theta_yD*(xD2-alpha_yD);
    if abs(-ty1+alpha_xD)<1e-1
        delta_Dx = 0;
    end
    if abs(tx1-alpha_yD)<1e-1
        delta_Dy = 0;
    end

    est_Dx = xD2-delta_Dx;
    est_Dy = yD2-delta_Dy;
    % scatter(est_Dx,est_Dy,50,[0.5,0.7,0],"square",'filled','+');

    xA = est_Ax;
    yA = est_Ay;
    xD = est_Dx;
    yD = est_Dy;

    scatter(xA,yA,20,[0,0.8,0.3],'filled');
    scatter(xD,yD,20,[0.5,0.7,0],'filled');

    i = i+1;
end
draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
plot(tx_lst,ty_lst,'LineWidth',1,'Color',[0.5,0.5,0.5]);
scatter(tx_lst,ty_lst,20,[0.5,0.5,0.5],'filled');
scatter(tx,ty,50,'yellow','filled','pentagram');
sqrt((op_tx-tx)^2+(op_ty-ty)^2)


function [ax,ay]=simulate_2_to_tf(tx,ty,tx2,ty2,x0,y0,vx0,vy0,u,mu,t_f,t_E)
    t = 0;
    dt = 0.1;
    theta0 = atan2(vy0,vx0);
    theta_t = ddi_target_theta(x0,y0,vx0,vy0,u,tx,ty,mu);
    theta_t2 = ddi_target_theta(x0,y0,vx0,vy0,u,tx2,ty2,mu);
    ftheta0 = ddi_ftheta(u,sqrt(vx0^2+vy0^2),theta_t,theta0);
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;
    x2 = x0;
    y2 = y0;
    vx2 = vx0;
    vy2 = vy0;
    ibx = 0;
    iby = 0;
    iax = 0;
    iay = 0;
    while true
        [x,y] = get_XY(dt,theta_t,u,x,y,vx,vy,mu);
        [vx,vy] = get_vXY(dt,theta_t,u,vx,vy,mu);
        [x2,y2] = get_XY(dt,theta_t2,u,x2,y2,vx2,vy2,mu);
        [vx2,vy2] = get_vXY(dt,theta_t2,u,vx2,vy2,mu);
        t = t + dt;
        v = sqrt(vx^2+vy^2);
        v2 = sqrt(vx2^2+vy2^2);
        v_theta = atan2(vy,vx);
        v_theta2 = atan2(vy2,vx2);
        ftheta = ddi_ftheta(u,v,theta_t,v_theta);
        ftheta2 = ddi_ftheta(u,v2,theta_t2,v_theta2);
        fv = ddi_fv(u,mu,v,theta_t,v_theta);
        fv2 = ddi_fv(u,mu,v2,theta_t2,v_theta2);
        ibx = ibx + dt*cos(v_theta2)*(v2-v);
        iby = iby + dt*sin(v_theta2)*(v2-v);
        iax = iax + dt*(ftheta2-ftheta)*y2;
        iay = iay + dt*(ftheta2-ftheta)*x2;
        if t>t_f
            ax = (iax+ibx)/(v_theta2-v_theta);
            ay = (iay+iby)/(v_theta2-v_theta);
            return
        end
        if abs(t-t_E)<=dt
            ax = (iax+ibx)/(v_theta2-v_theta);
            ay = (iay+iby)/(v_theta2-v_theta);
        end
    end
end
