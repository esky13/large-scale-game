% close all
mu = 1;
xA0 = -25; yA0 = 10;
xD0 = -1; yD0 = 2;
vAx0 = 0; vAy0 = 1; 
vDx0 = 1; vDy0 = 1;

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
xA = xA0;
yA = yA0;
vAx = vAx0;
vAy = vAy0;
xD = xD0;
yD = yD0;
vDx = vDx0;
vDy = vDy0;
i = 1;
I = 6;
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
alpha_xA = ((tx1-tx0)/delta_theta_trueA+ty1);
alpha_yA = (tx1-(ty1-ty0)/delta_theta_trueA);

[vDx1,vDy1] = get_vXY(Te,thetaD1,uD,vDx0,vDy0,mu);
[vDx2,vDy2] = get_vXY(Te,thetaD2,uD,vDx0,vDy0,mu);
delta_theta_trueD = atan2(vDy2,vDx2)-atan2(vDy1,vDx1);
alpha_xD = ((tx1-tx0)/delta_theta_trueD+ty1);
alpha_yD = (tx1-(ty1-ty0)/delta_theta_trueD);

alpha_xA_Te = ((xA-xA1)/delta_theta_trueA+yA);
alpha_yA_Te = (xA-(yA-yA1)/delta_theta_trueA);
alpha_xD_Te = ((xD-xD1)/delta_theta_trueD+yD);
alpha_yD_Te = (xD-(yD-yD1)/delta_theta_trueD);

xA2 = xA;
yA2 = yA;
xD2 = xD;
yD2 = yD;
tx_f = tx1;
ty_f = ty1;
% find terminal
while true
    if true
        tx_lst(end+1) = tx_f;
        ty_lst(end+1) = ty_f;
        thetaA_lst(end+1) = thetaA2;
        thetaD_lst(end+1) = thetaD2;
    end

    if i>=I
        break
    end
    [~,~,tx,ty]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
    [thetaA2,tA] = ddi_target_theta(xA0,yA0,vAx0,vAy0,uA,tx,ty,mu);
    [thetaD2,tD] = ddi_target_theta(xD0,yD0,vDx0,vDy0,uD,tx,ty,mu);
    
    % [TeA,~,~] = simulate_to_Te(tx,ty,xA0,yA0,vAx0,vAy0,uA,mu);
    % [TeD,~,~] = simulate_to_Te(tx,ty,xD0,yD0,vDx0,vDy0,uD,mu);
    % Te = min(TeA,TeD);

    [vAx,vAy] = get_vXY(Te,thetaA2,uA,vAx0,vAy0,mu);
    delta_theta_trueA = atan2(vAy1,vAx1)-atan2(vAy,vAx);

    delta_theta_xA = (tx1-tx)/(-ty1+alpha_xA);
    delta_theta_yA = (ty1-ty)/(tx1-alpha_yA); 
    avg_delta_theta_A = (delta_theta_xA+delta_theta_yA)/2;

    delta_Ax = avg_delta_theta_A*(-yA2+alpha_xA_Te);
    delta_Ay = avg_delta_theta_A*(xA2-alpha_yA_Te);

    est_Ax = xA2-delta_Ax;
    est_Ay = yA2-delta_Ay;
    scatter(est_Ax,est_Ay,50,[0,0.8,0.3],"square",'filled','+');
    %
    [vDx,vDy] = get_vXY(Te,thetaD2,uD,vDx0,vDy0,mu);
    delta_theta_trueD = atan2(vDy1,vDx1)-atan2(vDy,vDx);

    delta_theta_xD = (tx1-tx)/(-ty1+alpha_xD);
    delta_theta_yD = (ty1-ty)/(tx1-alpha_yD); 
    avg_delta_theta_D = (delta_theta_xD+delta_theta_yD)/2;

    delta_Dx = avg_delta_theta_D*(-yD2+alpha_xD_Te);
    delta_Dy = avg_delta_theta_D*(xD2-alpha_yD_Te);

    est_Dx = xD2-delta_Dx;
    est_Dy = yD2-delta_Dy;

    scatter(est_Dx,est_Dy,50,[0.5,0.7,0],"square",'filled','+');
    
    xA = est_Ax;
    yA = est_Ay;
    xD = est_Dx;
    yD = est_Dy;

    tt = min(tA,tD);
    draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
    i = i+1;
end
draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
plot(tx_lst,ty_lst,'LineWidth',1,'Color',[0.5,0.5,0.5]);
scatter(tx_lst,ty_lst,20,[0.5,0.5,0.5],'filled');
scatter(tx,ty,50,'yellow','filled','pentagram');
sqrt((op_tx-tx)^2+(op_ty-ty)^2)
figure, hold on
plot(thetaA_lst,'red');
yline(op_thetaA,'red');
plot(thetaD_lst,'blue');
yline(op_thetaD,'blue');

% % go with simple motion 
% while true
%     [thetaD,thetaA,tx,ty]=find_optimal_RA(xD,yD,xA,yA,uA/uD);
%     tD = find_to(xD,yD,vDx,vDy,uD,tx,ty,0,0,0,mu);
%     tA = find_to(xA,yA,vAx,vAy,uA,tx,ty,0,0,0,mu);
%     [rDc, xDc, yDc] = getIsochrones(tD(1),xD,yD,vDx,vDy,uD,mu);
%     [rAc, xAc, yAc] = getIsochrones(tA(1),xA,yA,vAx,vAy,uA,mu);
%     thetaA2 = atan2(ty-yAc,tx-xAc);
%     thetaD2 = atan2(ty-yDc,tx-xDc);
%     [tx2,ty2,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA_double_integrator(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
%     if mod(i,10) == 0
%         tx_lst(end+1) = tx;
%         ty_lst(end+1) = ty;
%         tx2_lst(end+1) = tx2;
%         ty2_lst(end+1) = ty2;
%         xA_lst(end+1) = xA;
%         yA_lst(end+1) = yA;
%         xD_lst(end+1) = xD;
%         yD_lst(end+1) = yD;
%         thetaA_lst(end+1) = thetaA2;
%         thetaD_lst(end+1) = thetaD2;
%     end
%     [xA,yA,vAx,vAy] = double_integrator_move(xA,yA,vAx,vAy,dt,mu,thetaA2,uA);
%     [xD,yD,vDx,vDy] = double_integrator_move(xD,yD,vDx,vDy,dt,mu,thetaD2,uD);
%     if sqrt(xA^2+yA^2)<0.5 || sqrt((xA-xD)^2+(yA-yD)^2)<0.05
%         break
%     end
%     t = t + dt;
%     i = i+1;
% end
% plot(tx_lst,ty_lst,'LineWidth',1,'Color',[0.5,0.5,0.5]);
% % scatter(tx_lst,ty_lst,20,[0.5,0.5,0.5],'filled');
% plot(tx2_lst,ty2_lst,'LineWidth',1,'Color',[1,0,1]);
% % scatter(tx2_lst,ty2_lst,20,[1,0,1],'filled','pentagram');
% plot(xA_lst,yA_lst,'LineWidth',1,'Color',[0.6,0.3,0.3]);
% plot(xD_lst,yD_lst,'LineWidth',1,'Color',[0.3,0.3,0.6]);

% [tx,ty,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA_double_integrator(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
% draw_ADR(ti,to,xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu,ax1,[0 0 0],false);
% draw_trajectory(op_thetaA,op_thetaD,t_minimal,mu,uA,xA,yA,vAx,vAy,uD,xD,yD,vDx,vDy,ax1)

% figure, hold on
% plot(thetaA_lst,'red');
% yline(op_thetaA,'red');
% plot(thetaD_lst,'blue');
% yline(op_thetaD,'blue');
% figure
% plot((tx_lst-tx2_lst).^2+(ty_lst-ty2_lst).^2)

function [nx,ny,nvx,nvy] = double_integrator_move(x,y,vx,vy,dt,mu,theta,u)
    nvx = vx + dt*(u*cos(theta)-mu*vx);
    nvy = vy + dt*(u*sin(theta)-mu*vy);
    nx = x + dt*vx;
    ny = y + dt*vy;
end

function [thetaP,thetaE,tx,ty] = find_optimal_RA(xP,yP,xE,yE,alpha)
    xc = 1/(1-alpha^2)*(xE-alpha^2*xP);
    yc = 1/(1-alpha^2)*(yE-alpha^2*yP);
    rc = alpha/(1-alpha^2)*(sqrt((xE-xP)^2+(yE-yP)^2));
    tx = xc-cos(atan2(yc,xc))*rc;
    ty = yc-sin(atan2(yc,xc))*rc;
    thetaP = atan2(ty-yP,tx-xP);
    thetaE = atan2(ty-yE,tx-xE);
end

