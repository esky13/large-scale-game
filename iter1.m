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

[op_tx,op_ty,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA_double_integrator(xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
[op_vx,op_vy] = get_vXY(t_minimal,op_thetaA,uA,vAx0,vAy0,mu);
draw_ADR(ti,to,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu,ax1,[0 0 0],true);
draw_trajectory(op_thetaA,op_thetaD,t_minimal,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
scatter(op_tx,op_ty,'magenta','filled','pentagram')
tx_lst = [];
ty_lst = [];
xA = xA0;
yA = yA0;
vAx = vAx0;
vAy = vAy0;
xD = xD0;
yD = yD0;
vDx = vDx0;
vDy = vDy0;
i = 1;
I = 50;
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

tx_f = tx1;
ty_f = ty1;
sx2 = tx1;
sy2 = ty1;
kkk = 1;
while true
    if true
        tx_lst(end+1) = tx_f;
        ty_lst(end+1) = ty_f;
    end
    if i>=I
        break
    end
    [~,~,sx1,sy1] = find_optimal_RA(xD,yD,xA,yA,uA/uD);
    tx_f = tx_f+kkk*(sx1-sx2);
    ty_f = ty_f+kkk*(sy1-sy2);

    [thetaA2,tA] = ddi_target_theta(xA0,yA0,vAx0,vAy0,uA,tx_f,ty_f,mu);
    [thetaD2,tD] = ddi_target_theta(xD0,yD0,vDx0,vDy0,uD,tx_f,ty_f,mu);
    [xA,yA] = get_XY(Te,thetaA2,uA,xA0,yA0,vAx0,vAy0,mu);
    [xD,yD] = get_XY(Te,thetaD2,uD,xD0,yD0,vDx0,vDy0,mu);
    sx2 = sx1;
    sy2 = sy1;

    tt = min(tA,tD);
    draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
    i = i+1;
end
draw_trajectory(thetaA2,thetaD2,tt,mu,uA,xA0,yA0,vAx0,vAy0,uD,xD0,yD0,vDx0,vDy0,ax1)
plot(tx_lst,ty_lst,'LineWidth',1,'Color',[0.5,0.5,0.5]);
scatter(tx_lst,ty_lst,20,[0.5,0.5,0.5],'filled');
scatter(tx,ty,50,'yellow','filled','pentagram');
sqrt((op_tx-tx_f)^2+(op_ty-ty_f)^2)

