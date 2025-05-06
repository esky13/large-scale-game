close all
mu = 1;
xA0 = -15; yA0 = 10;
xD0 = -1; yD0 = 2;
vAx0 = 0; vAy0 = 1; 
vDx0 = 1; vDy0 = 1;

uA = 1; uD = 2;
alpha = uA/uD;

ax1 = figure;hold on 
axis equal
xlabel('x');
ylabel('y');
set(gca,'FontName','Times New Roman','FontSize',20)
scatter(0,0,'green','filled')
scatter(xA0,yA0,'red','filled')
scatter(xD0,yD0,'blue','filled')

RA_ddi_problem = ddi_problem(xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
sm_solver = simple_motion_solver();

t = 0;
i = 0;
dt = 0.001;
xA = xA0;
yA = yA0;
vAx = vAx0;
vAy = vAy0;
xD = xD0;
yD = yD0;
vDx = vDx0;
vDy = vDy0;
traj_xA_lst = [];
traj_yA_lst = [];
traj_xD_lst = [];
traj_yD_lst = [];

op_tx_lst = [];
op_ty_lst = [];
approx_tx_lst = [];
approx_ty_lst = [];
simple_tx_lst = [];
simple_ty_lst = [];

% strategy: 1-optimal, 2-large-scale-approximate, 3-simple motion
strategyA = 1
strategyD = 1
while true
    if sqrt(xA^2+yA^2)<1 || sqrt((xA-xD)^2+(yA-yD)^2)<1
        break
    end
    % [tx,ty,approx_thetaA,approx_thetaD] = approximate_opti_control_RA(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
    [op_tx,op_ty,op_thetaA,op_thetaD,t_minimal,ti,to] = RA_ddi_problem.find_optimal_RA();
    [simple_thetaD,simple_thetaA,simple_tx,simple_ty] = sm_solver.find_optimal_RA(xD,yD,xA,yA,alpha);
    if strategyA == 1
        thetaA = op_thetaA;
    elseif strategyA == 2
        thetaA = approx_thetaA;
    else
        thetaA = simple_thetaA;
    end
    if strategyD == 1
        thetaD = op_thetaD;
    elseif strategyD == 2
        thetaD = approx_thetaD;
    else
        thetaD = simple_thetaD;
    end
    RA_ddi_problem = RA_ddi_problem.step(dt,thetaA,thetaD);
    [xA,yA,vAx,vAy,xD,yD,vDx,vDy] = RA_ddi_problem.get_state();

    traj_xA_lst(end+1) = xA;
    traj_yA_lst(end+1) = yA;
    traj_xD_lst(end+1) = xD;
    traj_yD_lst(end+1) = yD;

    op_tx_lst(end+1) = op_tx;
    op_ty_lst(end+1) = op_ty;
    % approx_tx_lst(end+1) = tx;
    % approx_ty_lst(end+1) = ty;
    simple_tx_lst(end+1) = simple_tx;
    simple_ty_lst(end+1) = simple_ty;
    if mod(i,100) == 0
        % t
        % scatter(tx,ty,50,'black','filled','pentagram')
        scatter(op_tx,op_ty,50,'magenta','filled','pentagram')
        scatter(simple_tx,simple_ty,50,[0.55,0.25,0],'filled','pentagram')
        % scatter(xA,yA,'red','filled')
    end
    i = i+1;
    t = t+dt
end
plot(traj_xA_lst,traj_yA_lst,'red')
plot(traj_xD_lst,traj_yD_lst,'blue')

% figure, hold on
% plot(sqrt((op_tx_lst-approx_tx_lst).^2+(op_ty_lst-approx_ty_lst).^2),'black')
% plot(sqrt((op_tx_lst-simple_tx_lst).^2+(op_ty_lst-simple_ty_lst).^2),'Color',[0.55,0.25,0])
% cost=sqrt(xA^2+yA^2)