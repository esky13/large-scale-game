function [tx,ty,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA_double_integrator(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu)
    a = find_to(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
    b = find_ti(xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
    ti = b(1);
    to = a(1);
    to1 = -1; to2 = -1;
    if length(a) == 3
        if a(3) < ti
            to = a(3); to1 = a(1); to2 = a(2);
        end
    end

    [t_minimal,op_thetaA,op_thetaD,valid] = attacker_win_strategy(to,ti,xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu);
    tx = 0;
    ty = 0;
    if ~valid
        [t_minimal,op_thetaA,op_thetaD,dis_min,tx,ty] = get_dis_minimal(to,ti,xA,yA,vAx,vAy,uA,xD,yD,vDx,vDy,uD,mu,a(1));
    end
end