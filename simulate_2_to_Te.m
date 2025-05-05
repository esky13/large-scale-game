function [Te,bx,by] = simulate_2_to_Te(tx,ty,tx2,ty2,x0,y0,vx0,vy0,u,mu)
    Te = 0;
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
    while true
        [x,y] = get_XY(dt,theta_t,u,x,y,vx,vy,mu);
        [vx,vy] = get_vXY(dt,theta_t,u,vx,vy,mu);
        [x2,y2] = get_XY(dt,theta_t2,u,x2,y2,vx2,vy2,mu);
        [vx2,vy2] = get_vXY(dt,theta_t2,u,vx2,vy2,mu);
        Te = Te + dt;
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
        if abs(ftheta) < 1e-3 || abs(ftheta) < abs(ftheta0/100)
            bx = ibx/(v_theta2-v_theta);
            by = iby/(v_theta2-v_theta);
            return
        end
    end
end
