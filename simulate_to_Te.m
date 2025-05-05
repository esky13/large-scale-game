function [Te,avg_x,avg_y] = simulate_to_Te(tx,ty,x,y,vx,vy,u,mu)
    Te = 0;
    dt = 0.1;
    theta0 = atan2(vy,vx);
    theta_t = ddi_target_theta(x,y,vx,vy,u,tx,ty,mu);
    ftheta0 = ddi_ftheta(u,sqrt(vx^2+vy^2),theta_t,theta0);
    sum_x = 0;
    sum_y = 0;
    i = 0;
    while true
        i = i+1;
        sum_x = sum_x + x;
        sum_y = sum_y + y;
        [x,y] = get_XY(dt,theta_t,u,x,y,vx,vy,mu);
        [vx,vy] = get_vXY(dt,theta_t,u,vx,vy,mu);
        Te = Te + dt;
        ftheta = ddi_ftheta(u,sqrt(vx^2+vy^2),theta_t,atan2(vy,vx));
        if abs(ftheta) < 1e-3 || abs(ftheta) < abs(ftheta0/100)
            avg_x = sum_x/i;
            avg_y = sum_y/i;
            return
        end
    end
end

