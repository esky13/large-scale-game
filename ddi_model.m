classdef ddi_model
    %DDI_MODEL 带阻尼双积分模型
    % methods: 
    % get_theta
    % get_v
    % ddi_ftheta
    % ddi_fv
    
    properties
        mu  % 阻尼稀疏
        x0  % 初始位置
        y0
        vx0 % 初始速度
        vy0
        u   % 加速度
        x   % 位置
        y
        vx  % 速度
        vy
    end
    
    methods
        function obj = ddi_model(x0,y0,vx0,vy0,u,mu)
            %构造此类的实例，输入初始位置
            obj.x0 = x0; obj.y0 = y0;
            obj.vx0 = vx0; obj.vy0 = vy0;
            obj.mu = mu;
            obj.u = u;
            obj.x = x0; obj.y = y0;
            obj.vx = vx0; obj.vy = vy0;
        end

        function obj = step(obj,dt,theta)
            [obj.x,obj.y] = obj.get_XY(dt,theta,obj.x,obj.y,obj.vx,obj.vy);
            [obj.vx,obj.vy] = obj.get_vXY(dt,theta,obj.vx,obj.vy);
        end

        function [x,y,vx,vy] = get_state(obj)
            x = obj.x;
            y = obj.y;
            vx = obj.vx;
            vy = obj.vy;
        end

        function theta = get_theta(obj)
            % 速度方向
            theta = atan2(obj.vy,obj.vx);
        end

        function v = get_v(obj)
            % 速度大小
            v = sqrt(obj.vx^2+obj.vy^2);
        end

        function ftheta = ddi_ftheta(obj,phi)
            % 当前状态的角度变化率
            % phi: 控制量
            ftheta = obj.u/obj.get_v()*sin(phi-obj.get_theta());
        end

        function ftheta = ddi_ftheta_state(obj,v,theta,phi)
            % 给定状态的角度变化率
            % v: 速度大小, theta: 速度朝向, phi: 控制量
            ftheta = obj.u/v*sin(phi-theta);
        end
        
        function fv = ddi_fv(obj,phi)
            % 速度变化率
            % phi: 控制量
            fv = obj.u*cos(phi-obj.get_theta())-obj.mu*obj.get_v();
        end

        function tc = cal_tc(obj,ts)
            % 计算t1=t2=t3的时刻
            % ts: 最后存在多个时刻的时间
            v02 = (obj.get_v())^2;
            tc = ts/2;
            eps = 1e-5;
            l = 0; r = ts;
            while true
                f = obj.u^2/obj.mu*tc...
                    -exp(-obj.mu*tc)*(...
                                obj.u^2/obj.mu^2 ...
                                -exp(-obj.mu*tc)*(obj.u^2/obj.mu^2-v02));
                if abs(f) < eps
                    break
                end
                if f > 0
                    r = tc;
                    tc = (tc+l)/2;
                else
                    l = tc;
                    tc = (tc+r)/2;
                end
            end
        end

        function ts = cal_ts(obj)
            ts = 1/obj.mu*log((obj.mu*sqrt(obj.vy^2+obj.vx^2)+obj.u)/obj.u);
        end

        function [x,y] = get_XY_from0(obj,t,theta)
            ExpMu = exp(-obj.mu*t);
            x = 1/obj.mu*obj.vx0*(1-ExpMu)+obj.u*cos(theta)/obj.mu.*(t-(1-ExpMu)/obj.mu)+obj.x0;
            y = 1/obj.mu*obj.vy0*(1-ExpMu)+obj.u*sin(theta)/obj.mu.*(t-(1-ExpMu)/obj.mu)+obj.y0;
        end

        function [vx,vy] = get_vXY_from0(obj,t,theta)
            ExpMu = exp(-obj.mu*t);
            vx = obj.vx0*ExpMu+obj.u*cos(theta)/obj.mu*(1-ExpMu);
            vy = obj.vy0*ExpMu+obj.u*sin(theta)/obj.mu*(1-ExpMu);
        end

        function [x,y] = get_XY(obj,dt,theta,x,y,vx,vy)
            % 给定状态得到下一步状态
            ExpMu = exp(-obj.mu*dt);
            x = 1/obj.mu*vx*(1-ExpMu)+obj.u*cos(theta)/obj.mu.*(dt-(1-ExpMu)/obj.mu)+x;
            y = 1/obj.mu*vy*(1-ExpMu)+obj.u*sin(theta)/obj.mu.*(dt-(1-ExpMu)/obj.mu)+y;
        end

        function [vx,vy] = get_vXY(obj,dt,theta,vx,vy)
            % 给定状态得到下一步状态
            ExpMu = exp(-obj.mu*dt);
            vx = vx*ExpMu+obj.u*cos(theta)/obj.mu*(1-ExpMu);
            vy = vy*ExpMu+obj.u*sin(theta)/obj.mu*(1-ExpMu);
        end

        function [rc,xc,yc] = getIsochrones(obj,t,x0,y0,vx0,vy0,u)
            rc = u/obj.mu*(t+exp(-obj.mu*t)/obj.mu-1/obj.mu);
            xc = x0+vx0/obj.mu*(1-exp(-obj.mu*t));
            yc = y0+vy0/obj.mu*(1-exp(-obj.mu*t));
        end

        function [rc,xc,yc] = getmyIsochrones(obj,t)
            % 得到从当前状态出发的等时线
            rc = obj.u/obj.mu*(t+exp(-obj.mu*t)/obj.mu-1/obj.mu);
            xc = obj.x+obj.vx/obj.mu*(1-exp(-obj.mu*t));
            yc = obj.y+obj.vy/obj.mu*(1-exp(-obj.mu*t));
        end

        function index = judge_t(obj,t,theta)
            ts = obj.cal_ts();
            tc = obj.cal_tc(ts);
            [vx_t,vy_t] = obj.get_vXY_from0(t,theta);
            if t > ts
                index = zeros(size(theta));
                return
            end
            index = ones(size(theta));
            if t > tc
                index = index*3;
            end
            index(vy_t.*sin(theta)+vx_t.*cos(theta) < 0) = 2;
        
        end

        function delta_o = cal_delta_o(obj,t,xo,yo,vxo,vyo,uo,out)
            [rAc, xAc, yAc] = obj.getIsochrones(t,obj.x,obj.y,obj.vx,obj.vy,obj.u);
            [rDc, xDc, yDc] = obj.getIsochrones(t,xo,yo,vxo,vyo,uo);
            if out
                delta_o = (xAc-xDc).^2+(yAc-yDc).^2-(rAc+rDc).^2;
            else
                delta_o = (xAc-xDc)^2+(yAc-yDc)^2-(rAc-rDc)^2;
            end
        end
        
        function d_delta_o = cal_d_delta_o(obj,t,delta_x,delta_v,uo,out)
            e = exp(-obj.mu*t);
            if out
                d_delta_o = e*(delta_v'*delta_v/obj.mu*(1-e)+delta_x'*delta_v)... 
                            -(obj.u+uo)^2/obj.mu^2*(t-(1-e)/obj.mu)*(1-e);
            else
                d_delta_o = e*(delta_v'*delta_v/obj.mu*(1-e)+delta_x'*delta_v)... 
                            -(obj.u-uo)^2/obj.mu^2*(t-(1-e)/obj.mu)*(1-e);
            end
        end
        
        function dd_delta_o = cal_dd_delta_o(obj,t,delta_x,delta_v,uo,out)
            e = exp(-obj.mu*t);
            if out
                dd_delta_o = e*(delta_v'*delta_v*(2*e-1)-obj.mu*delta_x'*delta_v)... 
                            -(obj.u+uo)^2/obj.mu^2*((1-e)^2+obj.mu*e*(t-(1-e)/obj.mu));
            else
                dd_delta_o = e*(delta_v'*delta_v*(2*e-1)-obj.mu*delta_x'*delta_v)... 
                            -(obj.u-uo)^2/obj.mu^2*((1-e)^2+obj.mu*e*(t-(1-e)/obj.mu));
            end
        end

        function t = find_t(obj,xo,yo,vxo,vyo,uo,out)
            zero_thr = 1e-5;
            minimal_interval = 0;
            delta_x = [obj.x-xo;obj.y-yo];
            delta_v = [obj.vx-vxo;obj.vy-vyo];
            func = @(t)obj.cal_delta_o(t,xo,yo,vxo,vyo,uo,out);
            d_func = @(t)obj.cal_d_delta_o(t,delta_x,delta_v,uo,out);
            dd_func = @(t)obj.cal_dd_delta_o(t,delta_x,delta_v,uo,out);
            if delta_v'*delta_v+obj.mu*delta_x'*delta_v < 0
                t = find_zero_divide(func,0,obj.mu*norm(delta_x)/(obj.u+uo)+1/obj.mu,minimal_interval,zero_thr);
            elseif delta_x'*delta_v > 0
                t = find_zero_divide(func,0,log(2)/obj.mu,minimal_interval,zero_thr);
            else
                t_d_peak = find_zero_divide(dd_func,0,log(2)/obj.mu,minimal_interval,zero_thr);
                j = d_func(t_d_peak);
                if j < 0
                    t = find_zero_divide(func,0,1,minimal_interval,zero_thr);
                else
                    t_d_zero1 = find_zero_divide(d_func,0,t_d_peak,minimal_interval,zero_thr);
                    t_d_zero2 = find_zero_divide(d_func,t_d_peak,2*t_d_peak,minimal_interval,zero_thr);
                    
                    local_min = func(t_d_zero1);
                    local_max = func(t_d_zero2);
                    if local_min > 0
                        t = find_zero_divide(func,t_d_zero2,2*t_d_zero2,minimal_interval,zero_thr);
                    elseif local_max < 0
                        t = find_zero_divide(func,0,t_d_zero1,minimal_interval,zero_thr);
                    else
                        t1 = find_zero_divide(func,0,t_d_zero1,minimal_interval,zero_thr);
                        t2 = find_zero_divide(func,t_d_zero1,t_d_zero2,minimal_interval,zero_thr);
                        t3 = find_zero_divide(func,t_d_zero2,2*t_d_zero2,minimal_interval,zero_thr);
                        t = [t1,t2,t3];
                    end
                end
            end
        end

        function to = find_to(obj,xo,yo,vxo,vyo,uo)
            to = obj.find_t(xo,yo,vxo,vyo,uo,1);
        end

        function ti = find_ti(obj,xo,yo,vxo,vyo,uo)
            ti = obj.find_t(xo,yo,vxo,vyo,uo,0);
        end

        function Te = simulate_to_Te(obj,tx,ty)
            % 从当前状态仿真到收敛时刻
            Te = 0;
            dt = 0.1;
            theta_t = obj.ddi_target_theta(tx,ty);
            ftheta0 = obj.ddi_ftheta(theta_t);
            i = 0;
            vx_t = obj.vx;
            vy_t = obj.vy;
            while true
                i = i+1;
                Te = Te + dt;
                [vx_t,vy_t] = obj.get_vXY(dt,theta_t,vx_t,vy_t);
                ftheta = obj.ddi_ftheta_state(sqrt(vx_t^2+vy_t^2),atan2(vy_t,vx_t),theta_t);
                if abs(ftheta) < 1e-3 || abs(ftheta) < abs(ftheta0/100)
                    return
                end
            end
        end

        function [theta,t] = ddi_target_theta(obj,tx,ty)
            t = obj.find_to(tx,ty,0,0,0);
            t = t(1);
            [~, xc, yc] = obj.getmyIsochrones(t);
            theta = atan2(ty-yc,tx-xc);
        end

    end
end

