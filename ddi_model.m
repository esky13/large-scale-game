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
        function obj = ddi_model(x0,y0,vx0,vy0,mu,u)
            %构造此类的实例，输入初始位置
            obj.x0 = x0; obj.y0 = y0;
            obj.vx0 = vx0; obj.vy0 = vy0;
            obj.mu = mu;
            obj.u = u;
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
            % 角度变化率
            % phi: 控制量
            ftheta = obj.u/obj.get_v()*sin(phi-obj.get_theta());
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

        function obj = state_update(obj,dt,theta)
            ExpMu = exp(-obj.mu*dt);
            obj.x = 1/obj.mu*obj.vx*(1-ExpMu)+obj.u*cos(theta)/obj.mu.*(t-(1-ExpMu)/obj.mu)+obj.x;
            obj.y = 1/obj.mu*obj.vy*(1-ExpMu)+obj.u*sin(theta)/obj.mu.*(t-(1-ExpMu)/obj.mu)+obj.y;
            obj.vx = obj.vx*ExpMu+obj.u*cos(theta)/obj.mu*(1-ExpMu);
            obj.vy = obj.vy*ExpMu+obj.u*sin(theta)/obj.mu*(1-ExpMu);
        end

        function delta_o = cal_delta_o(t,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu,out)
            [rAc, xAc, yAc] = getIsochrones(t,xA0,yA0,vAx0,vAy0,uA,mu);
            [rDc, xDc, yDc] = getIsochrones(t,xD0,yD0,vDx0,vDy0,uD,mu);
            if out
                delta_o = (xAc-xDc).^2+(yAc-yDc).^2-(rAc+rDc).^2;
            else
            end
        end
        
        function d_delta_o = cal_d_delta_o(t,delta_x,delta_v,uA,uD,mu,out)
            e = exp(-mu*t);
            d_delta_o = e*(delta_v'*delta_v/mu*(1-e)+delta_x'*delta_v)... 
                -(uA+uD)^2/mu^2*(t-(1-e)/mu)*(1-e);
        end
        
        function dd_delta_o = cal_dd_delta_o(t,delta_x,delta_v,uA,uD,mu,out)
            e = exp(-mu*t);
            dd_delta_o = e*(delta_v'*delta_v*(2*e-1)-mu*delta_x'*delta_v)... 
                -(uA+uD)^2/mu^2*((1-e)^2+mu*e*(t-(1-e)/mu));
        end

        function to = find_to(xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu)
            delta_x = [xA0-xD0;yA0-yD0];
            delta_v = [vAx0-vDx0;vAy0-vDy0];
            func = @(t)cal_delta_o(t,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu);
            d_func = @(t)cal_d_delta_o(t,delta_x,delta_v,uA,uD,mu);
            dd_func = @(t)cal_dd_delta_o(t,delta_x,delta_v,uA,uD,mu);
            if delta_v'*delta_v+mu*delta_x'*delta_v < 0
                to = find_zero_divide(func,0,mu*norm(delta_x)/(uA+uD)+1/mu);
            elseif delta_x'*delta_v > 0
                to = find_zero_divide(func,0,log(2)/mu);
            else
                t_d_peak = find_zero_divide(dd_func,0,log(2)/mu);
                j = d_func(t_d_peak);
                if j < 0
                    to = find_zero_divide(func,0,1);
                else
                    t_d_zero1 = find_zero_divide(d_func,0,t_d_peak);
                    t_d_zero2 = find_zero_divide(d_func,t_d_peak,2*t_d_peak);
                    
                    local_min = func(t_d_zero1);
                    local_max = func(t_d_zero2);
                    if local_min > 0
                        to = find_zero_divide(func,t_d_zero2,2*t_d_zero2);
                    elseif local_max < 0
                        to = find_zero_divide(func,0,t_d_zero1);
                    else
                        t1 = find_zero_divide(func,0,t_d_zero1);
                        t2 = find_zero_divide(func,t_d_zero1,t_d_zero2);
                        t3 = find_zero_divide(func,t_d_zero2,2*t_d_zero2);
                        to = [t1,t2,t3];
                    end
                end
            end
        end
    end
end

