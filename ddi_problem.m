classdef ddi_problem
    %DDI_PROBLEM 理论求解带阻尼双积分微分博弈问题
    %   此处显示详细说明
    
    properties
        ddi_model_A
        ddi_model_D
    end
    
    methods
        function obj = ddi_problem(xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu)
            %DDI_PROBLEM 构造此类的实例
            %   此处显示详细说明
            obj.ddi_model_A = ddi_model(xA0,yA0,vAx0,vAy0,uA,mu);
            obj.ddi_model_D = ddi_model(xD0,yD0,vDx0,vDy0,uD,mu);
        end

        function obj = step(obj,dt,thetaA,thetaD)
            obj.ddi_model_A = obj.ddi_model_A.step(dt,thetaA);
            obj.ddi_model_D = obj.ddi_model_D.step(dt,thetaD);
        end

        function [xA,yA,vAx,vAy,xD,yD,vDx,vDy] = get_state(obj)
            [xA,yA,vAx,vAy] = obj.ddi_model_A.get_state();
            [xD,yD,vDx,vDy] = obj.ddi_model_D.get_state();
        end
        
        function to = find_to(obj)
            to = obj.ddi_model_A.find_to(obj.ddi_model_D.x,obj.ddi_model_D.y,obj.ddi_model_D.vx,obj.ddi_model_D.vy,obj.ddi_model_D.u);
        end

        function ti = find_ti(obj)
            ti = obj.ddi_model_A.find_ti(obj.ddi_model_D.x,obj.ddi_model_D.y,obj.ddi_model_D.vx,obj.ddi_model_D.vy,obj.ddi_model_D.u);
        end

        function [bounderx,boundery]=draw_ADR(obj,ti,to,ax1,c,filled)
            t_step = [to:(ti-to)/500:to+(ti-to)/10,to+(ti-to)/10:(ti-to)/50:ti-(ti-to)/20,ti-(ti-to)/20:(ti-to)/500:ti];
            bounder1x = [];bounder1y = [];
            bounder2x = [];bounder2y = [];
            indexD1 = []; indexD2 = []; indexA1 = []; indexA2 = [];
        
            for t = t_step
                [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t);
                [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t);
                [hasroot,x1,y1,x2,y2] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
                if hasroot
                    thetaD1 = atan2(y1-yDc,x1-xDc); thetaD2 = atan2(y2-yDc,x2-xDc);
                    thetaA1 = atan2(y1-yAc,x1-xAc); thetaA2 = atan2(y2-yAc,x2-xAc);
                    indexD1(end+1) = obj.ddi_model_D.judge_t(t,thetaD1);
                    indexD2(end+1) = obj.ddi_model_D.judge_t(t,thetaD2);
                    indexA1(end+1) = obj.ddi_model_A.judge_t(t,thetaA1);
                    indexA2(end+1) = obj.ddi_model_A.judge_t(t,thetaA2);
                    bounder1x(end+1) = x1;bounder1y(end+1) = y1;
                    bounder2x(end+1) = x2;bounder2y(end+1) = y2;
                end
            end
            bounderx = [flip(bounder1x),bounder2x];
            boundery = [flip(bounder1y),bounder2y];
            indexD = [flip(indexD1),indexD2];
            indexA = [flip(indexA1),indexA2];
        
            figure(ax1)
            
            if filled
                fill(bounderx,boundery,[0.9290 0.6940 0.1250],'FaceAlpha',0.5,'EdgeColor','none');
            else
                plot(bounderx,boundery,'Color',c,'LineWidth',2,'HandleVisibility','off');
                % plot(bounderx(and(indexD==2,indexA==1)),boundery(and(indexD==2,indexA==1)),'Color',[1,0,0],'LineWidth',2,'DisplayName','12');
                % plot(bounderx(and(indexD==1,indexA==1)),boundery(and(indexD==1,indexA==1)),'Color',[1,1,0],'LineWidth',2,'DisplayName','11');
                % plot(bounderx(and(indexD==2,indexA==2)),boundery(and(indexD==2,indexA==2)),'Color',[0,1,0],'LineWidth',2,'DisplayName','22');
                % plot(bounderx(and(indexD==3,indexA==2)),boundery(and(indexD==3,indexA==2)),'Color',[0.2,0.5,1],'LineWidth',2,'DisplayName','32');
                % plot(bounderx(and(indexD==1,indexA==2)),boundery(and(indexD==1,indexA==2)),'Color',[1,0,1],'LineWidth',2,'DisplayName','21');
                % plot(bounderx(and(indexD==3,indexA==1)),boundery(and(indexD==3,indexA==1)),'Color',[0,1,1],'LineWidth',2,'DisplayName','13');
                % plot(bounderx(and(indexD==2,indexA==0)),boundery(and(indexD==2,indexA==0)),'Color',[1,0,0],'LineWidth',2,'DisplayName','02');
                % plot(bounderx(and(indexD==1,indexA==0)),boundery(and(indexD==1,indexA==0)),'Color',[1,1,0],'LineWidth',2,'DisplayName','01');
                % legend
            end
        end

        function draw_trajectory(thetaA,thetaD,tf,ax)
            t = 0:tf/100:tf;
            sA = zeros([length(t),2]);
            sD = zeros([length(t),2]);
            [sA(:,1),sA(:,2)] = obj.ddi_model_A.get_XY(t,thetaA,obj.ddi_model_A.x,obj.ddi_model_A.y, ...
                                                        obj.ddi_model_A.vx,obj.ddi_model_A.vy);
            [sD(:,1),sD(:,2)] = obj.ddi_model_D.get_XY(t,thetaD,obj.ddi_model_D.x,obj.ddi_model_D.y, ...
                                                        obj.ddi_model_D.vx,obj.ddi_model_D.vy);
            figure(ax);
            plot(sA(:,1),sA(:,2),'r','LineWidth',1);
            plot(sD(:,1),sD(:,2),'b','LineWidth',1);
        end

        function [hasroot,x1,y1,x2,y2] = getADIntersection(obj,rAc,xAc,yAc,rDc,xDc,yDc)
            a = 4*((xAc-xDc)^2+(yAc-yDc)^2);
            if abs(xAc-xDc) < abs(yAc-yDc)
                b = 4*((rAc-rDc)*(rAc+rDc)*(xAc-xDc)-(xAc+xDc)*((xAc-xDc)^2+(yAc-yDc)^2));
                c = ((rAc-rDc)*(rAc+rDc)-(xAc-xDc)*(xAc+xDc))^2+(yAc-yDc)^2*((yAc-yDc)^2+2*xAc^2+2*xDc^2-2*rAc^2-2*rDc^2);
                delta = b^2-4*a*c;
                x1 = (-b+sqrt(delta))/2/a;
                x2 = (-b-sqrt(delta))/2/a;
                y1 = -(xAc-xDc)/(yAc-yDc)*x1+(xAc^2-xDc^2+yAc^2-yDc^2-rAc^2+rDc^2)/2/(yAc-yDc);
                y2 = -(xAc-xDc)/(yAc-yDc)*x2+(xAc^2-xDc^2+yAc^2-yDc^2-rAc^2+rDc^2)/2/(yAc-yDc);
            else
                b = 4*(yAc-yDc)*(rAc-rDc)*(rAc+rDc)-4*(yAc+yDc)*((xAc-xDc)^2+(yAc-yDc)^2);
                c = (rAc^2-rDc^2-yAc^2+yDc^2)^2+(xAc-xDc)^2*((xAc-xDc)^2+2*(yAc^2+yDc^2)-2*(rAc^2+rDc^2));
                delta = b^2-4*a*c;
                y1 = (-b-sqrt(delta))/2/a;
                y2 = (-b+sqrt(delta))/2/a;
                x1 = -(yAc-yDc)/(xAc-xDc)*y1+(xAc^2-xDc^2+yAc^2-yDc^2-rAc^2+rDc^2)/2/(xAc-xDc);
                x2 = -(yAc-yDc)/(xAc-xDc)*y2+(xAc^2-xDc^2+yAc^2-yDc^2-rAc^2+rDc^2)/2/(xAc-xDc);
            end
            hasroot = delta >= 0;
            if (x1-xAc)*(yDc-yAc)-(y1-yAc)*(xDc-xAc) < 0
                temp = x1; x1 = x2; x2 = temp;
                temp = y1; y1 = y2; y2 = temp;
            end
        end

        function [tx,ty,op_thetaA,op_thetaD,t_minimal,ti,to] = find_optimal_RA(obj)
            a = obj.find_to();
            b = obj.find_ti();
            ti = b(1);
            to = a(1);
            if length(a) == 3
                if a(3) < ti
                    to = a(3);
                end
            end
        
            [t_minimal,op_thetaA,op_thetaD,valid] = obj.attacker_win_strategy(to,ti);
            tx = 0;
            ty = 0;
            if ~valid
                [t_minimal,op_thetaA,op_thetaD,tx,ty] = obj.get_dis_minimal(to,ti);
            end
        end

        function [tx,ty,op_thetaA,op_thetaD,ti] = find_optimal_PE(obj)
            b = obj.find_ti();
            ti = b(1);
        
            [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(ti);
            [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(ti);
            [~,x1,y1,x2,y2] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
            tx = real(x1+x2)/2;
            ty = real(y1+y2)/2;
            op_thetaD = atan2(ty-yDc,tx-xDc);
            op_thetaA = atan2(ty-yAc,tx-xAc);
        end

        function [t_reach,theta_reach,thetaD,valid] = attacker_win_strategy(obj,tmin,tmax)
            tf = obj.ddi_model_A.find_to(0,0,0,0,0);
            
            t_reach = 0;
            theta_reach = 0;
            thetaD = 0;
            valid = false;
            max_dis = -Inf;
            for t = tf
                t_reach = t;%
                if t > tmax
                    continue
                elseif t < tmin
                    valid = true;
                else
                    [~, xc, yc] = obj.ddi_model_A.getmyIsochrones(t);
                    thetaA = atan2(-yc,-xc);
                    step = (t-tmin)/100; t_cur = t;
                    while true
                        [xAt,yAt] = obj.ddi_model_A.get_XY(t_cur,thetaA,obj.ddi_model_A.x,obj.ddi_model_A.y,obj.ddi_model_A.vx,obj.ddi_model_A.vy);
                        [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t_cur);
                        dis_cur = (xDc-xAt)^2+(yDc-yAt)^2-rDc^2;
                        if dis_cur < 0
                            break
                        end
                        if t_cur < tmin
                            valid = true;
                            break
                        end
                        t_cur = t_cur - step;
                    end
                end
                if valid
                    [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t);
                    dis = sqrt(xDc^2+yDc^2)-rDc;
                    if dis > max_dis
                        max_dis = dis;
                        thetaD = atan2(-yDc,-xDc);
                        t_reach = t;
                        [~, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t);
                        theta_reach = atan2(-yAc,-xAc);
                    end
                end
            end
        end

        function [t_minimal,thetaA,thetaD,sx1,sy1] = get_dis_minimal(obj,tmin,tmax)
            delta = (tmax-tmin)/100;
            minimal_interval = 1e-10;
            prevs1 = 0; prevs2 = 0;
            dis_min = Inf;

            t_step = tmin:delta:tmax;
            for t = t_step
                [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t);
                [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t);
                [hasroot,x1,y1,x2,y2] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
                if ~hasroot
                    continue
                end
                if x1^2+y1^2<dis_min
                    dis_min = x1^2+y1^2;
                    thetaA = atan2(y1-yAc,x1-xAc);
                    thetaD = atan2(y1-yDc,x1-xDc);
                    t_minimal = t;
                    sx1=x1;sy1=y1;
                end
                if x2^2+y2^2<dis_min
                    dis_min = x2^2+y2^2;
                    thetaA = atan2(y2-yAc,x2-xAc);
                    thetaD = atan2(y2-yDc,x2-xDc);
                    t_minimal = t;
                    sx1=x2;sy1=y2;
                end
                
                cal_dis_derivative1 = @(tt)obj.cal_dis_derivative(tt,1);
                s1 = cal_dis_derivative1(t);
                if prevs1*s1 < 0
                    left = t-delta; right = t;
                    t_cur = find_zero_divide(cal_dis_derivative1, left, right, minimal_interval, 1e-5);
                    [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t_cur);
                    [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t_cur);
                    [~,x1,y1,~,~] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
                    if x1^2+y1^2 < dis_min
                        t_minimal = t_cur;
                        dis_min = x1^2+y1^2;
                        thetaA = atan2(y1-yAc,x1-xAc);
                        thetaD = atan2(y1-yDc,x1-xDc);
                        sx1=x1;sy1=y1;
                    end
                end
                    
                cal_dis_derivative2 = @(tt)obj.cal_dis_derivative(tt,0);
                s2 = cal_dis_derivative2(t);
                if prevs2*s2 < 0
                    left = t-delta; right = t;
                    t_cur = find_zero_divide(cal_dis_derivative2, left, right, minimal_interval, 1e-5);
                    [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t_cur);
                    [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t_cur);
                    [~,~,~,x2,y2] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
                    if x2^2+y2^2 < dis_min
                        thetaA = atan2(y2-yAc,x2-xAc);
                        thetaD = atan2(y2-yDc,x2-xDc);
                        t_minimal = t_cur;
                        dis_min = x2^2+y2^2;
                        sx1=x2; sy1=y2;
                    end
                end
            
                prevs1 = s1; prevs2 = s2;
            end
        end

        function s = cal_dis_derivative(obj,t,root1)
            [rDc, xDc, yDc] = obj.ddi_model_D.getmyIsochrones(t);
            [rAc, xAc, yAc] = obj.ddi_model_A.getmyIsochrones(t);
            if root1
                [~,xf,yf,~,~] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
            else
                [~,~,~,xf,yf] = obj.getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
            end
            thetaA = atan2(yf-yAc,xf-xAc);
            thetaD = atan2(yf-yDc,xf-xDc);
            [vAx,vAy] = obj.ddi_model_A.get_vXY(t,thetaA,obj.ddi_model_A.vx,obj.ddi_model_A.vy);
            [vDx,vDy] = obj.ddi_model_D.get_vXY(t,thetaD,obj.ddi_model_D.vx,obj.ddi_model_D.vy);
            phi = atan2(yf,xf);
            s = (sin(thetaA-phi)*(vDx*cos(thetaD)+vDy*sin(thetaD)) ...
                -sin(thetaD-phi)*(vAx*cos(thetaA)+vAy*sin(thetaA)))/sin(thetaA-thetaD);
        end

    end
end

