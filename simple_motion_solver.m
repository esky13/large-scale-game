classdef simple_motion_solver
    
    properties
        
    end
    
    methods
        function obj = simple_motion_solver()
        end
        
        function [thetaP,thetaE,tx,ty] = find_optimal_RA(obj,xP,yP,xE,yE,alpha)
            xc = 1/(1-alpha^2)*(xE-alpha^2*xP);
            yc = 1/(1-alpha^2)*(yE-alpha^2*yP);
            rc = alpha/(1-alpha^2)*(sqrt((xE-xP)^2+(yE-yP)^2));
            tx = xc-cos(atan2(yc,xc))*rc;
            ty = yc-sin(atan2(yc,xc))*rc;
            thetaP = atan2(ty-yP,tx-xP);
            thetaE = atan2(ty-yE,tx-xE);
        end

        function [thetaP,thetaE,tx,ty] = find_optimal_PE(obj,xP,yP,xE,yE,alpha)
            xc = 1/(1-alpha^2)*(xE-alpha^2*xP);
            yc = 1/(1-alpha^2)*(yE-alpha^2*yP);
            rc = alpha/(1-alpha^2)*(sqrt((xE-xP)^2+(yE-yP)^2));
            tx = xc+cos(atan2(yE-yP,xE-xP))*rc;
            ty = yc+sin(atan2(yE-yP,xE-xP))*rc;
            thetaP = atan2(ty-yP,tx-xP);
            thetaE = atan2(ty-yE,tx-xE);
        end
    end
end

