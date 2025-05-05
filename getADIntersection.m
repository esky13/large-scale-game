function [hasroot,x1,y1,x2,y2] = getADIntersection(rAc,xAc,yAc,rDc,xDc,yDc)
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