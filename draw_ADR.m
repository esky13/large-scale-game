function [bounderx,boundery]=draw_ADR(ti,to,xA0,yA0,vAx0,vAy0,uA,xD0,yD0,vDx0,vDy0,uD,mu,ax1,c,filled)
    t_step = [to:(ti-to)/500:to+(ti-to)/10,to+(ti-to)/10:(ti-to)/50:ti-(ti-to)/20,ti-(ti-to)/20:(ti-to)/500:ti];
    bounder1x = [];bounder1y = [];
    bounder2x = [];bounder2y = [];
    indexD1 = []; indexD2 = []; indexA1 = []; indexA2 = [];

    tDs = 1/mu*log((mu*sqrt(vDy0^2+vDx0^2)+uD)/uD);
    tDc = cal_tc(tDs,vDx0,vDy0,uD,mu);
    tAs = 1/mu*log((mu*sqrt(vAy0^2+vAx0^2)+uA)/uA);
    tAc = cal_tc(tAs,vAx0,vAy0,uA,mu);
    for t = t_step
        [rDc, xDc, yDc] = getIsochrones(t,xD0,yD0,vDx0,vDy0,uD,mu);
        [rAc, xAc, yAc] = getIsochrones(t,xA0,yA0,vAx0,vAy0,uA,mu);
        [hasroot,x1,y1,x2,y2]=getADIntersection(rDc,xDc,yDc,rAc,xAc,yAc);
        if hasroot
            thetaD1 = atan2(y1-yDc,x1-xDc); thetaD2 = atan2(y2-yDc,x2-xDc);
            thetaA1 = atan2(y1-yAc,x1-xAc); thetaA2 = atan2(y2-yAc,x2-xAc);
            indexD1(end+1) = judge_t(t,thetaD1,tDc,tDs,uD,vDx0,vDy0,mu);
            indexD2(end+1) = judge_t(t,thetaD2,tDc,tDs,uD,vDx0,vDy0,mu);
            indexA1(end+1) = judge_t(t,thetaA1,tAc,tAs,uA,vAx0,vAy0,mu);
            indexA2(end+1) = judge_t(t,thetaA2,tAc,tAs,uA,vAx0,vAy0,mu);
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