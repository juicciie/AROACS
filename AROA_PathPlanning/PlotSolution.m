function PlotSolution(sol,model)
    Colors = hsv(length(model.robs));

    xs=model.xs;
    ys=model.ys;
    xt=model.xt;
    yt=model.yt;
    xobs=model.xobs;
    yobs=model.yobs;
    robs=model.robs;
    
    XS=sol.XS;
    YS=sol.YS;
    xx=sol.xx;
    yy=sol.yy;
    
    theta=linspace(0,2*pi,100);
    for k=1:numel(xobs)
        Color=Colors(k,:);
        White=[1 1 1];
        
        Color=0.4*Color+0.6*White;
        fill(xobs(k)+robs(k)*cos(theta),yobs(k)+robs(k)*sin(theta),Color);
        hold on;
    end
    plot(xx,yy,'k','LineWidth',2);
    plot(XS,YS,'ro');
    plot(xs,ys,'bs','MarkerSize',12,'MarkerFaceColor','y');
    plot(xt,yt,'kp','MarkerSize',16,'MarkerFaceColor','g');
    hold off;
    grid on;
    axis equal;

end