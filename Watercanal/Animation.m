function Animation(V,omega,t)

pause(1)
dt = t(2)-t(1);
PlotParameters
% close all
figure
set(gcf,'units','points','position',[x0screen,y0screen,WidthScreen,HeightScreen])
%Creating A380 object and its trace line 
a380=A380;
traceline=line(0,0,0,'LineWidth',2);
titre = title(sprintf('$t=%.1f$',t(1)),'Interpreter','latex','fontsize',fs);
%Animating yaw motion
for i=1:length(t)
    a380.MoveForward(V(i)*dt);
    a380.AddPitch(omega(i)*dt);
    set(traceline,'XData',[get(traceline,'XData'),a380.cgPos(1)]);
    set(traceline,'YData',[get(traceline,'YData'),a380.cgPos(2)]);
    set(traceline,'ZData',[get(traceline,'ZData'),a380.cgPos(3)]);
    a380.ShowStarboardSide;
    axis([a380.cgPos(1)-2*a380.overallLength,...
          a380.cgPos(1)+2*a380.overallLength,...
          -a380.overallLength,a380.overallLength,...
          a380.cgPos(3)-2*a380.overallLength,...
          a380.cgPos(3)+2*a380.overallLength])
      set(titre,'String',sprintf('$t=%.1f$',t(i)))
%     pause(t(i)-t(i-1));
pause(dt*0.01);
end

end