%%
t = Simulation.simout.Time';
X = Simulation.simout.Data;
X = X(:,:);

h_d = X(1:N,:);
u_d = X(N+1:2*N,:);

Y = Ctrue*X;

lw = 4;
fs = 18;
x0screen=100;y0screen=100;WidthScreen=1200;HeightScreen=600;

figure
plot(t,Y,'LineWidth',lw)
set(gca,'FontSize',fs)
grid on
xlabel('time $t\, (s)$','Interpreter','latex')
legend({'$h(0,t)$','$h(L,t)$'},'Interpreter','latex')




figure
set(gcf,'units','points','position',[x0screen,y0screen,WidthScreen,HeightScreen])
i  = 1;
hold on
ZETA = [zeta(1:end-1);zeta(2:end)];
ZETA = [ZETA;flipud(ZETA)];
y1 = h_d(:,i)';
y2 = h_d(:,i)'*0;
H_d =[y1(1:end-1);y1(2:end); y2(1:end-1);y2(2:end)]; 

CC   = u_d(1:end-1,i)';

p_h = patch(ZETA,H_d,CC,'EdgeColor','none','CDataMapping','scaled');
% p_b = plot([zeta(end),zeta(end)],[U(i),100],'LineWidth',lw*4,'Color',[0,0,0])
set(gca,'FontSize',fs)
grid on
xlabel('space $\zeta$','Interpreter','latex')
xlim([a,b])
% ylim([35,45])
ylim([0,5])
colorbar
caxis([-1,1]);

uiwait(msgbox("Start Simulation"));
for i = 1:length(t)
    y1 = h_d(:,i)';
    H_d =[y1(1:end-1);y1(2:end); y2(1:end-1);y2(2:end)]; 
    CC   = u_d(1:end-1,i)';    
    set(p_h,'YData',H_d,'CData',CC)
    % set(p_b,'YData',[U(i),100])
    pause(dt*4)
    % if i == 1;
    %     f = msgbox("Start Simulation")
    %     % pause
    % end

end