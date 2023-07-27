function plotOut(tspan,yFOM,yROM_POD,yROM_NLMM)

figure
subplot(1,2,1)
hold on
box on
grid on
plot(tspan, yFOM(:,1),'b','Linewidth',4)
%plot(tspan, yROM_POD(:,1),'g:','Linewidth',4)
plot(tspan, yROM_NLMM(:,1),'r--','Linewidth',4)
xlim([0 tspan(end)])
xlabel('t(sec)','Interpreter','LaTeX')
ylabel('Avg. $\delta$ (rad)','Interpreter','LaTeX')
legend('FOM','NLMM-DMD','AutoUpdate','On','Interpreter','LaTeX');
set(gca,'FontSize',20,'TickLabelInterpreter','Latex')
if size(yFOM,2) >= 2
subplot(1,2,2)
hold on
grid on
box on
plot(tspan, yFOM(:,2),'b','Linewidth',4)
%plot(tspan, yROM_POD(:,2),'g:','Linewidth',4)
plot(tspan, yROM_NLMM(:,2),'r--','Linewidth',4)
xlabel('t(sec)','Interpreter','LaTeX')
ylabel('Avg. $\omega$ (rad/sec)','Interpreter','LaTeX')
set(gca,'FontSize',20,'TickLabelInterpreter','Latex')
legend('FOM','NLMM-DMD','AutoUpdate','On','Interpreter','LaTeX');

end


end