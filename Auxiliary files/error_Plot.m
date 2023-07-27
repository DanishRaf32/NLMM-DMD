function error_Plot(tspan,Error)

figure()
subplot(1,2,1)
hold on
set(gca, 'YScale', 'log')
%plot(tspan, Error.errorPODout1,'g','Linewidth',3)
plot(tspan, Error.errorNLMMOut1,'b','Linewidth',3)
title('Avg. delta','Interpreter','LaTeX')
xlim([0 tspan(end)])
xlabel('$t$','Interpreter','LaTex')
ylabel('$e(t)$','Interpreter','LaTex')
%legend('POD','NLMM-DMD','Interpreter','LaTex');
grid on
box on
set(gca,'FontSize',20,'TickLabelInterpreter','Latex')
subplot(1,2,2)
hold on
set(gca, 'YScale', 'log')
%plot(tspan, Error.errorPODout2,'g','Linewidth',3)
plot(tspan, Error.errorNLMMOut2,'b','Linewidth',3)
title('Avg. omega','Interpreter','LaTeX')
xlim([0 tspan(end)])
xlabel('$t$','Interpreter','LaTex')
ylabel('$e(t)$','Interpreter','LaTex')
%legend('POD','NLMM-DMD','Interpreter','LaTex');
grid on
box on
set(gca,'FontSize',20,'TickLabelInterpreter','Latex')
end