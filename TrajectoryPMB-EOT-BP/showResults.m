function showResults(X,Y)

figure
box on
hold on

for i = 1:length(X.tVec)
    plot(X.xState(1,X.tVec(i),i),X.xState(2,X.tVec(i),i),'kx','MarkerSize',10,'LineWidth',2);
    duration = X.tVec(i):X.tVec(i)+X.iVec(i)-1;
    h1 = plot(X.xState(1,duration,i),X.xState(2,duration,i),'k','LineWidth',2);
end

for i = 1:length(Y.tVec)
    plot(Y.xState(1,Y.tVec(i),i),Y.xState(2,Y.tVec(i),i),'rx','MarkerSize',10,'LineWidth',2);
    duration = Y.tVec(i):Y.tVec(i)+Y.iVec(i)-1;
    h2 = plot(Y.xState(1,duration,i),Y.xState(2,duration,i),'r','LineWidth',2);
end

legend([h1 h2],'Ground truth','Estimates','Interpreter','latex')

xlim([-150 150])
ylim([-150 150])

xlabel('$e_1~[m]$','Interpreter','latex')
ylabel('$e_2~[m]$','Interpreter','latex')

set(gca,'TickLabelInterpreter','latex','FontSize',16)

end