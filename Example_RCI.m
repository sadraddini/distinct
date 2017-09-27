clear
clc
X=importdata('state.txt');
R=importdata('RCI.txt');
N=(length(R))-1;
hold on
[x,y]=meshgrid(1:-2/N:-1);
surf(x,y,R','FaceAlpha',0.25,'LineWidth',0.01,'EdgeColor','none')
plot(X(:,1),X(:,2),'Linewidth',1)
plot(X(:,1),X(:,2),'o','Linewidth',3)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('$\rho^*=0.02,\eta=0.1,\epsilon=0.01,K=6$','Interpreter','latex')
hold off