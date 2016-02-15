

close all;

%load Bus30.mat
load 30BusExp1
figure;
 hold on;
%title('Objective values')
plot(F_t(:,2),'r<-','LineWidth',2.5,'MarkerSize',6);
%plot(F_t(:,2),'>-','LineWidth',2.5,'MarkerSize',6);
plot(F_t(:,3),'g*-','LineWidth',2.5,'MarkerSize',6);
plot(F_t(:,4),'o-','LineWidth',2.5,'MarkerSize',6);
hold off
box on; grid on;
xlim([0,100])


set(gca,'FontSize',32,'FontName', 'Arial')
xlabel('Number of iteration')
ylabel('F-norm of $$X-\hat{X}$$','Interpreter','Latex')
h = legend('\rho = 10','\rho = 20','\rho = 50');
set(h,'FontSize',32)


figure;
hold on;
%title('Norm of resisual')
plot(R_t(:,2),'r<-','LineWidth',2.5,'MarkerSize',6);
%plot(R_t(:,2),'>-','LineWidth',2.5,'MarkerSize',6);
plot(R_t(:,3),'g*-','LineWidth',2.5,'MarkerSize',6);
plot(R_t(:,4),'o-','LineWidth',2.5,'MarkerSize',6);
hold off
box on; grid on;
xlim([0,100])

set(gca,'FontSize',32,'FontName', 'Arial')
xlabel('Number of iteration')
ylabel('F-norm of $$X^k-Y^k$$','Interpreter','Latex')
h = legend('\rho = 10','\rho = 20','\rho = 50');
set(h,'FontSize',32)


