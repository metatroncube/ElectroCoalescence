x=0.2;
y=0;

%(2^1.5-(2^0.75+1))/(2^1.5-2^0.75)

plot(4,-4,'.k','MarkerSize',20)
%text(4+x,-4+y,'100%','Color','black')
hold on

plot(5.8,1.7,'.r','MarkerSize',20)
%text(6,1.6,'100%','Color','red')
hold on

%text(4.7+x,-1+y,'12.8%','Color','blue')
%text(5.7+x,-1+y,'87.2%','Color','blue')
plot(5,2,'.b','MarkerSize',20)
text(3.8,2+y,'10.6%','Color','blue')
hold on
plot(5,4,'.k','MarkerSize',20)
%text(5+x,4+y,'100%','Color','black')
hold on
plot(5,3,'.b','MarkerSize',20)
text(5+x,3+y,'2.2%','Color','blue')
hold on
plot(6,1,'.b','MarkerSize',20)
text(6+x,1+y,'30.0%','Color','blue')
hold on
plot(6,2,'.b','MarkerSize',20)
text(6+x,2+y,'57.2%','Color','blue')
%1, 1.68, 2.83 1.414
grid on
axis([0,12,-5,5])
set(gca,'XTick',0:1:15);
set(gca,'FontSize',10);
xlabel('radius bins, 4log_2(r/2)')
ylabel('charge bins, ¡Àlog_2|4q/r^2|')
legend('droplets before collision-coalescence','droplets formed after collision-coalescence','droplets redistributed into new bins')