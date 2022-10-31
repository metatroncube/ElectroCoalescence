R1=[30,40,10,20];
figure(1)
for k=1:length(R1)
if k==3
    figure(2)
end
draw=load ( strcat('D:\Fortran\collision_efficiency_charge\r1=',num2str(R1(k)),'\draw.csv'  ));
i=1;
j=2;
left=2;
color=[[0,0,0];[255,0,0];[255,155,0] ;[0,200,0];[255,0,255];[0,0,255]]/255;
r1=draw(1,1);
subplot(2,1,mod(k-1,2)+1)

n1=1;n2=1;
i=1;
for j=2:1000
    %if draw(i,1)<draw()
    if j>=length(draw) || draw(j+1,1)<draw(j,1)      
        n1=n2;
        n2=j;
        xx1 = 1:0.1:0.9*r1;
        yy1 = interp1(draw(n1+1:n2,1),draw(n1+1:n2,2),xx1,'spline');
        semilogy(xx1 ,yy1,'Color',color(i,:),'LineWidth',1);%
        i=i+1;
        hold on
    end
    j=j+1;
    if j>length(draw)
        break
    end 
end
axis([0,r1,0.01,1.5])
xlabel('r_2 (¦Ìm)')
ylabel('Collision Efficiency E')
if k==1 || k==3
    title(strcat('(a) r_1=',num2str(r1),'¦Ìm'))
else
    title(strcat('(b) r_1=',num2str(r1),'¦Ìm'))
end
 if k==2
    h=legend('(1) no charge','(2) q= +32r^2, no field ','(3) q= +32r^2, Field=400V cm^{-1}','(4) q= -32r^2, Field=400V cm^{-1}','(5) q= ¡À32r^2, no field ','(6) q= ¡À32r^2, Field=400V cm^{-1}')
    set(h,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
 end
set(gca,'FontSize',10);
set (gcf,'Position',[600,100,600,900])
end
hold off