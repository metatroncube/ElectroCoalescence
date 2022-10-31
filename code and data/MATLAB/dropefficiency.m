draw=load ('D:\Fortran\collision_efficiency\draw.csv');
c40=[[4.7602	0.73616
];[7.1961	1.38254
];[7.9922	1.50735
];[8.7552	1.62198
]];

c55=[[4.425	1.09945
];[5.5199	1.24978
];[9.8168	1.78247
];[10.9955	1.8261
]];
cbeard=[[4.0586	1.22384
];[4.3215	1.56192
];[4.7359	1.66627
];[5.5487	1.78092
];[8.0056	1.89109
]];
cbeardradius=[73,113,155,250,250];
ckinzer=[[8.0402	1.73608
]];
cjonas=[9.0204	1.68561];
i=1;
j=2;
left=2;
color=[ [9,0,0];[90,0,0]; [145,0,0];[200,0,0];[255,0,0];[255,122,0]; [205,205,0];[100,205,0];[0,205,0];[0,188,199];[0,0,255];[122,0,255]  ]/255;
n=draw(1,1);

r1=[10,15,20,25,30,40,50,60,70,104,142,305];%10,15,20,25,30,
legendCell=cell(1,length(r1));
for i=1:length(r1)
    legendCell{i} = num2str(r1(i));
end

n1=1;n2=1;
i=1;
for j=2:1000
    %if draw(i,1)<draw()
    if j>=length(draw) || draw(j+1,1)<draw(j,1)      
        n1=n2;
        n2=j;
        xx1 = 1.2:0.05:min(13.5,0.9*r1(i));
        yy1 = spline(draw(n1+1:n2,1),draw(n1+1:n2,2),xx1);
        semilogy(xx1 ,yy1,'Color',color(i,:),'LineWidth',1);%
        i=i+1;
        hold on
    end
    j=j+1;
    if j>length(draw)
        break
    end 
end

plot(c40(:,1),0.01*10.^c40(:,2),'.','MarkerSize',13,'MarkerFace','k','MarkerEdge','k','LineWidth',2)
plot(c55(:,1),0.01*10.^c55(:,2),'^','MarkerSize',5,'MarkerFace','k','MarkerEdge','k','LineWidth',0.1)
plot(cbeard(:,1),0.01*10.^cbeard(:,2),'x','MarkerSize',8,'MarkerFace','k','MarkerEdge','k','LineWidth',2)
plot(ckinzer(:,1),0.01*10.^ckinzer(:,2),'s','MarkerSize',5,'MarkerFace','k','MarkerEdge','k','LineWidth',.1)
plot(cjonas(:,1),0.01*10.^cjonas(:,2),'o','MarkerSize',5,'MarkerEdge','k','LineWidth',.1)
for i=1:length(c40)
    text(c40(i,1)-0.1,0.01*10.^(c40(i,2)-0.05),'40')
end
for i=1:length(c55)
    text(c55(i,1)-0.1,0.01*10.^(c55(i,2)-0.05),'55')
end
for i=1:length(cbeard)
    text(cbeard(i,1)-0.1,0.01*10.^(cbeard(i,2)-0.05),num2str(cbeardradius(i)))
end
text(ckinzer(1,1)-0.1,0.01*10.^(ckinzer(1,2)-0.05),'73')
text(cjonas(1,1)-0.1,0.01*10.^(cjonas(1,2)-0.05),'40')
axis([0,14,0.007,1])
xlabel('r_2 (¦Ìm)')
ylabel('Collision Efficiency E')
legendCell{end+1} = 'Picknet(1960)';
legendCell{end+1} = 'Woods & Mason(1964)';
legendCell{end+1} = 'Beard & Pruppacher(1971)';
legendCell{end+1} = 'Kinzer & Cobb(1958)';
legendCell{end+1} = 'Jonas & Goldsmith(1972)';
legend(legendCell)
%legend('30','40','50','60','70','104','142','305','Picknet(1960)','Woods & Mason(1964)','Beard & Pruppacher(1971)','Kinzer & Cobb(1958)','Jonas & Goldsmith(1972)')