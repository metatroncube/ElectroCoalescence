
record=load ('D:\Fortran\collision\record.txt');
%record(:,2)=-record(:,2);
r1=record(1,1);
r2=record(1,2)+record(1,1);
[hang,lie]=size(record);
p=1;
q=1;
for i=2:(hang-1)
    if abs(record(i+1,2)-record(i,2))>=(r1+r2) || i==hang-1
        p=q+2;
        q=i-1;
        plot(record(p:q,1),-record(p:q,2))
        hold on
    end
end




t=0:pi/100:2*pi;
x=r1*sin(t);y=r1*cos(t);
plot(x,y,'k','LineWidth',1.2)
hold on
x=r2*sin(t);y=r2*cos(t);
plot(x,y,':k','LineWidth',0.3)
hold on

y=-700:10:+290;
x=y.*0;
plot(x,y,':k','LineWidth',0.2)
axis equal
hold off
%delete('D:\´ó¶þ\Fortran\collision\record.txt')