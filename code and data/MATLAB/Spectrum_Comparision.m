clear;close all;clc
num=4;
totaltime=120;
path=strcat('D:\Fortran\stochastic\',num2str(totaltime),'\spectrum')

logy=1;Eng=1;Eng=Eng+1;
load ('D:\大二\Fortran\collision\scales.txt');
blockqs=2*scales(5)+1;
scale=scales(2)+1;
factor=4/log(2)*1e-6;
width=1;
color=[[0,0,255];[255,0,0]; [0,195,0];[255,0,255];[0,0,0];[122,122,122]]/255;
%color=color.*0+1
r=2*2.^((0:1:37)/4);
m=4/3*pi*r.^3;
rbar=9;
mbar=4/3*pi*rbar^3;

factor2=factor*(1./m)*1e9;


figure(1)

title({'Mean radius 15μm, evolution time 20min';'Spectrum mass distribution'})
name=strcat(path,num2str(1),'.txt');
%name=strcat(name,);
spectrum0=load (name);       
panels=size(spectrum0,1);   
spec=zeros(num,panels,scale);
for n=1:num

    name=strcat(path,num2str(n));
    name=strcat(name,'.txt');
    spectrum0=load (name);    
    numbers=zeros(blockqs,scale);
    x=1:1:scale;
    xx=scales(1)*2.^((x-1)/4);
    y=1:1:blockqs;
    z=zeros(panels,scale);
    %z2=zeros(1,scale);
    
    add=zeros(panels,scale);
    %y1=spectrum0(1,:)./1e6;
    %y2=spectrum0(2,:)./1e6;
    for i=0:blockqs-1
        add=spectrum0(:,1+i*scale:scale+i*scale);
        z=z+add;
        %numbers(i+1,:)=numbers(i+1,:)+add;
    end
    spec(n,:,:)=z;
end


hold on
for k=1:panels-1

subplot(4,2,2*k-1)
t=k*totaltime/(panels-1);
    semilogx(xx,factor.*z(1,:),':','LineWidth',width*0.75,'Color',[0 0 0])   
    axis( [ 1,1200,0,1.6])
set(gca,'XTick',[1,1e1,1e2,1e3]);
set(gca,'YTick',[0,0.5,1,1.5]);
if k==1
    title(strcat('(a) spectrum after',32,num2str(t),' min'))
elseif k==2
    title(strcat('(b) spectrum after',32,num2str(t),' min'))
elseif k==3
    title(strcat('(c) spectrum after',32,num2str(t),' min'))
elseif k==4
    title(strcat('(d) spectrum after',32,num2str(t),' min'))
end
ylabel('M(lnr)  (g m^{-3})')
if k==4
    xlabel(' r (μm)')
    
end
hold on
subplot(4,2,2*k)
    loglog(xx,factor2.*z(1,:),':','LineWidth',width*0.75,'Color',[0 0 0]) 
    axis( [ 1,1200,1e-6,1e6])
set(gca,'XTick',[1,1e1,1e2,1e3]);
set(gca,'YTick',[1e-5,1,1e5]);
title(strcat('spectrum after',32,num2str(t),' min'))
ylabel('n(lnr) ')

if k==4
    xlabel(' r (μm)')
    
end
hold on
for n=1:num
	subplot(4,2,2*k-1)

        semilogx(xx,factor.*reshape(spec(n,k+1,:),[1,38]),'LineWidth',width,'Color',color(n,:))
%     if k==1 && n==num
%         word=legend('t=0','No charge','Charged','Charged+ Field 200V/cm','Charged+ Field 400V/cm');
%         set(word,'FontSize',9)
%     end
    if k==1 && n==num
        word=legend('t=0','No charge','Charged','Charged, 200V cm^{-1}','Charged, 400V cm^{-1}','Location','northeast');
        set(word,'FontSize',6)
    end
    hold on
    subplot(4,2,2*k)
        loglog(xx,factor2.*reshape(spec(n,k+1,:),[1,38]),'LineWidth',width,'Color',color(n,:))


    hold on
    set(gcf,'unit','normalized','position',[0.2,-0.0,0.6,2.8]);
end
end
savename=strcat('D:\大二\本科生科研\图片\',num2str(totaltime),'.png')
saveas(gcf,[savename])

