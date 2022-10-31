clear;close all;clc
r0='60'
casee=[1 2 3 4];sequence=5;
figure(1)
myjet=[
        0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0];

for cas=casee
    

loadname=strcat('D:\Fortran\stochastic\',r0,'\spectrum',mat2str(cas),'.txt');
loading=load(loadname);
y2=loading(sequence,:);
Z=zeros(15,38);
scale=38;
for i=0:14
    Z(i+1,:)=y2(1+i*scale:scale+i*scale);      
end
x=2*2.^((0:1:37)/4);
xxishu=2*2.^((0:9));
%x=0:1:37;
ybins=-7:1:7;%
y=[-32,-16,-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8,16,32];
hAxes = gca;
subplot(ceil(length(casee)/2),min(2,length(casee)),1+cas-casee(1))
imagesc(x,y,Z.*1e-6)

labs = sprintfc('%d',xxishu);
set(gca,'XTick',linspace(2,1197,10),'XTickLabels',labs);
labs = sprintfc('%d',y(2:2:14));
set(gca,'YTick',linspace(-27,27,7),'YTickLabels',labs);
colormap(myjet)
colorbar
set(gca,'FontSize',13);
xlabel('radius bins, r (¦Ìm)')
ylabel('charge bins, q/r^2 (e ¦Ìm^{-2})')
charged=' charged droplets under';
if cas==1
    field=' no field';
    charged=' uncharged droplets under';
elseif cas==2
    field=' no field';
elseif cas==3
    field=' field 200V m^{-1}';
elseif cas==4
    field=' field 400V m^{-1}';
end
title1=strcat('(',char(96+cas-casee(1)+1),') ',charged,field);
title2=strcat( 'mass distribution (g m^{-3})  after ',[' ' mat2str((sequence-1)*15)],' min');
title({title1,title2})

end