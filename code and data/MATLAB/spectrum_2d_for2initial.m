clear;close all;clc
casee=[1 2];sequence=1;
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
    

loadname=strcat('D:\´ó¶þ\Fortran\stochastic\60\spectrum',mat2str(cas),'.txt');
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
subplot(ceil(length(casee)),min(1,length(casee)),1+cas-casee(1))
imagesc(x,y,Z.*1e-6)

labs = sprintfc('%d',xxishu);
set(gca,'XTick',linspace(2,1197,10),'XTickLabels',labs);
labs = sprintfc('%d',y(2:2:14));
set(gca,'YTick',linspace(-27,27,7),'YTickLabels',labs);
colormap(myjet)
colorbar

xlabel('radius bins, r (¦Ìm)')
ylabel('charge bins, q/r^2 (e ¦Ìm^{-2})')
charged=' charged cloud spectrum';
if cas==1
    field=' no field';
    charged=' uncharged cloud spectrum';
end
title1=strcat('(',char(96+cas-casee(1)+1),') ',charged);
title2=strcat( ' initial mass distribution (g m^{-3})');
title(strcat(title1,title2))

end