mat=load('D:\Fortran\collision_terminal\terminal_400 1.csv')./1e4;
x=1:1:41;
x=1.*2.^((x-1)/4);

figure(1)
subplot(2,1,1)
    cut=10
    loglog( x(cut:end),abs(mat(cut:end,1) ),'LineWidth',1)
    hold on
    loglog( x,abs(mat(:,15) ),'LineWidth',1)
    hold on
    loglog( x,abs(mat(:,7) ),'LineWidth',1)
    hold on
    axis([1,1000,1e-2,1e3])
    %set(gca,'xtick',[])
    ylabel('terminal velocity (cm s^{-1})')
    title('Electric field (400 V cm^{-1}) effect on terminal velocity')
    
    legend('q=-32e r^2 ','no charge','q=+32e r^2 ')
subplot(2,1,2)
    plot( x(1:cut+1),(mat(1:cut+1,1) ),'LineWidth',1)
    hold on
    axis([1,10,-1e-1,1e-2])
    yticks([-10^-1 -0.05 -10^-2  0 10^-2])
    %set(gca,'YDir','reverse');
xlabel('droplet radius r (¦Ìm)')
 ylabel('terminal velocity (cm s^{-1})')


%q=-32e (r/¦Ìm)^2