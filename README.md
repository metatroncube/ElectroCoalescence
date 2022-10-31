# ElectroCoalescence
The enhancement of droplet collision by electric charges and atmospheric electric fields

I have to admit that my program codes may be poor in readability.
The folder "collision" and "stochastic" are Fortran programs.
The folder "MATLAB" consists of the MATLAB programs used for plotting those figures.


(A)
".../collision/main.f95" is for collision efficiencies. It takes several days to calculate all the collision efficencies and kernels on a normal computer with linux system.

Ternimal velocities have been computed in advance, and saved in "terminal_000.csv", "terminal_200.csv" and "terminal_400.csv".

When computing collision efficencies and kernels with different electric field, please rename above csv files to "terminal.csv".Then, set the parameter e0 in the 4th line in main.f95:

for field = 0 V/m, please set
【Line 4】real :: ... e0=-160*( 00000 ) ...
【Line 83】open(unit=512,file="terminal_000.csv",status="old")

for field = 200 V/m, please set
【Line 4】real :: ... e0=-160*( 20000 ) ...
【Line 83】open(unit=512,file="terminal_200.csv",status="old")

for field = 400 V/m, please set
【Line 4】real :: ... e0=-160*( 40000 ) ...
【Line 83】open(unit=512,file="terminal_400.csv",status="old")

(nondimensionalize unit"1": time= 1 second, length= 1μm, charge= 1 elementary charge e, mass= 1×10^-15kg, force= 1×10^-21N, field= 1×6.25×10^-3 V/m...)

After computing, the output is in folder "matrix", where "matrix.csv" is collision efficiencies and "matrixkernel.csv" is collision kernels. They should be moved to ".../stochastic/E000", ".../stochastic/E200", or ".../stochastic/E400" accordingly.


I should declare that some of  program codes in ".../collision/main.f95" is redundant. Because I have tried different methods to compute droplet motions (like Hamielec 1962 ,and Beard 1976...), and finally decide to choose the former one, since it fits relatively well with all the sizes of droplets with different Reynolds numbers (when there's no electric charge). Other methods have been abandoned, while they can still be activated by changing some parameters.

(B)
".../stochastic/main.f95" is for solving the stochastic collision equation (SCE),based on the efficiencies & kernels already derived.
Coalescence efficiencies are got in in advance, and saved in "coal_efficiency_beard.csv".

You can just type following text in windows cmd, under the folder ".../stochastic/" :
-------------------------------------------------
gfortran main.f95 -o main.exe
main
1 1
main
1 2
main
1 3
main
1 4
main
2 1
main
2 2
main
2 3
main
2 4
main
3 1
main
3 2
main
3 3
main
3 4
______________________
where the first number means initial average radius r=15, 9, 6.5μm
and the second number means different initial electric condition: 1=uncharged cloud, 2=charged without field, 3=charged under 200V/m, 4=charged under 400V/m.
All the results of cloud spectrum will be output into these folders named as "30", "60", "120"

(C)
Then you can use MATLAB to plot figures based on these data. 
For example, in ".../MATLAB/Spectrum_Comparision.m":
______________________
totaltime=120;
path=strcat('D:\...\Fortran\stochastic\',num2str(totaltime),'\spectrum')
...
______________________
you can choose "30", "60" or "120", to generate Figure 7,8,10 in the paper separately.

It's important that paths in those MATLAB files, such as "('D:\...\stochastic\')", should be changed manually to fit in with your computer!
