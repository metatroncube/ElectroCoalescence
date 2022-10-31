Module pub
real :: dt,characteristictime,k=230093,lambda0=0.066,g=9800000,pi=3.1415926,viscosity=1.718E4,airdensity=0.00129
!      步长  特征时间    k单位化静电常数  分子自由程   重力加速度      圆周率      空气粘性
real :: temp=283, pressure=900,watervis,e0=-160*( 20000 ),sig=-1,wake=8,area=10!!!!!
   !      温度           压强
real::a10=-2.029,a11=0.8222,a12=-0.02253,a20=-1.666,a21=0.632,a30=-0.91,a31=0.2619,a32=0.04214! Beard 1976
real::AA1,AA2,FF1,FF2,bb1,bb2,yy1,yy2,re11,re22
real::g1,g2,g3,g4,g5,g6,delta,o1,o2,alp1,alp2,qq,f01,f02,sint,cost,ep1,ep2 ,dist=30,sig2=1!1973流
real A1,A2,A3,A4,B1,B2,B3,B4!flowcorrect2的方法 1962Canadian
integer :: firstordercorrect=1,relative=1,renoldcorrect=1,rungekutta=0,exactelectric=1,outfield=1,ignore=1
      !    是否一阶修正    是否相对位移平移 对stokes流修正 龙格库塔法  电荷严格解       外电场
integer::relaflow=1,flowcorrect=2,pruppacher=1,testingtime=800000,l=0
            !相对流    流场修正   随时计算拖曳系数 开始显示       计算回合数
integer :: flowinitialization=0, txtoutput=0 ,matrixoutput=1,screenoutput=1, autorefine=1,formatrix=1
         !初速度考虑流场     是否向文件输出轨迹 ,是否屏幕输出每次判断, 自动细化
integer:: electric_terminal=1 !0末速度不考虑电场 1末速度考虑电场
real f_low(24),f_up(10),k_low(10),k_up(10,24),beta(10)!1964davis 电荷严格解
real:: draw_right=28,draw_numbers=20
type drop
    real mass,r,q, x,y, vx,vy, ax,ay,vt,density ,c1,c2   ,renold
    !  质量半径电荷位置，速度，加速度， 密度，  修正系数  雷诺数
end type
real:: ini_r1=35,ini_r2=4, &
ini_q1=-0*(30**2),ini_q2=0*(12.6**2)

End Module pub

program main
  use pub
  implicit none
  real,external:: s_fun,t_fun,u_fun,ss_fun
  integer ,external:: near
  integer :: blocks_q=7,blocks_r=37,orders
  integer :: updown=1! 大小滴那个在上哪个在下
  real:: r0=2,q0=0.5,log_dr=1.5/6.0,log_dq=1.0,qmaxindex=2
  integer :: collision=0,fail=0 !判断指标
  real:: x,y,t,r,f=0,fex=0,fey=0,u1x=0,u1y=0,u2x=0,u2y=0
  !       位置      电作用力    induced流场,空气密度
  real X1K1,X1K2,X1K3,X1K4,X1L1,X1L2,X1L3,X1L4
  real Y1K1,Y1K2,Y1K3,Y1K4,Y1L1,Y1L2,Y1L3,Y1L4
  real X2K1,X2K2,X2K3,X2K4,X2L1,X2L2,X2L3,X2L4
  real Y2K1,Y2K2,Y2K3,Y2K4,Y2L1,Y2L2,Y2L3,Y2L4
  real stokesnumber ,lambda,dividing
  integer :: i=0,j=0,p1=11,p2=10,q1=7,q2=1,hori,vert
  real:: efficiency=1
  real targetradius,rmax,rmin !二分法
  real:: t1=1,t2=1,rlast=0,rlast2=0
  type(drop) drop1,drop2,drop1correct,drop2correct,drop1save,drop2save

  real,allocatable::scale_radius(:),scale_charge(:,:),matrix_effi(:,:),matrix_kernel(:,:),terminal(:,:)
  integer,allocatable::register(:,:)
  !以下三行仅供本研展示!!
  real::ranges(1)
  integer:: rangenow=1

  ranges=(/31.4/)!0,40,50,60,70,104,142,305/)
  !!!!!!!!!!!!!!
  f_low=0.0
  f_up=0.0
  k_low=0.0
  k_up=0.0
  beta=0.0
  beta(4)=-1.0
  beta(10)=-1.0
  if(outfield.eq.0)then
    e0=0
  end if

  if(formatrix.eq.1 .or. electric_terminal.eq.1)then
        open(unit=666,file="scales.txt",status="old")
        read(666,*),r0,blocks_r,log_dr,q0,blocks_q,log_dq,qmaxindex ,rlast2,rlast2,rlast2,rlast2,rlast2,rlast2
        log_dr=log_dr/3.0
  end if



  allocate( scale_radius(blocks_r),scale_charge(blocks_r,-blocks_q:blocks_q) )
  orders=blocks_r*(2*blocks_q+1)
  allocate( matrix_effi(orders,orders),matrix_kernel(orders,orders),register(orders,orders))

  if(electric_terminal.eq.1)then
    allocate( terminal(blocks_r,-blocks_q:blocks_q)  )
    open(unit=512,file="terminal_200.csv",status="old")! 这是末速度的文件
    do i=1,blocks_r,1
        read(512,*),terminal(i,:)
    end do
  end if

  scale_radius(1)=r0

  do i=2,blocks_r,1
    scale_radius(i)=scale_radius(i-1)*2**log_dr
  end do

  scale_charge(:,0)=0
  if(blocks_q.gt.0)then
  scale_charge(:,1)=q0*(scale_radius)**qmaxindex
  scale_charge(:,-1)=-q0*(scale_radius)**qmaxindex
  end if
  if(blocks_q.gt.1)then
  do i=2,blocks_q,1
    scale_charge(:,i)=scale_charge(:,i-1)*2**log_dq
    scale_charge(:,-i)=-scale_charge(:,i)
  end do
  end if

  matrix_effi=0
  matrix_kernel=0
  register=0
  if(txtoutput.eq.1)then
  open(unit=999,file="record.txt",status="replace")
  end if
  if(formatrix.eq.2)then
    open(unit=111,file="draw.csv",position='append')
    write(111,*) draw_numbers,draw_numbers
  end if
  if(matrixoutput.eq.1 .and. formatrix.eq.1)then
  open(unit=888,file="./matrix/matrix.csv",status="replace")
  open(unit=889,file="./matrix/matrixkernel.csv",status="replace")
  open(unit=777,file="./matrix/scales.txt",status="replace")
  open(unit=555,file="./matrix/temp.csv")
  write(777,*),r0,",",blocks_r,",",log_dr,",",q0,",",blocks_q,",",log_dq,","
  end if
  viscosity=viscosity*(temp/273.15)**1.5*(124+273.15)/(124+temp)
  airdensity=airdensity*(pressure/1013.25)*273.15/temp
  t=0

  lambda=lambda0*1013.25/pressure*temp/273.15
  drop1=drop(1, 1, 0,  0,   0,    0,0,   0,0,   0, 1     , 1.0 ,1.0,0)
  drop2=drop(1, 1, 0, 10,900,    0,0,   0,0,   0, 1     , 1.0 ,1.0,0)
!  初始化       r  q    x    y    vx vy  ax ay  vt density c         re
  drop1%r=ini_r1
  drop2%r=ini_r2
  drop1%q=ini_q1
  drop2%q=ini_q2
  drop1%density=1
  drop2%density=1

   call distance(drop1,drop2,x,y,r)
  if(formatrix.eq.1) then
    p1=2
    p2=1
    q1=-7
    q2=-7
  elseif(electric_terminal.eq.1)then
    read(*,*),p1,p2,q1,q2
  end if

4 continue ! matrix循环
if(matrixoutput.eq.1 .and. formatrix.eq.1)then
    close(555)
    open(unit=555,file="./matrix/temp.csv",status="old",position='append')
end if
    !print*,"p1=",p1,"q1=",q1,"p2=",p2,"q2=",q2
    vert=p2*(2*blocks_q+1)-blocks_q+q2
    hori=p1*(2*blocks_q+1)-blocks_q+q1

  if(register(hori-2*q1,vert-2*q2).eq.1 .and. outfield.eq.0)then
    matrix_effi(hori,vert)=matrix_effi(hori-2*q1,vert-2*q2)
    matrix_kernel(hori,vert)=matrix_kernel(hori-2*q1,vert-2*q2)
    goto 18
  end if

  efficiency=1.0
  l=0
  i=0
  dividing=300

  if(formatrix.eq.1 .or. electric_terminal.eq.1)then
  !
  drop1%r=scale_radius(p1)
  drop2%r=scale_radius(p2)
  !
  drop1%q=scale_charge(p1,q1)
  drop2%q=scale_charge(p2,q2)

  !drop1%q=sign(28800,q1)!
  end if



  if( (drop1%r+drop2%r).lt.5)then
    dividing=1000
  end if
 if(txtoutput.eq.1 .and.  formatrix.eq.0)then
   write (999,*) drop1%r,drop2%r
  end if
  rmax=(drop1%r+drop2%r)*5
  rmin=0
targetradius=rmax/15.0
  call givemass(drop1)
  call givemass(drop2)


    drop1%c1=1.0/(1+1.26*lambda/drop1%r)

    call caculatec2(drop1)

    drop2%c1=1.0/(1+1.26*lambda/drop2%r)

    call caculatec2(drop2)
if(drop1%r.gt.10)drop1%c1=1.0
if(drop2%r.gt.10)drop2%c1=1.0
if(formatrix.eq.0) print*,"c1=",drop1%c1,drop2%c1
  if(drop2%r.le.20)then
    drop2%c2=1.0
  end if

  if(renoldcorrect.eq.0)then
    drop1%c2=1
    drop2%c2=1
  end if
do while( (rmax-rmin).gt.(targetradius/1000.0) .and. efficiency.gt.1e-12 .and. l.lt.testingtime) !不同入射半径，二分法循环
    i=i+1
!再次初始化
  fail=0
  collision=0
  drop1%x=0
  drop2%x=targetradius*sig2
  drop1%y=0
  drop2%y=dist*(drop2%r+drop1%r)
  drop1%vx=0
  drop1%vy=0
  drop2%vx=0
  drop2%vy=0
  drop1%ax=0
  drop1%ay=0
  drop2%ax=0
  drop2%ay=0
  drop1save=drop1
  drop2save=drop2

2 continue  !对发散的重新细化改正
  t=0

  drop1=drop1save
  drop2=drop2save
  drop1%vt=2.0*drop1%density*g*(1-airdensity/drop1%density)*drop1%r**2.0/9.0/viscosity/drop1%c1/drop1%c2
  drop2%vt=2.0*drop2%density*g*(1-airdensity/drop2%density)*drop2%r**2.0/9.0/viscosity/drop2%c1/drop2%c2

  if(electric_terminal.eq.1)then
    if(e0.eq.0 .or. outfield.eq.0)then
        !print*,terminal(p1,0)/drop1%vt
        drop1%vt=terminal(p1,0)
        drop2%vt=terminal(p2,0)
    else if(e0.lt.0)then
        drop1%vt=terminal(p1,q1)
        drop2%vt=terminal(p2,q2)
    else
        drop1%vt=terminal(p1,-q1)
        drop2%vt=terminal(p2,-q2)
    end if
  end if

  updown=1
  if (drop1%vt.lt.drop2%vt)then
        updown=-1
        drop2%y=-drop2%y
  end if

  call distance(drop1,drop2,x,y,r)
  drop1%vy=drop1%vt
  drop2%vy=drop2%vt
  u1x=0
  u1y=0
  u2x=0
  u2y=0
if(flowcorrect.eq.3)then
    flowcorrect=2
    call stokesflow(drop1,drop2,u1x,u1y,u2x,u2y,x,y,r)
    flowcorrect=3
else
    call stokesflow(drop1,drop2,u1x,u1y,u2x,u2y,x,y,r)
end if
fex=0
fey=0


if(flowinitialization.eq.1 )then
  drop2%vx=drop2%vx+u1x
  drop2%vy=drop2%vy+u1y
  drop1%vx=drop1%vx+u2x
  drop1%vy=drop1%vy+u2y
  end if

  stokesnumber=(drop1%vy-drop2%vy)*drop2%vt/g/2.0/drop1%r


  do
    rlast2=rlast
    rlast=r

  if (r.le.(drop1%r+1.0*drop2%r) ) then
    collision=1
    goto 10
    end if


    if ( updown*drop2%y.lt.(updown*drop1%y-wake*drop1%r-wake*drop2%r) .or. r.gt.dist*area*(drop2%r+drop1%r)) then
    fail=1
    goto 11
  end if

   call distance(drop1,drop2,x,y,r)

  call electricforce(drop1,drop2,f,fex,fey)

   if(flowcorrect.ne.3) call stokesflow(drop1,drop2,u1x,u1y,u2x,u2y,x,y,r)


  call accelerate (drop1,drop2,fex,fey,u1x,u1y,u2x,u2y)


  !t1=sqrt( (drop1%vy-drop2%vy)**2/(drop2%ax**2+drop2%ay**2))/dividing
  !t2=drop2%vt**2/g
  !dt=min(t1,t2)
  !print*,t1,t2
  t1=sqrt( (drop2%vy**2+drop2%vx**2)/(drop2%ax**2+drop2%ay**2))
  t2=abs( r / sqrt( (drop1%vy-drop2%vy)**2+(drop1%vx-drop2%vx)**2 ) )
  characteristictime=min(t1,t2)
  dt=characteristictime/dividing

    l=l+1


  if(firstordercorrect.eq.0 .and. rungekutta.eq.0) then
    call moving(drop1,drop2)

  else if (rungekutta.eq.0) then
    drop1correct=drop1
    drop2correct=drop2
    call virtualmoving(drop1,drop2)
    call distance(drop1,drop2,x,y,r)
    call electricforce(drop1,drop2,f,fex,fey)
     if(flowcorrect.ne.3) call stokesflow(drop1,drop2,u1x,u1y,u2x,u2y,x,y,r)
    call accelerate (drop1,drop2,fex,fey,u1x,u1y,u2x,u2y)
    drop1correct%ax=drop1%ax
    drop1correct%ay=drop1%ay
    drop2correct%ax=drop2%ax
    drop2correct%ay=drop2%ay
    drop1=drop1correct
    drop2=drop2correct
    call moving(drop1,drop2)

  else ! Runge-Kutta大法好！
    drop1correct=drop1
    drop2correct=drop2
    x1k1=drop1%vx
    x2k1=drop2%vx
    y1k1=drop1%vy
    y2k1=drop2%vy
    x1L1=drop1%ax
    x2L1=drop2%ax
    y1L1=drop1%ay
    y2L1=drop2%ay
    !一阶
    x1k2=drop1%vx+0.5*dt*x1L1
    x2k2=drop2%vx+0.5*dt*x2L1
    y1k2=drop1%vy+0.5*dt*y1L1
    y2k2=drop2%vy+0.5*dt*y2L1
    call teleport(drop1correct,drop2correct,0.5*dt*x1k1,0.5*dt*y1k1,0.5*dt*x2k1,0.5*dt*y2k1)
    drop1correct%vx=x1k2
    drop2correct%vx=x2k2
    drop1correct%vy=y1k2
    drop2correct%vy=y2k2
    call distance(drop1correct,drop2correct,x,y,r)
    call electricforce(drop1correct,drop2correct,f,fex,fey)
    if(flowcorrect.ne.3)  call stokesflow(drop1correct,drop2correct,u1x,u1y,u2x,u2y,x,y,r)
    call accelerate(drop1correct,drop2correct,fex,fey,u1x,u1y,u2x,u2y)
    x1L2=drop1correct%ax
    x2L2=drop2correct%ax
    y1L2=drop1correct%ay
    y2L2=drop2correct%ay
    !二阶
    drop1correct=drop1
    drop2correct=drop2
    x1k3=drop1%vx+0.5*dt*x1L2
    x2k3=drop2%vx+0.5*dt*x2L2
    y1k3=drop1%vy+0.5*dt*y1L2
    y2k3=drop2%vy+0.5*dt*y2L2
    call teleport(drop1correct,drop2correct,0.5*dt*x1k2,0.5*dt*y1k2,0.5*dt*x2k2,0.5*dt*y2k2)
    drop1correct%vx=x1k3
    drop2correct%vx=x2k3
    drop1correct%vy=y1k3
    drop2correct%vy=y2k3
    call distance(drop1correct,drop2correct,x,y,r)
    call electricforce(drop1correct,drop2correct,f,fex,fey)
   if(flowcorrect.ne.3)   call stokesflow(drop1correct,drop2correct,u1x,u1y,u2x,u2y,x,y,r)
    call accelerate(drop1correct,drop2correct,fex,fey,u1x,u1y,u2x,u2y)
    x1L3=drop1correct%ax
    x2L3=drop2correct%ax
    y1L3=drop1correct%ay
    y2L3=drop2correct%ay
    !三阶
    drop1correct=drop1
    drop2correct=drop2
    x1k4=drop1%vx+dt*x1L3
    x2k4=drop2%vx+dt*x2L3
    y1k4=drop1%vy+dt*y1L3
    y2k4=drop2%vy+dt*y2L3
    call teleport(drop1correct,drop2correct,dt*x1k3,dt*y1k3,dt*x2k3,dt*y2k3)
    drop1correct%vx=x1k4
    drop2correct%vx=x2k4
    drop1correct%vy=y1k4
    drop2correct%vy=y2k4
    call distance(drop1correct,drop2correct,x,y,r)
    call electricforce(drop1correct,drop2correct,f,fex,fey)
     if(flowcorrect.ne.3) call stokesflow(drop1correct,drop2correct,u1x,u1y,u2x,u2y,x,y,r)
    call accelerate(drop1correct,drop2correct,fex,fey,u1x,u1y,u2x,u2y)
    x1L4=drop1correct%ax
    x2L4=drop2correct%ax
    y1L4=drop1correct%ay
    y2L4=drop2correct%ay
    !四阶
    !结算！
    drop1%x=drop1%x+dt/6*(x1k1+2*x1k2+2*x1k3+x1k4)
    drop1%y=drop1%y+dt/6*(y1k1+2*y1k2+2*y1k3+y1k4)
    drop2%x=drop2%x+dt/6*(x2k1+2*x2k2+2*x2k3+x2k4)
    drop2%y=drop2%y+dt/6*(y2k1+2*y2k2+2*y2k3+y2k4)

    drop1%vx=drop1%vx+dt/6*(x1L1+2*x1L2+2*x1L3+x1L4)
    drop1%vy=drop1%vy+dt/6*(y1L1+2*y1L2+2*y1L3+y1L4)
    drop2%vx=drop2%vx+dt/6*(x2L1+2*x2L2+2*x2L3+x2L4)
    drop2%vy=drop2%vy+dt/6*(y2L1+2*y2L2+2*y2L3+y2L4)

    end if



  if (relative.eq.1)then
  drop2%x=drop2%x-drop1%x
  drop2%y=drop2%y-drop1%y
  drop1%x=0
  drop1%y=0
  end if

 if(txtoutput.eq.1 .and. formatrix.eq.0)then
 write (999,*) drop2%x,drop2%y
 end if


    x=drop2%x-drop1%x
    y=drop2%y-drop1%y
    r=sqrt(x**2+y**2)

if (l.ge.testingtime+1)then
    print*,"run time error!","p1=",p1,"q1=",q1,"p2=",p2,"q2=",q2
    exit
end if
if( (isnan(drop2%ax) .or. isnan(drop2%ay))  )then

    if(rlast2.lt.1.2*(drop1%r+drop2%r) .or. rlast.lt.1.2*(drop1%r+drop2%r))then
        collision=1
   !
        exit
    else
        efficiency=-100
        goto 19
    end if

end if
    if(abs(drop2%vy-drop1%vy) .gt. abs(1000000*(drop2%vt+drop1%vt)) .and. autorefine.eq.1 ) then
        dividing=dividing*1.5
        if(screenoutput.eq.1)then
            print*,abs(drop2%vy-drop1%vy),abs(1000000*(drop2%vt+drop1%vt)),"unstable! auto-refine..."
            print*,"p1=",p1,"q1=",q1,"p2=",p2,"q2=",q2
            print*,"drop1 r,q=",drop1%r,drop1%q,"drop2 r,q=",drop2%r,drop2%q
            print*,"vt, v=",drop1%vt,drop1%vy,"2vt, v=",drop2%vt,drop2%vy
        end if
        goto 2
    end if


  end do
  10 continue
  if (collision.eq.1 .and. screenoutput.eq.1)then
  print *,"collision","totaltimes=",l
  end if
  11 continue
  if (fail.eq.1 .and. screenoutput.eq.1)then
  print *,"fail","totaltimes=",l
  end if

if(fail.eq.1)then
    rmax=targetradius
end if
if(collision.eq.1)then
    rmin=targetradius
end if
  targetradius=(rmax+rmin)/2
efficiency=(targetradius**2) /( (drop1%r+drop2%r)**2  )
end do

if(efficiency.le.1e-12 .or.  l.ge.testingtime)then
    efficiency=0
    targetradius=0
end if
19 continue
if(formatrix.ne.1)then
print*,"p1=",p1,"q1=",q1,"p2=",p2,"q2=",q2
print*,"drop1 r,q=",drop1%r,drop1%q,"drop2 r,q=",drop2%r,drop2%q
print*,"vt, v=",drop1%vt,drop1%vy,"2vt, v=",drop2%vt,drop2%vy
print*,"c2=",drop1%c2,drop2%c2
print*,"rounds and times=",i,l
print*,"efficiency=",efficiency
print*
end if
 close(999)

if(formatrix.eq.2)then
    write(111,*) drop2%r,efficiency
    drop2%r=drop2%r+(draw_right - ini_r2)/draw_numbers
    drop2%q=ini_q2*(drop2%r/ini_r2)**2
    if(drop2%r.gt.draw_right)then
        drop2%r=ini_r2
        drop2%q=ini_q2
        rangenow=rangenow+1
        if(rangenow.gt.size(ranges)) stop
        drop1%r=ranges(rangenow)
    end if
    goto 4
end if


if(formatrix.ne.1) goto 5

matrix_effi(hori,vert)=efficiency
matrix_kernel(hori,vert)=pi*targetradius**2*abs(drop1%vt-drop2%vt)
18 continue

if(blocks_r.ge.2 .and. ignore.eq.1 .and. q2.eq.(-blocks_r+1) .and. near( matrix_effi(hori,vert),matrix_effi(hori,vert-1)).eq.1)then
  if(efficiency.gt.1e-6)then
    matrix_effi(hori,vert-1 : vert-1+2*blocks_q)=efficiency
    matrix_kernel(hori,vert-1 : vert-1+2*blocks_q)=matrix_kernel(hori,vert)
    register(hori,vert-1 : vert-1+2*blocks_q)=1
    q2=blocks_q
  end if
end if

 write(555,*),p1,",",q1,",",p2,",",q2,",",matrix_kernel(hori,vert)

 register(hori,vert)=1
    if(p1.eq.blocks_r .and. p2.eq.blocks_r-1 .and. q1.eq.blocks_q.and.q2.eq.blocks_q)then
        goto 5
    end if

    if(q2.lt.blocks_q)then
        q2=q2+1
        else if((p2+1).lt.p1)then
            q2=-blocks_q
            p2=p2+1
        else if((p2+1).eq.p1 .and. q1.lt.blocks_q)then
            q1=q1+1
            q2=-blocks_q
            p2=1
        else if(q1.eq.blocks_q)then
            p1=p1+1
            p2=1
            q1=-blocks_q
            q2=-blocks_q
        end if
        goto 4
5 continue
print*,"ended"

!7 format(1x,f6.4,",")
8 format(1x,e10.4,",")
if(formatrix.eq.1 .and. matrixoutput.eq.0)then
    do i=1,orders,1

        do j=1,orders,1
        write(*,8,advance='no'),matrix_effi(j,i)
        end do
     write(*,*)

    end do
   do i=1,orders,1
        do j=1,orders,1
        write(*,8),matrix_kernel(j,i)
        end do

    end do
end if
if(matrixoutput.eq.1 .and. formatrix.eq.1)then
    do i=1,orders,1


        do j=1,orders,1
        write(888,8,advance='no'),matrix_effi(j,i)
        end do
     write(888,*)
    end do
   do i=1,orders,1
        do j=1,orders,1
        write(889,8),matrix_kernel(j,i)
        end do

    end do

end if
close(777)
close(888)
close(889)
close(555)


end


subroutine distance(drop1,drop2,x,y,r) !给出相对位置信息
    use pub
    real x,y,r
    type(drop) drop1,drop2
    x=drop2%x-drop1%x
    y=drop2%y-drop1%y
    r=sqrt(x**2+y**2)
end subroutine

subroutine electricforce(drop1,drop2,f,fex,fey) !计算电相互作用
    use pub
    real x,y,r,f1,f2,f3,s1,s2,f,fex,fey
    type(drop) drop1,drop2
    x=drop2%x-drop1%x
    y=drop2%y-drop1%y
    r=sqrt(x**2+y**2)
    if(abs(drop1%q).le.1e-6 .and. abs(drop2%q).le.1e-6 .and. outfield.eq.0)then
        f=0
        fex=0
        fey=0
    elseif(exactelectric.eq.0)then
    f1=drop1%q*drop2%q/r**2
    s1=drop2%r*drop1%q**2*( 1.0/r**3-r/(r**2-drop2%r**2)**2 )
    s2=drop1%r*drop2%q**2*( 1.0/r**3-r/(r**2-drop1%r**2)**2 )
    f3=1.0/r**4+1.0/(r**2-drop1%r**2-drop2%r**2)**2 -1.0/(r**2-drop2%r**2)**2-1.0/(r**2-drop1%r**2)**2
    f2=f3*drop1%q*drop2%q*drop2%r*drop1%r
    f=k*(f1+f2+s1+s2)

    fex=f*x/r
    fey=f*y/r
    else
    call electric_accurate(drop1,drop2,f,fex,fey)
    end if



end subroutine

subroutine givemass(drop1) !按照半径赋予质量
    use pub
    type(drop) drop1
    drop1%mass=4.0*(drop1%r**3)*pi*drop1%density/3.0

end subroutine
subroutine accelerate (drop1,drop2,fex,fey,u1x,u1y,u2x,u2y) !瞬时加速度
    use pub
    type(drop) drop1,drop2
    real fex,fey,flow1x,flow1y,flow2x,flow2y,u1x,u1y,u2x,u2y,v1,v2,u1,u2,x,y,r,re1,re2,lg

if(pruppacher.eq.1 .and. flowcorrect.ne.3)then
    v1=sqrt((drop1%vx-u2x)**2+(drop1%vy-u2y)**2)
    v2=sqrt((drop2%vx-u1x)**2+(drop2%vy-u1y)**2)
    re1=2*airdensity/viscosity*v1*drop1%r
    re2=2*airdensity/viscosity*v2*drop2%r

    lg=log(re1)
    if(re1.lt.20)then
        flow1x=6*pi*viscosity*drop1%r*v1*(1+exp(a10+a11*lg+a12*lg**2))/drop1%mass
    elseif(re1.lt.258)then
        flow1x=6*pi*viscosity*drop1%r*v1*(1+exp(a20+a21*lg))/drop1%mass
    else
        flow1x=6*pi*viscosity*drop1%r*v1*(1+exp(a30+a31*lg+a32*lg**2))/drop1%mass
    end if
    flow1y=flow1x*(drop1%vy-u2y)/v1
    flow1x=flow1x*(drop1%vx-u2x)/v1

    lg=log(re2)
    if(re2.lt.20)then
        flow2x=6*pi*viscosity*drop2%r*v2*(1+exp(a10+a11*lg+a12*lg**2))/drop2%mass
    elseif(re2.lt.258)then
        flow2x=6*pi*viscosity*drop2%r*v2*(1+exp(a20+a21*lg))/drop2%mass
    else
        flow2x=6*pi*viscosity*drop2%r*v2*(1+exp(a30+a31*lg+a32*lg**2))/drop2%mass
    end if
    flow2y=flow2x*(drop2%vy-u1y)/v2
    flow2x=flow2x*(drop2%vx-u1x)/v2
elseif(flowcorrect.ne.3)then
    flow1x=g*(drop1%vx-u2x)/drop1%vt
    flow1y=g*(drop1%vy-u2y)/drop1%vt
    flow2x=g*(drop2%vx-u1x)/drop2%vt
    flow2y=g*(drop2%vy-u1y)/drop2%vt
else

    x=drop2%x-drop1%x
    y=drop2%y-drop1%y
    r=sqrt(x**2+y**2)
    cost=y/r
    sint=x/r
    alp1=drop1%r/r
    alp2=drop2%r/r
    u1=drop1%vy
    v1=drop1%vx
    u2=drop2%vy
    v2=drop2%vx
    re1=airdensity/viscosity*sqrt(u1**2+v1**2)*drop1%r*2
    re2=airdensity/viscosity*sqrt(u2**2+v2**2)*drop2%r*2
    re11=airdensity/viscosity*sqrt(u1**2+v1**2)*r
    re22=airdensity/viscosity*sqrt(u2**2+v2**2)*r
call klett1973(u1,u2,re1,re2,re11,re22,alp1,alp2,sint,cost,drop1%r,drop2%r,AA1,AA2,FF1,FF2)

    flow1y=6*pi*viscosity*drop1%r*u1*AA1/drop1%mass
    flow1x=6*pi*viscosity*drop1%r*u1*FF1/drop1%mass
    flow2y=6*pi*viscosity*drop2%r*u2*AA2/drop2%mass
    flow2x=6*pi*viscosity*drop2%r*u2*FF2/drop2%mass

call klett1973(v1,v2,re1,re2,re11,re22,alp1,alp2,cost,sint,drop1%r,drop2%r,AA1,AA2,FF1,FF2)
    flow1y=(flow1y-6*pi*viscosity*drop1%r*v1*FF1/drop1%mass)*drop1%c1
    flow1x=(-flow1x+6*pi*viscosity*drop1%r*v1*AA1/drop1%mass)*drop1%c1
    flow2y=(flow2y-6*pi*viscosity*drop2%r*v2*FF2/drop2%mass)*drop2%c1
    flow2x=(-flow2x+6*pi*viscosity*drop2%r*v2*AA2/drop2%mass)*drop2%c1
end if

    drop2%ax=fex/drop2%mass-flow2x
    drop2%ay=fey/drop2%mass+g-flow2y
    drop1%ax=(-fex)/drop1%mass-flow1x
    drop1%ay=(-e0*(drop1%q+drop2%q)-fey)/drop1%mass+g-flow1y

    if( mod(l,50)==0 )then
    x=drop2%x-drop1%x
    y=drop2%y-drop1%y
    r=sqrt(x**2+y**2)
    cost=y/r
    sint=x/r
    end if
    !call cout(drop1)
    !call cout(drop2)
    !pause
end subroutine

subroutine klett1973(u1,u2,re1,re2,re11,re22,alp1,alp2,sint,cost,r1,r2,AA1,AA2,FF1,FF2)

real::AA1,AA2,FF1,FF2,bb1,bb2,yy1,yy2,re11,re22
real::g1,g2,g3,g4,g5,g6,delta,o1,o2,alp1,alp2,qq,f01,f02,sint,cost,ep1,ep2,u1,u2,re1,re2,r1,r2
    ep1=exp(-0.5*re22*(1-cost))
    ep2=exp(-0.5*re11*(1+cost))
   ! ep2=ep1
    g1=(1-3*cost**2)/6.0*(ep1*alp1**2+alp2**2)-(1-ep1)*cost/re22+(1+cost)*ep1/2.0
    g2=(1-3*cost**2)/6.0*(ep2*alp2**2+alp1**2)+(1-ep2)*cost/re11+(1-cost)*ep2/2.0
    g3=sint*cost/2*(ep1*alp1**2+alp2**2)+(1-ep1)*sint/re22-sint*ep1/2.0
    g4=sint*cost/2*(ep2*alp2**2+alp1**2)-(1-ep2)*sint/re11+sint*ep2/2.0
    g5=(1-3*sint**2)/6.0*(ep1*alp1**2+alp2**2)+ep1*(1-cost)/2+(1-ep1)*(1+cost-cost**2)/re22/(1-cost)
    if(ep1.eq.1 .and. cost.eq.1)g5=(1-3*sint**2)/6.0*(ep1*alp1**2+alp2**2)+ep1*(1-cost)/2
    g6=(1-3*sint**2)/6.0*(ep2*alp2**2+alp1**2)+ep2*(1+cost)/2+(1-ep2)*(1-cost-cost**2)/re11/(1+cost)
    if(ep2.eq.1 .and. cost.eq.-1)g6=(1-3*sint**2)/6.0*(ep2*alp2**2+alp1**2)+ep2*(1+cost)/2
    qq=u2*r2/(u1*r1)
    if(u2.eq.0 .or. abs(u2/u1).lt.1e-10) qq=1e-10
    f01=1.5*qq*alp1
    f02=1.5/qq*alp2
    AA1=1-f01*g1-f01*f02*(g3*g4+g5*g6)+f01**2*f02*(g1*g5-g3**2)*g6
    AA2=1-f02*g2-f01*f02*(g3*g4+g5*g6)+f02**2*f01*(g2*g6-g4**2)*g5
    FF1= -f01*g3+f01*f02*(g2*g3+g4*g5)+f01**2*f02*(-g1*g5+g3**2)*g4
    FF2= -f02*g4+f01*f02*(g1*g4+g3*g6)+f02**2*f01*(-g2*g6+g4**2)*g3
    bb1=-(1.5+2.5*qq)+1.5*f01*g1*(1+qq)+0.5*f01*f02*qq*(3*g3*g4+2*g5*g6)
    bb2=-(1.5+2.5/qq)+1.5*f02*g2*(1+1.0/qq)+0.5*f01*f02/qq*(3*g3*g4+2*g5*g6)
    yy1=f01*g3*(1+1.5*qq)-0.5*f01*f02*qq*(3*g2*g3+2*g4*g5)
    yy2=f02*g4*(1+1.5/qq)-0.5*f01*f02/qq*(3*g1*g4+2*g3*g6)
    o1=(1+qq)*(2.5*(1-f01*f02*g3*g4)-0.5*f01*f02*(2*g5*g6+3*g1*g2) )
    o2=o1/qq
    delta=1-f01*f02*(g1*g2+2*g3*g4+g5*g6)+f01**2*f02**2*(g4**2-g2*g6)*(g3**2-g1*g5)
!print*,alp1,alp2,qq,f01,f02,sint,cost,ep1,ep2,u1,u2,re1,re2,r1,r2

    bb1=bb1/AA1+o1/delta
    bb2=bb2/AA2+o2/delta
    yy1=yy1/FF1+o1/delta
    yy2=yy2/FF2+o2/delta
    AA1=AA1/delta*(1+0.375*bb1*re1)
    AA2=AA2/delta*(1+0.375*bb2*re2)
    FF1=FF1/delta*(1+0.375*yy1*re1)
    FF2=FF2/delta*(1+0.375*yy2*re2)
  ! print*,AA1,AA2,FF1,FF2,bb1,bb2,yy1,yy2,re11,re22
end subroutine
subroutine moving(drop1,drop2) !步长移动 无修正
    use pub
    type(drop) drop1,drop2

    drop1%vx=drop1%vx+drop1%ax*dt
    drop1%vy=drop1%vy+drop1%ay*dt
    drop2%vx=drop2%vx+drop2%ax*dt
    drop2%vy=drop2%vy+drop2%ay*dt
    drop1%x=drop1%x+drop1%vx*dt-0.5*drop1%ax*dt*dt
    drop1%y=drop1%y+drop1%vy*dt-0.5*drop1%ay*dt*dt
    drop2%x=drop2%x+drop2%vx*dt-0.5*drop2%ax*dt*dt
    drop2%y=drop2%y+drop2%vy*dt-0.5*drop2%ay*dt*dt

end subroutine

subroutine cout(drop2) !输出液滴信息
    use pub
    type(drop) drop2
    write (*,*), 'drop2 x2=',drop2%x,drop2%y,"v2=",drop2%vx,drop2%vy,"a=",drop2%ax,drop2%ay
end subroutine

subroutine virtualmoving(drop1,drop2) !假设移动0.5dt
    use pub
    type(drop) drop1,drop2

    drop1%vx=drop1%vx+drop1%ax*dt/2.0
    drop1%vy=drop1%vy+drop1%ay*dt/2.0
    drop2%vx=drop2%vx+drop2%ax*dt/2.0
    drop2%vy=drop2%vy+drop2%ay*dt/2.0
    drop1%x=drop1%x+drop1%vx*dt/2.0-0.5*drop1%ax*dt*dt/4.0
    drop1%y=drop1%y+drop1%vy*dt/2.0-0.5*drop1%ay*dt*dt/4.0
    drop2%x=drop2%x+drop2%vx*dt/2.0-0.5*drop2%ax*dt*dt/4.0
    drop2%y=drop2%y+drop2%vy*dt/2.0-0.5*drop2%ay*dt*dt/4.0

end subroutine
subroutine caculatec2(drop1)
    use pub
    type(drop) drop1
    real x,y,z

    z=32*drop1%r**3*(drop1%density-airdensity)*airdensity*g/3.0/viscosity**2
    x=log(z)
    y=-3.18657+0.992696*x-0.00153193*x**2-(0.987059E-3)*x**3
    y=y-(0.578878E-3)*x**4+(0.855176E-4)*x**5-(0.327815E-5)*x**6
    drop1%renold=exp(y)
    drop1%c2=z/drop1%renold/24.0
    if(formatrix.eq.0)print*,"renolds=",drop1%renold
end subroutine

subroutine stokesflow(drop1,drop2,u1x,u1y,u2x,u2y,x,y,r)
    use pub

    type(drop) drop1,drop2
    real u1x,u1y,u2x,u2y,x,y,r,sinb,cosb,xx,yy,sin0,cos0,v,uxx,uyy,ratio,ur,utheta,re1,rn,prop
   ! print*,drop1%vy,drop2%vy
    if(flowcorrect.eq.3) return
if(relaflow.eq.0)then
  v=sqrt(drop1%vx**2+drop1%vy**2)!induced stokes flow 流场
  sinb=drop1%vx/v
  cosb=drop1%vy/v
else
  v=sqrt((drop1%vx-u2x)**2+( drop1%vy -u2y)**2)
  sinb=(drop1%vx-u2x)/v
  cosb=(drop1%vy-u2y)/v
end if
  xx=x*cosb-y*sinb
  yy=x*sinb+y*cosb
  sin0=xx/sqrt(xx**2+yy**2)
  cos0=yy/sqrt(xx**2+yy**2)
  ratio=drop1%r/r

  if(flowcorrect.eq.0)then
  uxx=v*sin0*cos0*0.75*(ratio-ratio**3)
  uyy=v*(0.75*ratio*(1+cos0**2)+0.25*ratio**3*(1-3*cos0**2))
  elseif(flowcorrect.eq.1)then
  re1=2*airdensity/viscosity*v*drop1%r
  rn=re1*3.0/16.0

  ur=-v*cos0*(1-ratio*3/2+ratio**3/2)*(1+rn)-v*rn*(cos0**2-sin0**2/2)*(1-ratio*3/2+ratio**2/2-ratio**3/2+ratio**4/2)
  utheta=v*sin0*((1-3/4*ratio-1/4*ratio**3)*(1+rn)+rn*cos0*(1-ratio*3/4+ratio**2/4-ratio**3/4-ratio**4/4)  )
 prop= max( min( (1/ratio-2)/3,1.0 ),0.0)
  uxx=(ur*sin0+utheta*cos0)*(1-prop)
  uyy=(ur*cos0-utheta*sin0)*(1-prop)+v
  ur=-v*cos0*(1+ratio**3/2)-3*v/re1*ratio**2*(((1-cos0)*re1/4/ratio+1)*exp(-(1+cos0)*re1/4/ratio)-1 )
  utheta=v*sin0*((1-ratio**3/4)-3*ratio/4*exp(-(1+cos0)*re1/4/ratio))
  uxx=uxx+(ur*sin0+utheta*cos0)*(prop)
  uyy=uyy+(ur*cos0-utheta*sin0)*(prop)

  elseif(flowcorrect.eq.2)then
  re1=2*airdensity/viscosity*v*drop1%r


if(re1.lt.800)then
  A1=-1.598e-18*re1**7 + 5.663e-15*re1**6 - 8.147e-12*re1**5 + &
6.05e-09*re1**4-2.403e-06*re1**3 + 0.0004472*re1**2 - 0.007474*re1 - 4.458
  B1=-2.27e-15*re1**6+6.26e-12*re1**5-6.619e-09*re1**4+3.228e-06*re1**3-0.0005917*re1**2-0.06869*re1-0.1604
elseif(re1.lt.2000)then
    A1= - 1.84e-07*re1**2 + 0.001061*re1 + 2.78
    B1=9.8e-7*re1**2-0.0057*re1-32
else
    A1= -2.888e-08*re1**2 + 0.0002851*re1 + 3.711
    B1=1.6e-07*re1**2-0.0016*re1-37
end if




  A2=-120.0/29.0-75*A1/29.0
  A3=(153+63*A1)/29.0
  A4=(-47.5-17*A1)/29.0
  B2=-69.0/27.0*B1
  B3=57.0/27.0*B1
  B4=-15.0/27.0*B1

  ur=-((ratio*A1+ratio**2*A2+ratio**3*A3+ratio**4*A4)*2*cos0 &
  +sig*(ratio*B1+ratio**2*B2+ratio**3*B3+ratio**4*B4)*(3*cos0**2-1) )*ratio**2*v
  utheta= ((-ratio**2*A1-2*ratio**3*A2-3*ratio**4*A3-4*ratio**5*A4)*sin0 &
  +sig*(-ratio**2*B1-2*ratio**3*B2-3*ratio**4*B3-4*ratio**5*B4)*sin0*cos0 )*ratio*v
  uxx=ur*sin0+utheta*cos0
  uyy=ur*cos0-utheta*sin0
  uxx=(uxx*re1+v*sin0*cos0*0.75*(ratio-ratio**3)/re1)/(re1+1/re1)
  uyy=(uyy*re1+( v*(0.75*ratio*(1+cos0**2)+0.25*ratio**3*(1-3*cos0**2)) )/re1)/(re1+1/re1)
  end if


  u1x=(uxx*cosb+uyy*sinb)*drop1%c1
  u1y=(-uxx*sinb+uyy*cosb)*drop1%c1



!print*,u1x,u1y
   !!!!!!!!!!!2
if(relaflow.eq.0)then
  v=sqrt(drop2%vx**2+drop2%vy**2)!induced stokes flow 流场
  sinb=drop2%vx/v
  cosb=drop2%vy/v
else
  v=sqrt((drop2%vx-u1x)**2+( drop2%vy -u1y)**2)
  sinb=(drop2%vx-u1x)/v
  cosb=(drop2%vy-u1y)/v
end if

  xx=x*cosb-y*sinb
  yy=x*sinb+y*cosb
  sin0=xx/sqrt(xx**2+yy**2)
  cos0=yy/sqrt(xx**2+yy**2)
  ratio=drop2%r/r

  if(flowcorrect.eq.0)then
  uxx=v*sin0*cos0*0.75*(ratio-ratio**3)
  uyy=v*(0.75*ratio*(1+cos0**2)+0.25*ratio**3*(1-3*cos0**2))
  elseif(flowcorrect.eq.1)then
  re1=2*airdensity/viscosity*v*drop2%r
  rn=re1*3.0/16.0


  ur=-v*cos0*(1-ratio*3/2+ratio**3/2)*(1+rn)-v*rn*(cos0**2-sin0**2/2)*(1-ratio*3/2+ratio**2/2-ratio**3/2+ratio**4/2)
  utheta=v*sin0*((1-3/4*ratio-1/4*ratio**3)*(1+rn)+rn*cos0*(1-ratio*3/4+ratio**2/4-ratio**3/4-ratio**4/4)  )
  prop= max( min( (1/ratio-2)/3,1.0 ),0.0)
  uxx=(ur*sin0+utheta*cos0)*(1-prop)
  uyy=(ur*cos0-utheta*sin0)*(1-prop)+v
  ur=-v*cos0*(1+ratio**3/2)-3*v/re1*ratio**2*(((1-cos0)*re1/4/ratio+1)*exp(-(1+cos0)*re1/4/ratio)-1 )
  utheta=v*sin0*((1-ratio**3/4)-3*ratio/4*exp(-(1+cos0)*re1/4/ratio))
  uxx=uxx+(ur*sin0+utheta*cos0)*(prop)
  uyy=uyy+(ur*cos0-utheta*sin0)*(prop)

  elseif(flowcorrect.eq.2)then
  re1=2*airdensity/viscosity*v*drop2%r

if(re1.lt.800)then
  A1=-1.598e-18*re1**7 + 5.663e-15*re1**6 - 8.147e-12*re1**5 + &
6.05e-09*re1**4-2.403e-06*re1**3 + 0.0004472*re1**2 - 0.007474*re1 - 4.458
  B1=-2.27e-15*re1**6+6.26e-12*re1**5-6.619e-09*re1**4+3.228e-06*re1**3-0.0005917*re1**2-0.06869*re1-0.1604
elseif(re1.lt.2000)then
    A1= - 1.84e-07*re1**2 + 0.001061*re1 + 2.78
    B1=9.8e-7*re1**2-0.0057*re1-32
else
    A1= -2.888e-08*re1**2 + 0.0002851*re1 + 3.711
    B1=1.6e-07*re1**2-0.0016*re1-37
end if

  A2=-120.0/29.0-75*A1/29.0
  A3=(153+63*A1)/29.0
  A4=(-47.5-17*A1)/29.0
  B2=-69.0/27.0*B1
  B3=57.0/27.0*B1
  B4=-15.0/27.0*B1
  ur=-((ratio*A1+ratio**2*A2+ratio**3*A3+ratio**4*A4)*2*cos0 &
  +sig*(ratio*B1+ratio**2*B2+ratio**3*B3+ratio**4*B4)*(3*cos0**2-1) )*ratio**2*v
  utheta= ((-ratio**2*A1-2*ratio**3*A2-3*ratio**4*A3-4*ratio**5*A4)*sin0 &
  +sig*(-ratio**2*B1-2*ratio**3*B2-3*ratio**4*B3-4*ratio**5*B4)*sin0*cos0 )*ratio*v
  uxx=ur*sin0+utheta*cos0
  uyy=ur*cos0-utheta*sin0
  uxx=(uxx*re1+v*sin0*cos0*0.75*(ratio-ratio**3)/re1)/(re1+1/re1)
  uyy=(uyy*re1+( v*(0.75*ratio*(1+cos0**2)+0.25*ratio**3*(1-3*cos0**2)) )/re1)/(re1+1/re1)
  end if
  u2x=(uxx*cosb+uyy*sinb)*drop2%c1
  u2y=(-uxx*sinb+uyy*cosb)*drop2%c1

end subroutine

subroutine teleport(drop1,drop2,x1,y1,x2,y2)
    use pub
    type(drop) drop1,drop2
    real x1,y1,x2,y2
    drop1%x=drop1%x+x1
    drop2%x=drop2%x+x2
    drop1%y=drop1%y+y1
    drop2%y=drop2%y+y2
end subroutine

function s_fun(m,xi,b,precise)
    real xi,b,sumfun, y,n,t,w,add
    integer m,k,i,precise
    k=0
    sumfun=0
    do
        if(exp(-(2.0*k+1)).lt. 0.01 )then
            exit
        end if
        k=k+1
    end do

    do i=0,k-1,1
        add=(2*i+1)**m*exp((2.0*i+1.0)*xi)/( exp( (2.0*i+1)*b )-1 )
        sumfun=sumfun+add
    end do
        n=2*k+1
    do i=0,precise,1
        y=xi-b*i-b
        t=exp(2*y)
        w=1.0/(1-t)
        if(m.eq.0)then
            add=w*exp(n*y)
        elseif(m.eq.1)then
            add=(2*t*w**2+n*w)*exp(n*y)
        elseif(m.eq.2)then
            add=(8*t**2*w**3+4*t*w**2*(n+1)+n**2*w)*exp(n*y)
        elseif(m.eq.3)then
            add=(48*t**3*w**4+24*t**2*w**3*(n+2)+2*t*w**2*(3*n**2+6*n+4)+n**3*w)*exp(n*y)
        end if
        sumfun=sumfun+add

    end do
   s_fun=sumfun
end function

function ss_fun(m,xi,b,precise)
    real xi,b,sumfun, y,n,t,w,add
    integer m,k,i,precise
    k=0
    sumfun=0
    do
        if(exp(-(2.0*k+1)).lt. 0.01 )then
            exit
        end if
        k=k+1
    end do

    do i=0,k-1,1
        add=(-1)**i*(2*i+1)**m*exp((2*i+1.0)*xi)/( exp( (2*i+1)*b )-1 )
        sumfun=sumfun+add
    end do
        n=2*k+1
    do i=0,precise,1
        y=xi-b*i-b
        t=exp(2*y)
        w=1.0/(1+t)
        if(m.eq.0)then
            add=w*exp(n*y)
        elseif(m.eq.1)then
            add=(-2*t*w**2+n*w)*exp(n*y)
        elseif(m.eq.2)then
            add=(8*t**2*w**3-4*t*w**2*(n+1)+n**2*w)*exp(n*y)
        elseif(m.eq.3)then
            add=(-48*t**3*w**4+24*t**2*w**3*(n+2)-2*t*w**2*(3*n**2+6*n+4)+n**3*w)*exp(n*y)
        end if
        add=add*(-1)**k
        sumfun=sumfun+add
    end do
   ss_fun=sumfun
end function

function t_fun(m,xi,b,precise)
    real xi,b,sumfun, y,n,t,w,add
    integer m,k,i,precise
    k=0
    sumfun=0
    do
        if(exp(-(2.0*k+1)).lt. 0.01 )then
            exit
        end if
        k=k+1
    end do

    do i=0,k-1,1
        sumfun=sumfun+(2*i+1)**m*exp((2*i+1.0)*xi)/( exp( (2*i+1)*b )-1 )**2
    end do
        n=2*k+1
    do i=0,precise,1
        y=xi-b*i-2*b
        t=exp(2*y)
        w=1.0/(1-t)
        if(m.eq.0)then
            add=w*exp(n*y)
        elseif(m.eq.1)then
            add=(2*t*w**2+n*w)*exp(n*y)
        elseif(m.eq.2)then
            add=(8*t**2*w**3+4*t*w**2*(n+1)+n**2*w)*exp(n*y)
        elseif(m.eq.3)then
            add=(48*t**3*w**4+24*t**2*w**3*(n+2)+2*t*w**2*(3*n**2+6*n+4)+n**3*w)*exp(n*y)
        end if
        add=add*(i+1)
        sumfun=sumfun+add
    end do
   t_fun=sumfun
end function

function u_fun(m,xi,b,precise)
    real xi,b,sumfun, y,n,t,w,add,multiply
    integer m,k,i,j,precise
    k=0
    sumfun=0
    do
        if(exp(-(2.0*k+1)).lt. 0.01 )then
            exit
        end if

        k=k+1
    end do

    do i=0,k-1,1
        add=(2*i+1.0)**m*exp((2*i+1.0)*xi)/((exp((2*i+1)*b)-1)*(exp((2*i+3)*b)-1))
        sumfun=sumfun+add
     !   print*,i," =",add
    end do
        n=2*k+1
    do i=0,precise,1
        y=xi-b*i-2*b
        t=exp(2*y)
        w=1.0/(1-t)
        if(m.eq.0)then
            add=w*exp(n*y)
        elseif(m.eq.1)then
            add=(2*t*w**2+n*w)*exp(n*y)
        elseif(m.eq.2)then
            add=(8*t**2*w**3+4*t*w**2*(n+1)+n**2*w)*exp(n*y)
        elseif(m.eq.3)then
            add=(48*t**3*w**4+24*t**2*w**3*(n+2)+2*t*w**2*(3*n**2+6*n+4)+n**3*w)*exp(n*y)
        end if
        multiply=0
        do j=1,i+1,1
            multiply=multiply+exp(-2*j*b)
        end do
        add=add*multiply
       ! print*,"add=",add
        sumfun=sumfun+add
    end do

   u_fun=sumfun
end function

subroutine electric_accurate(drop1,drop2,f,fex,fey)
    use pub
    type(drop) drop1,drop2
    real x,y,r,r1,r2,q1,q2,a,miu1,miu2,alpha,gama,d1,d2

    real p11,p12,p22,c11,c12,c22,det,tempora,v1,v2,qq1,qq2,cosphi,sinphi
    real::fper=0
    integer i,precise,j
    precise=1
    y=drop2%y-drop1%y
    x=drop2%x-drop1%x
    r=sqrt(x**2+y**2)
    cosphi=y/r
    sinphi=x/r

    r1=drop1%r
    r2=drop2%r
    q1=drop1%q
    q2=drop2%q
    d1=(r**2+r1**2-r2**2)/2/r
    d2=(r**2+r2**2-r1**2)/2/r
    a=sqrt(d1**2-r1**2)
    miu1=-log(r1/(d1+a))
    miu2=-log(r2/(d2+a))

    alpha=((d1+a)/r1)**2
    gama=( ((d2+a)/r2)**2 +1)*0.5

    k_low(1)=(a/r2)**2
    k_low(5:7)=1.0/ k_low(1)
    k_low(6)=(r2/a)**2


    do i=1,4,1
        f_low(i)=t_fun(i-1,2*miu1+miu2,miu1+miu2,precise)
        f_low(i+4)=t_fun(i-1,miu1+miu2,miu1+miu2,precise)
        f_low(i+8)=t_fun(i-1,miu2,miu1+miu2,precise)
        f_low(i+12)=u_fun(i-1,2*miu1+miu2,miu1+miu2,precise)
        f_low(i+16)=u_fun(i-1,miu1+miu2,miu1+miu2,precise)
        f_low(i+20)=u_fun(i-1,miu2,miu1+miu2,precise)
    end do

    c11=2*s_fun(0,miu2,miu1+miu2,precise)
    c22=2*s_fun(0,miu1,miu1+miu2,precise)
    c12=-2*s_fun(0,0.0,miu1+miu2,precise)
    det=c11*c22-c12**2
    p11=c22/det
    p22=c11/det
    p12=-c12/det

    k_up(5,2)=p12**2
    k_up(6,2)=2*p12*p22
    k_up(7,2)=p22**2

    k_up(5,6)=-2*(p11*p12)
    k_up(6,6)=-2*(p11*p22+p12**2)
    k_up(7,6)=-2*(p22*p12)

    k_up(5,10)=p11**2
    k_up(6,10)=2*p12*p11
    k_up(7,10)=p12**2

    k_up(5,13:14)=-gama*alpha*p12**2
    k_up(6,13:14)=-2*gama*alpha*p12*p22
    k_up(7,13:14)=-gama*alpha*p22**2

    k_up(5,17:18)=gama*p11*p12*(alpha+1)
    k_up(6,17:18)=gama*(p11*p22+p12**2)*(alpha+1)
    k_up(7,17:18)=gama*p22*p12*(alpha+1)

    k_up(5,21:22)=-gama*p11**2
    k_up(6,21:22)=-2*gama*p12*p11
    k_up(7,21:22)=-gama*p12**2

if(outfield.eq.1)then
    qq1=2*(s_fun(1,miu2,miu1+miu2,precise)+s_fun(1,0.0,miu1+miu2,precise))
    qq2=-2*(s_fun(1,miu1,miu1+miu2,precise)+s_fun(1,0.0,miu1+miu2,precise))
    v1=-(p11*qq1+p12*qq2)
    v2=-(p12*qq1+p22*qq2)

    k_low(2)=k_low(1)*0.5
    k_low(3)=1
    k_low(4)=1
    k_low(8)=(a/r2)**2*( ((d2+a)/r2)**2 -1)/8.0
    k_low(9:10)=( ((d2+a)/r2)**2 -1)/4.0

    k_up(1,2)=v2**2
    k_up(1,3)=-2*v2
    k_up(1,4)=1
    k_up(1,6)=-2*v1*v2
    k_up(1,7)=2*(v1-v2)
    k_up(1,8)=2
    k_up(1,10)=v1**2
    k_up(1,11)=2*v1
    k_up(1,12)=1
    k_up(1,13)=-gama*alpha*v2*(v2-2)
    k_up(1,14)=-gama*alpha*(v2**2-4*v2+2)
    k_up(1,15)=-gama*alpha*(-2*v2+3)
    k_up(1,16)=-gama*alpha
    k_up(1,17)=-gama*(2*v1*alpha-v1*v2*(alpha+1)-2*v2 )
    k_up(1,18)=-gama*(v1*(3*alpha+1)-v2*(alpha+3)+(2-v1*v2)*(alpha+1))
    k_up(1,19)=-gama*(alpha+1)*(v1-v2+3)
    k_up(1,20)=-gama*(alpha+1)
    k_up(1,21)=-gama*v1*(v1+2)
    k_up(1,22)=-gama*(v1**2+4*v1+2)
    k_up(1,23)=-gama*(2*v1+3)
    k_up(1,24)=-gama

    k_up(2,2)=-1
    k_up(2,4)=1
    k_up(2,6)=2
    k_up(2,8)=-2
    k_up(2,10)=-1
    k_up(2,12)=1
    k_up(2,13)=3*gama*alpha
    k_up(2,14)=gama*alpha
    k_up(2,15)=-3*gama*alpha
    k_up(2,16)=-gama*alpha
    k_up(2,17)=-3*gama*(alpha+1)
    k_up(2,18)=-gama*(alpha+1)
    k_up(2,19)=3*gama*(alpha+1)
    k_up(2,20)=gama*(alpha+1)
    k_up(2,21)=3*gama
    k_up(2,22)=gama
    k_up(2,23)=-3*gama
    k_up(2,24)=-gama

    k_up(3,2)=2*p12*v2
    k_up(3,3)=-2*p12
    k_up(3,6)=-2*(p12*v1+p11*v2)
    k_up(3,7)=2*(p11-p12)
    k_up(3,10)=2*p11*v1
    k_up(3,11)=2*p11
    k_up(3,13)=-2*gama*alpha*p12*(v2-1)
    k_up(3,14)=-2*gama*alpha*p12*(v2-2)
    k_up(3,15)=2*gama*alpha*p12
    k_up(3,17)=-gama*(2*p11*alpha-(p12*v1+p11*v2)*(alpha+1)-2*p12 )
    k_up(3,18)=-gama*(p11*(3*alpha+1)-(p12*v1+p11*v2)*(alpha+1)-p12*(3+alpha) )
    k_up(3,19)=-gama*(alpha+1)*(p11-p12)
    k_up(3,21)=-2*gama*p11*(v1+1)
    k_up(3,22)=-2*gama*p11*(v1+2)
    k_up(3,23)=-2*gama*p11

    k_up(4,2)=2*p22*v2
    k_up(4,3)=-2*p22
    k_up(4,6)=-2*(p22*v1+p12*v2)
    k_up(4,7)=2*(p12-p22)
    k_up(4,10)=2*p12*v1
    k_up(4,11)=2*p12
    k_up(4,13)=-2*gama*alpha*p22*(v2-1)
    k_up(4,14)=-2*gama*alpha*p22*(v2-2)
    k_up(4,15)=2*gama*alpha*p22
    k_up(4,17)=-gama*(2*p12*alpha-(p22*v1+p12*v2)*(alpha+1)-2*p22 )
    k_up(4,18)=-gama*(p12*(3*alpha+1)-(p22*v1+p12*v2)*(alpha+1)-p22*(alpha+3.0) )
    k_up(4,19)=-gama*(alpha+1)*(p12-p22)
    k_up(4,21)=-2*gama*p12*(v1+1)
    k_up(4,22)=-2*gama*p12*(v1+2)
    k_up(4,23)=-2*gama*p12

    k_up(8,13)=2*alpha*(2*v2-1.0)
    k_up(8,14)=4*alpha*(v2-1)
    k_up(8,15)=-2*alpha

    k_up(8,17)=-v1*(1+3*alpha)-v2*(3+alpha)-2*(1-alpha)
    k_up(8,18)=-4*(alpha*v1+v2)+2*(1-alpha)
    k_up(8,19)=(v1-v2+6.0)*(1-alpha)
    k_up(8,20)=2*(1-alpha)
    k_up(8,21)=2*(2*v1+1)
    k_up(8,22)=4*(v1+1)
    k_up(8,23)=2.0

    k_up(9,13)=4*alpha*p12
    k_up(9,14)=4*alpha*p12
    k_up(9,17)=-p11*(1+3*alpha)-p12*(3+alpha)
    k_up(9,18)=-4*(alpha*p11+p12)
    k_up(9,19)=(p11-p12)*(1-alpha)
    k_up(9,21)=4*p11
    k_up(9,22)=4*p11

    k_up(10,13)=4*alpha*p22
    k_up(10,14)=4*alpha*p22
    k_up(10,17)=-p12*(1+3*alpha)-p22*(3+alpha)
    k_up(10,18)=-4*(alpha*p12+p22)
    k_up(10,19)=(p12-p22)*(1-alpha)
    k_up(10,21)=4*p12
    k_up(10,22)=4*p12
end if

if(outfield.eq.0)then
    do i=5,7,1
    tempora=0
    do j=1,24,1
        tempora=tempora+k_up(i,j)*f_low(j)
    end do
    f_up(i)=k_low(i)*tempora+beta(i)
    end do
    f=-(f_up(5)*q1**2+f_up(6)*q1*q2+f_up(7)*q2**2)*k/r2**2
else
    do i=1,10,1
    tempora=0
    do j=1,24,1
        tempora=tempora+k_up(i,j)*f_low(j)
    !    if(i.eq.8) print*,k_low(i),k_up(i,j),f_low(j),k_up(i,j)*f_low(j)
    end do
    f_up(i)=k_low(i)*tempora+beta(i)
    end do

    f=r2**2*e0**2/k*(f_up(1)*cosphi**2+f_up(2)*sinphi**2)+e0*cosphi*(f_up(3)*q1+f_up(4)*q2) &
    +(f_up(5)*q1**2+f_up(6)*q1*q2+f_up(7)*q2**2)*k/r2**2+e0*q2*cosphi
    f=-f
    fper=r2**2*e0**2/k*f_up(8)*2*cosphi*sinphi+e0*sinphi*(f_up(9)*q1+f_up(10)*q2)+e0*q2*sinphi
end if
    fex=f*sinphi+fper*cosphi
    fey=f*cosphi-fper*sinphi
 end subroutine

subroutine swap(a,b)
    real a,b,c
    c=a
    a=b
    b=c
end subroutine

function near(a,b)

    if( a.gt.b*0.99 .and. a.lt.b*1.01)then
        near=1
    else
        near=0
    end if

end function
