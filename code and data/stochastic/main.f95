
program main
real ::t=0,maxtime,massbar,watercontent
real:: r0=4,q0=10,log_dr=1.0/6.0,log_dq=0.5,qmaxindex=2.0,relaxationtime=7200 !弛豫时间
integer ::blocks=37,show=0,fileoutput=1,screeoutput=0,choosekernel=1,blocks_q=0,blocksqs,testing=0
integer i,j,k,k1,k2,l1,l2,eof,p1,p2,q1,q2
integer:: panels=4;
character(len=2) field !请输入-4 0 2 4
character(len=22) filename,folder !请输入文件序号
character condition !1无电荷2正负电荷3电场200+正负电荷4电场400+正负电荷5电场400+正电荷6电场400+负电荷
character initial ! 1 半径15μm 演化20分钟 2 半径9μm 演化250分钟 1 半径6μm 演化120分钟
!integer:: electric_terminal=0 !0末速度不考虑电场 1末速度考虑电场
real , allocatable:: numbers(:,:),mass(:),totalmass(:,:),radius(:),charge(:,:),kernel(:,:),n0(:),qallocate(:,:)
real , allocatable:: test(:),efficiency(:,:),terminal(:,:),matrix_coal(:,:)
real ::pi=3.1415926,dn,dni,dnj,dnk,dt=1,summass,addmass,sumcharge,addcharge,totalcharge1,totalcharge2
real addm1,addm2,addq1,addq2,addn1,addn2,qmean1,qmean2
real , allocatable:: coef(:), qbin(:)
real:: allnumbers=0,allmass=0,sum_dn=0
real::sigma=2.9,meancharge=0.0
!real , allocatable:: gaussian(:)
real,dimension(574):: temporal
integer::point=1
temporal(:)=0.0
  if(choosekernel.eq.1)then
      open(unit=77,file="scales.txt",status="old")
    read(77,*),r0,blocks,log_dr,q0,blocks_q,log_dq,qmaxindex,maxtime,massbar,watercontent
    log_dr=log_dr/3
  end if
  blocksqs=2*blocks_q+1
allocate(numbers(blocks+1,-blocks_q-1:blocks_q+1),mass(blocks+1),charge(blocks+1,-blocks_q-1:blocks_q+1))
allocate (totalmass(blocks+1,-blocks_q-1:blocks_q+1),radius(blocks+1)   ,test(blocks+1) )
allocate(kernel(blocks*blocksqs,blocks*blocksqs) ,n0(blocks+1),qallocate(blocks+1,-blocks_q:blocks_q) )
allocate(efficiency(blocks*blocksqs,blocks*blocksqs) )
allocate(matrix_coal(blocks*blocksqs,blocks*blocksqs) )
allocate( terminal(blocks,blocksqs)  )
allocate( coef(blocksqs),qbin(blocksqs) )

qbin=(/-32.0,-16.0,-8.0,-4.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,4.0,8.0,16.0,32.0/)
7 continue
if(testing.eq.1)read(*,*),sigma
do i=1,15,1
    coef(i)=1.0*exp( -(  abs(qbin(i))/sigma  )**2/2 )
end do
coef(8)=1
    meancharge=0.0
do i=9,15,1
    meancharge=meancharge+coef(i)*0.25*2**abs(i-8.0)

end do
    meancharge=meancharge/(sum(coef(9:15) )+0.5 )
if(testing.eq.1) then
    print*,"meancharge=",meancharge
    print*,coef
    goto 7
end if
  !



  print*,"please set initial number 1-3, condition number 1-6"
  read(*,*),initial,condition
  if(initial.eq."1")then
    maxtime=1800
    massbar=14137
    folder="30"
  elseif(initial.eq."2")then
    maxtime=3600
    massbar=3053
    folder="60"
  elseif(initial.eq."3")then
    maxtime=7200
    massbar=1150!1150是6.5  998是6.2  905是6 1098是6.4
    folder="120"
  end if

  if(condition.eq."1")then
    field="0"
    filename="spectrum1.txt"
    coef(1:15)=0
    coef(8)=1.0

  elseif(condition.eq."2")then
    field="0"
    filename="spectrum2.txt"


  elseif(condition.eq."3")then
    field="2"
    filename="spectrum3.txt"


  elseif(condition.eq."4")then
    field="4"
    filename="spectrum4.txt"


  elseif(condition.eq."5")then
    field="4"
    filename="spectrum5.txt"
    coef(9:15)=0

  elseif(condition.eq."6")then
    field="4"
    filename="spectrum6.txt"
    coef(1:7)=0



!  elseif(condition.eq."7")then
!    field="4"
!    filename="spectrum5.txt"
!    coef=(/2,2,0,0,0/)
!
!  elseif(condition.eq."8")then
!    field="4"
!    filename="spectrum6.txt"
!    coef=(/0,0,0,2,2/)
!
!  elseif(condition.eq."9")then
!    field="0"
!    filename="spectrum5.txt"
!    coef=(/2,2,0,0,0/)
!
!  elseif(condition.eq."0")then
!    field="0"
!    filename="spectrum6.txt"
!    coef=(/0,0,0,2,2/)
  end if
    print*,meancharge
    print*,coef
  filename=trim(folder)//"\"//filename

print*,initial," ",condition," ",filename
  if(fileoutput.eq.1)then
  open(unit=99,file=filename,status="replace")
  end if
  if(choosekernel.eq.1)then
  !print*,"请输入场的大小 百伏/cm"
  !read(*,*),field
  open(unit=33,file="E"//trim(field)//"00\matrixkernel.csv",status="old")
  open(unit=34,file="E"//trim(field)//"00\matrix.csv",status="old")
  open(unit=35,file="coal_efficiency_beard.csv",status="old")
  open(unit=512,file="terminal.csv",status="old")!末速度
  kernel=0

do i=1,blocks*blocksqs,1
    read(34,*,iostat=eof),efficiency(i,:)
    read(35,*,iostat=eof),matrix_coal(i,:)
end do

do i=1,blocks*blocksqs,1
    do j=1,blocks*blocksqs,1
        read(33,*,iostat=eof),kernel(i,j)
        kernel(i,j)=kernel(i,j)*min( max(matrix_coal(i,j),0.3 ),1.0 )
    end do

end do


do i=1,blocks,1

        read(512,*,iostat=eof),terminal(i,:)


  end do


!10 continue
  kernel=kernel*1e-12
  end if

    radius(1)=r0
    mass(1)=4.0/3.0*pi*r0**3
do i=2,blocks+1,1
    radius(i)=radius(i-1)*(2**log_dr )
    mass(i)=4.0/3.0*pi*radius(i)**3

end do

    charge(:,0)=0
    charge(:,1)=q0*(radius)**qmaxindex
    charge(:,-1)=-charge(:,1)

if(blocks_q.gt.0)then
do i=2,blocks_q+1,1
    charge(:,i)=charge(:,i-1)*(2**log_dq )
    charge(:,-i)=-charge(:,i)
end do
end if

numbers=0
totalmass=0
qallocate=0
do i=1,blocks+1,1
    qallocate(i,:)=coef
end do

!qallocate(:,0)=5!-min(1,blocks_q):min(1,blocks_q)
do i=1,blocks+1,1
    n0(i)=exp(-mass(i)/massbar)*mass(i)*watercontent/massbar**2*0.51986!*log_dr*3

    allnumbers=allnumbers+n0(i)
    allmass=allmass+n0(i)*mass(i)
end do

if(blocks_q.eq.0)then
test=qallocate(:,0)
else
test=sum(qallocate(:,-blocks_q:blocks_q),dim=2)
end if





do i=1,blocks+1,1
    numbers(i,-blocks_q:blocks_q)=qallocate(i,:)*n0(i)/test(i)  !个/立方厘米!!!!!!!!
    totalmass(i,-blocks_q:blocks_q)=numbers(i,-blocks_q:blocks_q)*mass(i)

end do

!!!!!!!!!!!!!!!!!!!!!



if(show.eq.1)then
print*,"numbers ",numbers(:,-blocks_q:blocks_q)
print*,"radius ",radius
print*,"mass  ",mass
end if

  if(fileoutput.eq.1)then


    totalcharge1=0.0
    totalcharge2=0.0
    do i=1,blocks,1
        do j=1,blocks_q,1

            totalcharge1=totalcharge1+numbers(i,j)*charge(i,j)
            totalcharge2=totalcharge2+numbers(i,-j)*charge(i,-j)
        end do
    end do

    allmass=0
    allnumbers=0
    n0=numbers(:,0)!-blocks_q:blocks_q
    do i=1,blocks+1,1
        do j=1,blocks_q,1
            n0(i)=n0(i)+numbers(i,j)
            n0(i)=n0(i)+numbers(i,-j)
        end do
        allnumbers=allnumbers+n0(i)
        allmass=allmass+n0(i)*mass(i)
    end do

    write(99,*),totalmass(:,-blocks_q:blocks_q),allnumbers,allmass,totalcharge1,totalcharge2
    print *,"numbers and mass=",allnumbers,allmass,totalcharge1,totalcharge2
  end if

    if(screeoutput.eq.1)then
      write(*,*),totalmass(:,-blocks_q:blocks_q),allnumbers,allmass,totalcharge1,totalcharge2

    end if



do while(t.lt.maxtime)

t=t+dt
p1=2
p2=1
q1=-blocks_q
q2=-blocks_q
do i=(p1*blocksqs-blocks_q+q1 ),blocks*blocksqs,1
    do j=1,(p1-1)*blocksqs,1


            dn=dt*numbers(p1,q1)*numbers(p2,q2)*kernel(j,i)


        if(numbers(p1,q1).le.0)then
            numbers(p1,q1)=0
        end if
        if(numbers(p2,q2).le.0)then
            numbers(p2,q2)=0
        end if
        if(dn.le.0)then
            goto 21
        end if
        if(dn.gt.min(numbers(p1,q1),numbers(p2,q2)))then
            dn=min(numbers(p1,q1),numbers(p2,q2))
        end if
        numbers(p1,q1)=numbers(p1,q1)-dn
        numbers(p2,q2)=numbers(p2,q2)-dn
        totalmass(p1,q1)=numbers(p1,q1)*mass(p1)
        totalmass(p2,q2)=numbers(p2,q2)*mass(p2)
        summass=mass(p1)+mass(p2)
        sumcharge=charge(p1,q1)+charge(p2,q2)
        addmass=summass*dn
        addcharge=sumcharge*dn
        k1=0
        k2=1
        l1=0
        l2=1
        if(summass.le.mass(blocks+1))then
        do k=1,blocks,1
            k1=k1+1
            k2=k2+1
            if(summass.le.mass(k2) .and. summass.ge.mass(k1))then

                exit
            end if

        end do
        else
            k1=blocks
            k2=blocks+1
            summass=mass(k2)
        end if

        addm1=addmass*(mass(k2)-summass)/(mass(k2)-mass(k1))
        addm2=addmass*(summass-mass(k1))/(mass(k2)-mass(k1))
        addq1=addcharge*(mass(k2)-summass)/(mass(k2)-mass(k1))
        addq2=addcharge*(summass-mass(k1))/(mass(k2)-mass(k1))
        addn1=addm1/mass(k1)
        addn2=addm2/mass(k2)
        qmean1=addq1/addn1
        qmean2=addq2/addn2

        !!!!!!!!!!!!!!!1
        l1=0
        l2=1
        if(abs(qmean1).le.charge(k1,blocks_q))then
        do l=0,blocks_q,1

            if(abs(qmean1).le.charge(k1,l2) .and. abs(qmean1).ge.charge(k1,l1))then
                exit
            end if
            l1=l1+1
            l2=l2+1
        end do
        if(qmean1.lt.0)then
            l1=-l1
            l2=-l2
        end if
        else
            if(qmean1.gt.0)then
            l1=blocks_q
            l2=blocks_q+1
            qmean1=charge(k1,l1)
            else
            l1=-blocks_q
            l2=-blocks_q-1
            qmean1=charge(k1,l1)
            end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!


        totalmass(k1,l1)=totalmass(k1,l1)+addm1*(charge(k1,l2)-qmean1)/(charge(k1,l2)-charge(k1,l1))
        totalmass(k1,l2)=totalmass(k1,l2)+addm1*(qmean1-charge(k1,l1))/(charge(k1,l2)-charge(k1,l1))
        numbers(k1,l1)=totalmass(k1,l1)/mass(k1)
        numbers(k1,l2)=totalmass(k1,l2)/mass(k1)

       !!!!!!!!!!!!!!!1
        l1=0
        l2=1
        if(abs(qmean2).le.charge(k2,blocks_q))then
        do l=0,blocks_q,1

            if(abs(qmean2).le.charge(k2,l2) .and. abs(qmean2).ge.charge(k2,l1))then
                exit
            end if
            l1=l1+1
            l2=l2+1
        end do
        if(qmean2.lt.0)then
            l1=-l1
            l2=-l2
        end if
        else
            if(qmean2.gt.0)then
            l1=blocks_q
            l2=blocks_q+1
            qmean2=charge(k2,l1)
            else
            l1=-blocks_q
            l2=-blocks_q-1
            qmean2=charge(k2,l1)
            end if
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!

        totalmass(k2,l1)=totalmass(k2,l1)+addm2*(charge(k2,l2)-qmean2)/(charge(k2,l2)-charge(k2,l1))
        totalmass(k2,l2)=totalmass(k2,l2)+addm2*(qmean2-charge(k2,l1))/(charge(k2,l2)-charge(k2,l1))
        numbers(k2,l1)=totalmass(k2,l1)/mass(k2)
        numbers(k2,l2)=totalmass(k2,l2)/mass(k2)


21 continue
    if(q2.lt.blocks_q)then  !!!p1p2q1q2指针移动
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
    end do

end do



if (relaxationtime.ge.0)then !电荷流失
    do i=1,blocks,1
        do j=1,blocks_q,1
            numbers(i,0)=numbers(i,0)+numbers(i,j)/relaxationtime
            numbers(i,j)=numbers(i,j)-numbers(i,j)/relaxationtime
            numbers(i,0)=numbers(i,0)+numbers(i,-j)/relaxationtime
            numbers(i,-j)=numbers(i,-j)-numbers(i,-j)/relaxationtime

            totalmass(i,0)=totalmass(i,0)+totalmass(i,j)/relaxationtime
            totalmass(i,j)=totalmass(i,j)-totalmass(i,j)/relaxationtime
            totalmass(i,0)=totalmass(i,0)+totalmass(i,-j)/relaxationtime
            totalmass(i,-j)=totalmass(i,-j)-totalmass(i,-j)/relaxationtime
        end do
    end do
end if


if (mod(t,maxtime/60).eq.0)then
    totalcharge1=0.0
    do i=1,blocks,1
        do j=1,blocks_q,1
            totalcharge1=totalcharge1+numbers(i,j)*charge(i,j)
        end do
    end do

    allnumbers=0
    n0=numbers(:,0)!-blocks_q:blocks_q
    do i=1,blocks+1,1
        do j=1,blocks_q,1
            n0(i)=n0(i)+numbers(i,j)
            n0(i)=n0(i)+numbers(i,-j)
        end do
        allnumbers=allnumbers+n0(i)
    end do
    temporal(point)=t/(maxtime/60)
    temporal(point+1)=allnumbers
    temporal(point+2)=totalcharge1
    point=point+3
end if


if (mod(t,maxtime/panels).eq.0)then
    totalcharge1=0.0
    totalcharge2=0.0
    do i=1,blocks,1
        do j=1,blocks_q,1

            totalcharge1=totalcharge1+numbers(i,j)*charge(i,j)
            totalcharge2=totalcharge2+numbers(i,-j)*charge(i,-j)
        end do
    end do

    allmass=0
    allnumbers=0
    n0=numbers(:,0)!-blocks_q:blocks_q
    do i=1,blocks+1,1
        do j=1,blocks_q,1
            n0(i)=n0(i)+numbers(i,j)
            n0(i)=n0(i)+numbers(i,-j)
        end do
        allnumbers=allnumbers+n0(i)
        allmass=allmass+n0(i)*mass(i)
    end do



      write(99,*),totalmass(:,-blocks_q:blocks_q),allnumbers,allmass,totalcharge1,totalcharge2
    print *,"numbers and mass=",allnumbers,allmass,totalcharge1,totalcharge2
end if

end do
 write(99,*),temporal
!print*,"result"
  if(fileoutput.eq.1)then
  !write(99,*),totalmass(:,-blocks_q:blocks_q)
  close(99)
  close(33)
  close(77)
  end if






deallocate(numbers,mass,totalmass,radius,charge,kernel   )

end

subroutine findchargeplace(sumcharge,l1,l2,blocks_q,charge)
    integer l1,l2,blocks_q
    real summcharge,charge
    l1=0
    l2=1

end subroutine
