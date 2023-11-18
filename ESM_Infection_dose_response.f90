program inf
!!This program calculate the total number of infection using dose response models
implicit none
include 'mpif.h'
!-------------------------------------------------------------------------------------
!Basic Input Parameters:
integer,parameter::Npm=300000,Npm2=10000000,Kpm=12
!Inputs related to Basic Input Parameters:
real::start,DIST,finish,IYMAX,IXMAX,IYMIN,IXMIN  
real::alpha,  deltat, b1, p1
integer::r,NP,NP1,NP2,I,k,N1,Nstep,p,i1,i2,j1,kk,nn,m,ii,jj,k1,j,l,Ni,kp,N_infect,ID,temp,Nin
real,dimension(Npm,3)::R0,R1
real,dimension(Npm,Kpm)::pmx1,pmx2
real,dimension(Kpm)::infection_distance,betalo, betahi, pmx, q1, q2 !Weibull parameters
real,dimension(Npm,Kpm,Kpm)::bsum
real,dimension(Npm,Kpm)::Pn, Qn, bs, PPn, QQn
integer,dimension(Npm)::itr,kc,index 
character,dimension(1)::Char, c1(npm)
real,dimension(Npm,1)::xp,yp,zp,x0,y0,z0
integer irank, isize, ierror, tag, status(MPI_STATUS_SIZE)
CHARACTER*100 MNO, MESSAGE, header
real    :: s
!--------------------------------------------------------------------------------------
call cpu_time(start)
!--------------------------------------------------------------------------------------

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierror)      

write(MESSAGE, *)irank
MNO =adjustl(trim(MESSAGE))
write(*,*), 'node',irank, MNO
!******************************Input%%%%%%%%%%%%%%%%%%
!Nin=1    !The number of potential index people
temp=1    !Initial ID for loop
IYMAX=720
IXMAX=115
IYMIN=-5
IXMIN=-5
!!!Z-coordinate???


!***************************************%%%%%%%%%%%%
! read the input parameters file:
open(unit=19,file='results.txt_'//MNO)
open(unit=15,file='prob_params_'//MNO) !Contains parameters for dose-response models
read(15,*)header
read(15,*)kp
read(15,*) alpha
read(15,*)(infection_distance(i), i=1,kp)
do i=1, kp
infection_distance(i)=infection_distance(i)*39.37
enddo
read(15,*) (betalo(i), i=1,kp)
read(15,*) (betahi(i), i=1,kp)
read(15,*) deltat
!read(15,*) (q1(i), i=1,kp)
!read(15,*) (q2(i), i=1,kp)
!read(15,*) pmx
read(15,*)(pmx(i), i=1,kp)
close (15)


!---------------------------------------------------

    Np1=0
    kk=0
    open(unit=2,file='xmol.xyz') ! Contain configurations and pedestrian positions
    read(2,*)Np
    read(2,*)
    do i=1, Np
    read(2,*)c1(i),(R0(i,N1),N1=1,3)
    kc(i)=INT(iachar(c1(i)))
    if(kc(i).eq.67)then
       Np1=Np1+1
    endif 
    enddo
    close(2)
    Nin=Np1
    Np2=Np-Np1
    write(19,*) Np, Np1, Np2, '1'
    write(*,*)'NP', Np, Np1, Np2, N_infect
!-----------Initialize summation array PPn QQn

do k1= 1, kp
do i1=1, Np1
   PPn(i1,k1)=0
   QQn(i1,k1)=0
enddo
enddo

!------------------------------------------
!Determine the number of steps:
Nstep = 0
    open(unit=2,file='xmol.xyz')
     DO
        READ (2,*, END=10)
         Nstep = Nstep + 1
          END DO
           10 CLOSE (2)
           Nstep=Nstep /(NP+2)
!-----------------------------------------------
Do ID=1,Nin
!call random_seed()
!call random_number(s)
!temp=floor(Np1*s)+1
print*,temp
    do i=1,Np1
        if (i .EQ. temp) then
        itr(i)=2
        else
        itr(i)=1
        endif
    enddo
    kk=1   !!!?????????/
    N_infect=1
    index(kk)=temp


    N_infect=kk
    kk=0

!----------------------------------loop start---------------------------------------


 
do k1= 1, kp
do i1=1, Np1
   Pn(i1,k1)=0       
   Qn(i1,k1)=0
do j1=1,N_infect
   bsum(i1,j1,k1)=0
enddo   
enddo
enddo


 
  open(unit=2,file='xmol.xyz')
!------------------------------ inner loop start-------------------------------------



do i1=1,Nstep 
! print*,Nstep 
  read (2,*) 
  read (2,*) 

  do j1=1,Np1
    read(2,*) Char,xp(j1,1),yp(j1,1),zp(j1,1), ii !read passengers data
  end do

  do j1=1,Np2
    read(2,*) Char,x0(j1,1),y0(j1,1),z0(j1,1),ii !read obstacle data
  end do

  
     !! DATA FOR INFECTIVE r
   ! Estimating the total number of contacts of each passenger with the imported infective r

   
 do kk=1,N_infect
   nn=1

  do m=1,Np1
    if (m.NE.index(kk))then
       if ((yp(m,1).LE. IYMAX).and.(xp(m,1).LE. IXMAX)) then
       if ((yp(m,1).GE. IYMIN).and.(xp(m,1).GE. IXMIN)) then
         DIST = sqrt( ((xp(m,1)-xp(index(kk),1))**2)+((yp(m,1)-yp(index(kk),1))**2))  ! +((zp(m,1)-zp(index(kk),1))**2) )

         do k1=1,kp
         if (DIST .LE. infection_distance(k1)) then

          b1=(deltat/60)*(1-DIST/infection_distance(k1))**alpha
          
          bsum(m,kk,k1)=bsum(m,kk,k1)+b1
          end if
          enddo
       end if
       end if 
       nn=nn+1
    end if
  end do 
 end do

end do

!-------------------------------------------- end of inner loop --------------------------------------------

!-------------------------------------------------------
close(2,STATUS='KEEP')




do k1=1,kp
do m=1,Np1
do kk=1,N_infect
if (itr(m).eq.1)then
 bs(m,k1)=bs(m,k1)+bsum(m,kk,k1)
endif
enddo
enddo
enddo


!Pn and Qn calculation

!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Dose response function @@@@@@@@@@@@@@@@@@@@
do k1=1,kp
do i1=1,Np1
        if (i1.GT.temp) then                                                     
                Pn(i1,k1)= (1- exp(-1*betalo(k1)*bs(i1,k1)))                     !Exponential model
                Qn(i1,k1)= (1- exp(-1*betahi(k1)*bs(i1,k1)))                     !betalo(k1)*bs(i1,k1)=dose 
!                Pn(i1,k1)= (1- (1+betalo(k1)*bs(i1,k1)/bt1(k1))**(-1*af1(k1)))    !Beta model
!                Qn(i1,k1)= (1- (1+betahi(k1)*bs(i1,k1)/bt1(k1))**(-1*af1(k1)))    
!                Pn(i1,k1)= (1- exp((-1*q1(k1)*(betalo(k1)*bs(i1,k1))**(q2(k1))) ) )   !Weibull model
!                Qn(i1,k1)= (1- exp((-1*q1(k1)*(betahi(k1)*bs(i1,k1))**(q2(k1))) ) )

        else
                Pn(i1,k1)=0
                Qn(i1,k1)=0   !!!!Set index case as 1
        endif
enddo
enddo

do k1=1,kp
do i1=1,Np1
PPn(i1,k1)= PPn(i1,k1) + Pn(i1,k1)/Nin
QQn(i1,k1)= QQn(i1,k1) + Qn(i1,k1)/Nin    
enddo
enddo
temp=temp+1
enddo !!Enddo for ID loop

write(19,*)header
write(19,*) 'inf_rad'
write(19,'(12f12.5)') (infection_distance(k1), k1=1,kp)
write(19,*) 'betalo, betahi'
write(19,'(25f12.5)')(betalo(k1),betahi(k1), k1=1,kp)
write(19,*) 'Average Prob inf: lo, hi....'
do i=1,Np1
write(19,'(2i5,25f12.5)') i,itr(i),(PPn(i,k1), k1=1,kp) 
enddo
write (19,*)'***************************************'
do i=1,Np1
write(19,'(2i5,25f12.5)') i,itr(i),(QQn(i,k1), k1=1,kp)
enddo

!-------------------------------------------------------------------------
call cpu_time(finish)
!-----------------------------------------------------
!write(19,*) 
write(19,*)'Simulation Elapsed Time = ', finish-start, 'seconds'

call MPI_FINALIZE(ierror)
!----------------------------------------------------------------

end program inf
