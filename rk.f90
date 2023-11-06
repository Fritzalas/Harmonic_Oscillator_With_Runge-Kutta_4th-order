program rk_solve
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8),allocatable       :: T(:),X1(:),X2(:)  !We have the arrays of our problem. T is the array of time,X1 is the array for space and X2 is the array for velocity
  real(8)                   :: Ti,Tf             !Initial time and final time
  real(8)                   :: X10,X20           !Initial values for space and velocity
  integer                   :: Nt                !The total steps from ti to tf
  integer                   :: i,j               !Useful for do loops
  !---------------------------------------------------------------------------------
  !User Interface:
  print *, "#Enter the number of spaces in [ti,tf]: "
  read  *, Nt
  print *, "#Enter the initial time ti and the final time tf: "
  read  *, Ti,Tf
  print *, "#Enter the initial values for space and velocity: "
  read  *, X10,X20
  print *, "#We have Nt,Ti,Tf,X10,X20 equal to: ",Nt,Ti,Tf,X10,X20
  !----------------------------------------------------------------------------------
  !Array Allocation:
  ALLOCATE(X1(Nt))
  ALLOCATE(X2(Nt))
  ALLOCATE(T(Nt))
  !----------------------------------------------------------------------------------
  !Mistake in given data:
  if (Tf<Ti) stop 'The final time is less than the initial time! Try Again!!!'
  !----------------------------------------------------------------------------------
  !Calculations:
  call RK(T,X1,X2,Ti,Tf,X10,X20,Nt)   !Useful subroutine to make the calculations
  !----------------------------------------------------------------------------------
  !Output file:
  !Calculation of time,position and velocity:
  open(unit=11,file='rk.dat')
  do i=1,Nt
     write(11,*) T(i),X1(i),X2(i)
  end do
  !Calculation of total energy:
  open(unit=12,file='energy.dat')
  do i=1,Nt
     write(12,*) T(i),(0.5D0 * X2(i)*X2(i))+(0.5D0 * 10.0D0 * X1(i) * X1(i))
  end do
  !Calculation of the difference between real value and caluclated value:
  open(unit=13,file='difference.dat')
  do i=1,Nt
     write(13,*) T(i), ABS(X1(i)-cos(sqrt(10.0D0)*T(i)))
  end do
  close(13)
  close(12)
  close(11)
  !----------------------------------------------------------------------------------
end program rk_solve
!====================================================================================
real(8) function f1(t,x1,x2)
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)     :: t,x1,x2 !Given time,space and velocity
  !----------------------------------------------------------------------------------
  !Calculations:
  f1=x2     !Because dx/dt=V=x2 from our problem
end function f1
!====================================================================================
real(8) function f2(t,x1,x2)
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)  :: t,x1,x2 !Given time,space and velocity
  !----------------------------------------------------------------------------------
  !Calculations:
  f2=-10.0D0 * x1 !We have an ideal oscillator
end function f2
!====================================================================================
subroutine RK(T,X1,X2,Ti,Tf,X10,X20,Nt)
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  integer                :: Nt
  real(8),dimension(Nt)  :: T,X1,X2
  real(8)                :: Ti,Tf,X10,X20
  integer                :: i,j
  real(8)                :: dt         !The step of the time dt=(Tf-Ti)/(Nt-1)
  real(8)                :: TS,X1S,X2S !Local variables that hold the time,space and velocity for every step
  !----------------------------------------------------------------------------------
  !Initial values:
  dt=(Tf-Ti)/(Nt-1)
  T(1)=Ti
  X1(1)=X10
  X2(1)=X20
  TS=Ti
  X1S=X10
  X2S=X20
  !----------------------------------------------------------------------------------
  !Calculations:
  do i=2,Nt
     call RKSTEP(TS,X1S,X2S,dt) !We call a subroutine to calculate the next time, next place and next velocity
     T(i)=TS
     X1(i)=X1S
     X2(i)=X2S   !New values to our arrays
  end do
  !----------------------------------------------------------------------------------
end subroutine RK
!====================================================================================
subroutine RKSTEP(TS,X1S,X2S,dt)
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)   :: TS,X1S,X2S,dt
  real(8)   :: f1,f2           !We declare our functions that we are going to use
  real(8)   :: k11,k12,k13,k14
  real(8)   :: k21,k22,k23,k24
  real(8)   :: h,h2,h6         !We have that h is dt, h2=h/2,h6=h/6
  !----------------------------------------------------------------------------------
  !Initial values:
  h  =  dt         !h =dt, integration step                                                                                                                                                                                               
  h2 =  0.5D0 * h  !h2=h/2                                                                                                                                                                                                                
  h6 =h/(6.0D0)    !h6=h/6
  !----------------------------------------------------------------------------------
  !Calculations:
  k11=f1(TS,X1S,X2S)
  k21=f2(TS,X1S,X2S)
  k12=f1(TS+h2,X1S+h2*k11,X2S+h2*k21)
  k22=f2(TS+h2,X1S+h2*k11,X2S+h2*k21)
  k13=f1(TS+h2,X1S+h2*k12,X2S+h2*k22)
  k23=f2(TS+h2,X1S+h2*k12,X2S+h2*k22)
  k14=f1(TS+h ,X1S+h *k13,X2S+h *k23)
  k24=f2(TS+h ,X1S+h *k13,X2S+h *k23)
  TS  =TS+h
  X1S =X1S+h6*(k11+2.0D0*(k12+k13)+k14)
  X2S =X2S+h6*(k21+2.0D0*(k22+k23)+k24) !Next steps
end subroutine RKSTEP
!====================================================================================
