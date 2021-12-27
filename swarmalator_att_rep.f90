     !!! program for swarmalators with attractive-repulsive phase coupling
      program swarmalators
      implicit none
      integer, parameter :: N = 100  !!! Number of swarmalators
      integer :: i, ISEED, IRAND, iteration, i1, j1, m,               &
      transient, m1
      integer, dimension(N) :: ind, ind1
      integer, dimension(N,N) :: LA, LR 
      double precision, dimension(N):: theta, phi, thetaa, dtheta, x, &
      y, z, dx, dy, xx, yy, k11, k12, k13, k21, k22, k23, k31, k32,   &
      k33, k41, k42, k43, k51, k52, k53, k61, k62, k63, dtheta1,      &
      dtheta2
      double precision, dimension(N,N) :: d1
      double precision :: g, dt, t0, J, A, B, KA, KR, pi, r


      ISEED = time()
      CALL SRAND (ISEED)

      open(1, file = 'att_rep_position.dat')
      open(11, file = 'att_rep_theta_time_series.dat')


      iteration = 50001   !!! number of iterations
      transient = iteration - 5000

      dt = 0.01d0   !!! integration step-size
      t0 = 0.d0
      A = 1.d0
      B = 1.d0
      
      r = 0.00001d0  !!! vision radius

      J = 0.1d0    !!! phase dependent spatial coupling strength
      KA = 0.5d0   !!! attractive phase coupling strength
      KR = -0.5d0  !!! repulsive phase coupling strength
      
      pi = 4.d0*datan(1.d0)      

!!!   initial conditions

      do i = 1,N
      x(i) = -1.d0 + 2.d0 * rand()
      y(i) = -1.d0 + 2.d0 * rand()
      theta(i) = -pi + 2.d0*pi*rand()
      end do
      

      iloop: do i = 1, iteration
      write(*,*)i
      
      
!!!   save the values of x, y, theta after each 100 iterations      
      if(mod(i-1,100) .eq. 0)then
      do m =1,N
      write(1,*)x(m)
      write(1,*)y(m)
      write(1,*)theta(m)
      end do
      end if 
      
      
!!!   save the values of theta after transients      
      if(i .gt. transient)then
      do m =1,N
      write(11,*)theta(m)
      end do  
      end if            


!!!   at the start of each integration step, initialize both the attractive and repulsive adjacency matrices to null matrix	
      do i1 = 1, N
      do j1 = 1, N
      LA(i1,j1) = 0     !!! adjacency matrix for attractive phase coupling
      LR(i1,j1) = 0     !!! adjacency matrix for repulsive phase coupling
      d1(i1,j1) = 0.d0  !!! spatial distance between the i1-th and j1-th swarmalator
      end do
      ind(i1) = 0      !!! number of swarmalators inside the vision range of the i1-th swarmalator
      ind1(i1) = 0     !!! number of swarmalators outside the vision range of the i1-th swarmalator
      end do 
      
      do i1 = 1, N
      do j1 = i1, N
      d1(i1,j1) = dist(x(i1), y(i1), x(j1), y(j1))
      d1(j1,i1) = d1(i1,j1)
      end do
      end do 
      
      do i1 = 1,N
      do j1 = 1,N
      if (d1(i1,j1) .le. r)then
      LA(i1,j1) = 1
      ind(i1) = ind(i1) + 1
      end if
      if (d1(i1,j1) .gt. r)then
      LR(i1,j1) = 1
      ind1(i1) = ind1(i1) + 1
      end if
      end do
      LA(i1,i1) = 0
      end do
      
      do i1 = 1,N
      if (ind(i1) .eq. 1)then  !!! this is done to avoid divison by zero
      ind(i1) = 2
      end if
      if (ind1(i1) .eq. 0)then  !!! this is done to avoid divison by zero
      ind1(i1) = 1  
      end if 
      end do                
	
!!!  integration by RK-4 method


      do i1 =1,N
      dx(i1) = 0.d0
      dy(i1) = 0.d0
      dtheta(i1) = 0.d0
      dtheta1(i1) = 0.d0
      dtheta2(i1) = 0.d0
      
      do 10 j1 = 1,N
      if (j1 .eq. i1) go to 10
      d1(i1,j1) = dist(x(i1),y(i1),x(j1),y(j1))
      dx(i1) = dx(i1) + ((x(j1)-x(i1)) * ( A + J*dcos(theta(j1)-theta(i1))) &
        / d1(i1,j1)) - (B * (x(j1)-x(i1)) / (d1(i1,j1)**2))
      dy(i1) = dy(i1) + ((y(j1)-y(i1)) * ( A + J*dcos(theta(j1)-theta(i1))) &
        / d1(i1,j1)) - (B * (y(j1)-y(i1)) / (d1(i1,j1)**2))
      dtheta1(i1) = dtheta1(i1) + LA(i1,j1) * KA * dsin(theta(j1)-theta(i1)) / d1(i1,j1)
      dtheta2(i1) = dtheta2(i1) + LR(i1,j1) * KR * dsin(theta(j1)-theta(i1)) / d1(i1,j1)
10    continue

      dx(i1) = dx(i1) / N
      dy(i1) = dy(i1) / N
      dtheta(i1) = dtheta(i1) + (dtheta1(i1) / (ind(i1)-1)) + (dtheta2(i1) / ind1(i1))

      k11(i1) = dx(i1) * dt
      k12(i1) = dy(i1) * dt
      k13(i1) = dtheta(i1) * dt
      end do

      do i1 =1,N
       xx(i1) = x(i1) + k11(i1) / 2.0
       yy(i1) = y(i1) + k12(i1) / 2.0
       thetaa(i1) = theta(i1) + k13(i1) / 2.0      
      end do
      
      do i1 =1,N
      dx(i1) = 0.d0
      dy(i1) = 0.d0
      dtheta(i1) = 0.d0
      dtheta1(i1) = 0.d0
      dtheta2(i1) = 0.d0      
      
      do 20 j1 = 1,N
      if (j1 .eq. i1) go to 20
      d1(i1,j1) = dist(xx(i1),yy(i1),xx(j1),yy(j1))
      dx(i1) = dx(i1) + ((xx(j1)-xx(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (xx(j1)-xx(i1)) / d1(i1,j1)**2)
      dy(i1) = dy(i1) + ((yy(j1)-yy(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (yy(j1)-yy(i1)) / d1(i1,j1)**2)       
      dtheta1(i1) = dtheta1(i1) + LA(i1,j1) * KA * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
      dtheta2(i1) = dtheta2(i1) + LR(i1,j1) * KR * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
20    continue

      dx(i1) = dx(i1) / N
      dy(i1) = dy(i1) / N
      dtheta(i1) = dtheta(i1) + (dtheta1(i1) / (ind(i1)-1)) + (dtheta2(i1) / ind1(i1))

      k21(i1) = dx(i1) * dt
      k22(i1) = dy(i1) * dt
      k23(i1) = dtheta(i1) * dt
      end do
      
      do i1 =1,N
       xx(i1) = x(i1) + k21(i1) / 2.0
       yy(i1) = y(i1) + k22(i1) / 2.0
       thetaa(i1) = theta(i1) + k23(i1) / 2.0      
      end do
      
      do i1 =1,N
      dx(i1) = 0.d0
      dy(i1) = 0.d0
      dtheta(i1) = 0.d0
      dtheta1(i1) = 0.d0
      dtheta2(i1) = 0.d0      
      
      do 30 j1 = 1,N
      if (j1 .eq. i1) go to 30
      d1(i1,j1) = dist(xx(i1),yy(i1),xx(j1),yy(j1))
      dx(i1) = dx(i1) + ((xx(j1)-xx(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (xx(j1)-xx(i1)) / d1(i1,j1)**2)
      dy(i1) = dy(i1) + ((yy(j1)-yy(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (yy(j1)-yy(i1)) / d1(i1,j1)**2)       
      dtheta1(i1) = dtheta1(i1) + LA(i1,j1) * KA * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
      dtheta2(i1) = dtheta2(i1) + LR(i1,j1) * KR * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
30    continue

      dx(i1) = dx(i1) / N
      dy(i1) = dy(i1) / N
      dtheta(i1) = dtheta(i1) + (dtheta1(i1) / (ind(i1)-1)) + (dtheta2(i1) / ind1(i1))

      k31(i1) = dx(i1) * dt
      k32(i1) = dy(i1) * dt
      k33(i1) = dtheta(i1) * dt
      end do
      

      do i1 =1,N
      xx(i1) = x(i1) + k31(i1)
      yy(i1) = y(i1) + k32(i1)
      thetaa(i1) = theta(i1) + k33(i1)      
      end do

      do i1 =1,N
      dx(i1) = 0.d0
      dy(i1) = 0.d0
      dtheta(i1) = 0.d0
      dtheta1(i1) = 0.d0
      dtheta2(i1) = 0.d0      
      
      do 40 j1 = 1,N
      if (j1 .eq. i1) go to 40
      d1(i1,j1) = dist(xx(i1),yy(i1),xx(j1),yy(j1))
      dx(i1) = dx(i1) + ((xx(j1)-xx(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (xx(j1)-xx(i1)) / d1(i1,j1)**2)
      dy(i1) = dy(i1) + ((yy(j1)-yy(i1)) * ( A + J*dcos(thetaa(j1)-         &
      thetaa(i1))) / d1(i1,j1)) - (B * (yy(j1)-yy(i1)) / d1(i1,j1)**2)       
      dtheta1(i1) = dtheta1(i1) + LA(i1,j1) * KA * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
      dtheta2(i1) = dtheta2(i1) + LR(i1,j1) * KR * dsin(thetaa(j1)-thetaa(i1)) / d1(i1,j1)
40    continue

      dx(i1) = dx(i1) / N
      dy(i1) = dy(i1) / N
      dtheta(i1) = dtheta(i1) + (dtheta1(i1) / (ind(i1)-1)) + (dtheta2(i1) / ind1(i1))

      k41(i1) = dx(i1) * dt
      k42(i1) = dy(i1) * dt
      k43(i1) = dtheta(i1) * dt
      end do
      
      do i1 =1,N
      x(i1) = x(i1) + (k11(i1) + 2.0 * (k21(i1) + k31(i1))                  &
      + k41(i1)) / 6.0
      y(i1) = y(i1) + (k12(i1) + 2.0 * (k22(i1) + k32(i1))                  &
      + k42(i1)) / 6.0
      theta(i1) = theta(i1) + (k13(i1) + 2.0 * (k23(i1) + k33(i1))          &
      + k43(i1)) / 6.0       
      end do
         

      t0 = t0 + dt


      end do iloop
      
!      do i = 1,N
!      phi(i) = datan(y(i)/x(i))
!      write(2,*)theta(i),phi(i)
!      end do
	

      stop
      
	contains


      function dist(a, b, c, d)
      implicit none
      double precision, intent(in) :: a, b, c, d
      double precision :: dist
      dist = dsqrt((a-c)**2.0 + (b-d)**2)
      end function dist
      end program swarmalators
