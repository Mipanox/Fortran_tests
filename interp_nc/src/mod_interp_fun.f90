!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for interpolation functions; partly inherited from NR  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_interp_fun

  use nrtype

  implicit none

  private

  interface spline
    module procedure spline
  end interface

  interface splint
    module procedure splint
  end interface

  interface splie2
    module procedure splie2
  end interface

  interface splin2
    module procedure splin2
  end interface

  interface splint_2d
    module procedure splint_2d
  end interface

  interface splint_3d
    module procedure splint_3d
  end interface

  interface interp_grid
    module procedure interp_grid1d
    module procedure interp_grid2d
    module procedure interp_grid3d
  end interface

  interface linint
    module procedure lin_1d
    module procedure bilin_2d
  end interface

  public :: locate
  public :: spline, splint
  public :: splie2, splin2, splint_2d
  public :: splint_3d

  public :: interp_grid

  public :: linint

contains

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  !-----      Auxiliary routines
  !--------------------------------------------------------------------- 
  function locate(xx,x)

    implicit none

    real(p_double), dimension(:), intent(in) :: xx
    real(p_double), intent(in) :: x
    integer :: locate
    ! Given an array xx(1:N), and given a value x, returns a value 
    ! j such that x is between xx(j) and xx(j+1). xx must be monotonic,
    ! either increasing or decreasing. j = 0 or j = N is returned 
    ! to indicate that x is out of range.

    integer :: n,jl,jm,ju
    logical :: ascnd

    n = size(xx)
    ascnd = (xx(n) >= xx(1)) ! True if ascending order of table,
                             ! false otherwise.
    jl = 0                   ! Initialize lower and upeer limits.
    ju = n+1

    do 
      if (ju-jl <= 1) exit   ! Repeat until this is satisfied
        
      jm = (ju+jl)/2         ! Compute a midpoint
      if (ascnd .eqv. (x >= xx(jm))) then
        jl = jm              ! and replace either the lower limit
      else
        ju = jm              ! or the upper limit, as appropriate.
      end if
    end do

    if (x == xx(1)) then     ! Then set the output, being careful with 
      locate = 1             ! the endpoints.
    else if (x == xx(n)) then
      locate = n-1
    else 
      locate = jl
    end if
  end function locate

  !---------------------------------------------------------------------
  subroutine tridag_ser(a,b,c,r,u)

    implicit none

    real(p_double), dimension(:), intent(in) :: a,b,c,r
    real(p_double), dimension(:), intent(out) :: u
    ! Solves for a vector u of size N the tridagonal linear set given by 
    ! equation (2.4.1) using a serial algorithm. Input vectors b (diagonal
    ! elements) and r (right-hand sides) have size N, while a and c 
    ! (off-diagonaml elements) are size N-1.

    real(p_double), dimension(size(b)) :: gam ! One vector of workspace, game is needed
    integer :: n,j
    real(p_double) :: bet

    n = size(b)

    bet = b(1)
    if (bet == 0.0) then
      print *, 'tridag_ser: Error at code stage 1'
      !call abort_program()
    end if

    u(1) = r(1)/bet
    do j=2,n                            ! Decomposition and forward substitution
      gam(j) = c(j-1) / bet
      bet = b(j) - a(j-1) * gam(j)

      if (bet == 0.0) then              ! Algorithm fails; see routine in Vol. 1
        print *, 'tridag_ser: Error at code stage 2'
        !call abort_program()
      end if

      u(j) = (r(j)-a(j-1)*u(j-1))/bet
    end do

    do j=n-1,1,-1                       ! Backsubstitution.
      u(j) = u(j) - gam(j+1) * u(j+1)
    end do
  end subroutine tridag_ser

  !---------------------------------------------------------------------
  recursive subroutine tridag_par(a,b,c,r,u)

    implicit none

    real(p_double), dimension(:), intent(in) :: a,b,c,r
    real(p_double), dimension(:), intent(out) :: u
    ! Solves for a vector u of size N the tridagonal linear set given by 
    ! equation (2.4.1) using a parallel algorithm. Input vectors b 
    ! (diagonal elements) and r (right-hand sides) have size N, while a 
    ! and c (off-diagonaml elements) are size N-1.

    integer, parameter :: NPAR_TRIDAG = 4 ! Determines when serial 
    integer :: n,n2,nm,nx,i               ! algorithm is invoked.
    real(p_double), dimension(size(b)/2) :: y,q,piva
    real(p_double), dimension(size(b)/2-1) :: x,z
    real(p_double), dimension(size(a)/2) :: pivc

    n = size(b)

    if (n < NPAR_TRIDAG) then
      call tridag_ser(a,b,c,r,u)
    else
      if (maxval(abs(b(1:n))) == 0.0) then     ! Algorithm fails
        print *, 'tridag_par: possible singular matrix'
        !call abort_program()
      end if
      
      n2 = size(y)
      nm = size(pivc)
      nx = size(x)

      piva = a(1:n-1:2) / b(1:n-1:2)           ! Zero the odd a's and even
      pivc = c(2:n-1:2) / b(3:n:2)             ! c's, giving x, y, z, q.

      y(1:nm) = b(2:n-1:2) - piva(1:nm)*c(1:n-2:2) - pivc*a(2:n-1:2)
      q(1:nm) = r(2:n-1:2) - piva(1:nm)*r(1:n-2:2) - pivc*r(3:n:2)

      if (nm < n2) then
        y(n2) = b(n) - piva(n2)*c(n-1)
        q(n2) = r(n) - piva(n2)*r(n-1)
      end if

      x = -piva(2:n2) * a(2:n-2:2)
      z = -pivc(1:nx) * c(3:n-1:2)

      call tridag_par(x,y,z,q,u(2:n:2))        ! Recurse and get even u's.
      
      u(1) = (r(1)-c(1)*u(2))/b(1)             ! Substitute and get odd u's.
      u(3:n-1:2) =  ( r(3:n-1:2)-a(2:n-2:2)*u(2:n-1:2) &
                    - c(3:n-1:2)*u(4:n:2) ) / b(3:n-1:2)
      
      if (nm == n2) u(n) = (r(n)-a(n-1)*u(n-1))/b(n)
      
      !!! NaN handling !!!
      !do i=1,n
      !  if (i==1) then
      !    if (isnan(u(1))) u(1)=0.0_p_double
      !  else
      !    if (isnan(u(i))) u(i)=u(i-1)
      !  end if
      !end do

    end if
  end subroutine tridag_par



  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  !-----      Cubic spline implementations
  !---------------------------------------------------------------------!
  !----- 1D ------------------------------------------------------------!
  !---------------------------------------------------------------------!
  
  subroutine spline(xx,y,yp1,ypn,y2)
    
    implicit none

    real(p_double), dimension(:), intent(in) :: xx,y
    real(p_double), intent(in) :: yp1, ypn
    real(p_double), dimension(:), intent(out) :: y2
    ! Given arrays x and y of length N containing a tabulated function,
    ! i.e., y_i = f(x_i), with x_1 < x_2 < ... < x_N, and given values 
    ! yp1 and ypn for the first derivative of the interpolating function
    ! at points 1 and N, respectively, this routine returns an array y2
    ! of length N that contains the second derivatives of the 
    ! interpolating function at the tabulated points x_i.
    ! If yp1 and/or ypn are equal to 1 x 10^30 or larger, the routine is
    ! signaled to set the corresponding boundary condition for a natural 
    ! spline, with zero second derivative on that boundary.

    integer :: n 
    real(p_double), dimension(size(xx)) :: a,b,c,r
    real(p_double), dimension(size(xx)) :: x
    
    n = size(xx)

    ! Smooth out oscillations
    x = xx*1000.0_p_double*size(xx)

    c(1:n-1) = x(2:n) - x(1:n-1)       ! Set up the tridiagonal equations
    r(1:n-1) = 6.0_p_double*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1) = r(2:n-1) - r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.0_p_double*(c(2:n-1)+a(2:n-1))

    b(1) = 1.0_p_double
    b(n) = 1.0_p_double

    if (yp1 > 0.99e30_p_double) then         ! The lower boundary condition is
      r(1) = 0.0_p_double                    ! set either to be "natural"
      c(1) = 0.0_p_double
    else                               ! or else to have a specified first 
      r(1) = (3.0_p_double / (x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) ! derivative
      c(1) = 0.5_p_double
    end if

    if (ypn > 0.99e30_p_double) then         ! Similarly for upper boundary
      r(n) = 0.0_p_double
      a(n) = 0.0_p_double
    else
      r(n) = (-3.0_p_double / (x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-yp1)
      a(n) = 0.5_p_double
    end if

    call tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))

  end subroutine spline

  !---------------------------------------------------------------------
  function splint(xa,ya,y2a,x)
    
    implicit none

    real(p_double), dimension(:), intent(in) :: xa,ya,y2a
    real(p_double), intent(in) :: x
    real(p_double) :: splint
    ! Given the arrays xa and ya, which tabulate a function (with the xa_i's)
    ! in increasing or decreasing order), and given the array y2a, which is 
    ! the output from spline above, and given a value of x, this routine 
    ! returns a cubic-spline interpolated value. The arrays xa, ya, and y2a
    ! are all of the same size.

    integer :: khi,klo,n
    real(p_double) :: a,b,h

    klo = max(min(locate(xa,x),n-1),1)
    ! We will find the right place in the table by means of locate's bisection
    ! algorithm. This is optimal if sequential calls to this routine are at 
    ! random values of x. If sequential calls are in order, and closely spaced, 
    ! one would do better to store previous values of klo and khi and test if 
    ! they remain appropriate on the next call.
    khi = klo + 1                      ! klo and khi now bracket the input x
    h = xa(khi) - xa(klo)

    if (h == 0.0) then
      print *,'bad xa input in splint' ! The xa's must be distinct
      !call abort_program()
    end if

    a = (xa(khi)-x) / h                ! Cubic spline polynomial is now evaluated.
    b = (x-xa(klo)) / h

    splint = a*ya(klo) + b*ya(khi) + &
             ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * (h**2) / 6.0_p_double

  end function splint



  !---------------------------------------------------------------------!
  !----- 2D ------------------------------------------------------------!
  !---------------------------------------------------------------------!
  subroutine splie2(x1a,x2a,ya,y2a)

    implicit none

    real(p_double), dimension(:), intent(in) :: x1a,x2a
    real(p_double), dimension(:,:), intent(in) :: ya
    real(p_double), dimension(:,:), intent(out) :: y2a
    ! Given an M x N tabulated function ya, and N tabulated independent
    ! variables x2a, this routine constructs one-dimensional natural Cubic
    ! splines of the rows of ya and returns the second derivatives in the 
    ! M x N array y2a. (The array x1a is included in the argument list
    ! merely for consestency with routine splin2.)

    integer :: j,m,ndum

    m    = size(x1a)
    ndum = size(x2a)

    do j=1,m
      call spline(x2a,ya(j,:),1.0e30_p_double,1.0e30_p_double,y2a(j,:))
      ! Values 1 x 10^30 signal a natural spline.
    end do

  end subroutine splie2

  !---------------------------------------------------------------------
  function splin2(x1a,x2a,ya,y2a,x1,x2)

    implicit none

    real(p_double), dimension(:), intent(in) :: x1a,x2a
    real(p_double), dimension(:,:), intent(in) :: ya,y2a
    real(p_double), intent(in) :: x1,x2
    real(p_double) :: splin2
    ! Given x1a, x2a, ya as described in splie2 and y2a as produced by 
    ! that routine; and given a desired interpolating point x1,x2;
    ! this routine returns an interpolated function value by bicubic 
    ! spline interpolation.
    integer :: j,m,ndum
    real(p_double), dimension(size(x1a)) :: yytmp,y2tmp2

    m    = size(x1a)
    ndum = size(x2a)

    do j=1,m
      yytmp(j) = splint(x2a,ya(j,:),y2a(j,:),x2)
      ! Perform m evaluations of the row splines constructed by splie2,
      ! using the one-dimensional spline evaluator splint.
    end do

    call spline(x1a,yytmp,1.0e30_p_double,1.0e30_p_double,y2tmp2)
    ! Construct the one-dimensional column spline and evaluate it.

    splin2 = splint(x1a,yytmp,y2tmp2,x1)

  end function splin2

  !---------------------------------------------------------------------
  !!! Wrap-up (integration) of 'splie2' and 'splin2' above
  function splint_2d(x1a,x2a,ya,x1,x2)
    
    implicit none

    real(p_double), dimension(:), intent(in) :: x1a,x2a
    real(p_double), dimension(:,:), intent(in) :: ya
    real(p_double), intent(in) :: x1,x2
    real(p_double), dimension(size(ya,1),size(ya,2)) :: y2a
    real(p_double) :: splint_2d

    integer :: j,m,ndum
    real(p_double), dimension(size(x1a)) :: yytmp,y2tmp2

    m    = size(x1a)
    ndum = size(x2a)

    ! first dimension
    do j=1,m ! Values 1 x 10^30 signal a natural spline.
      call spline(x2a,ya(j,:),1.0e30_p_double,1.0e30_p_double,y2a(j,:))
    end do
    
    do j=1,m
      yytmp(j) = splint(x2a,ya(j,:),y2a(j,:),x2)
      ! Perform m evaluations of the splines 
    end do
    
    ! second dimension
    call spline(x1a,yytmp,1.0e30_p_double,1.0e30_p_double,y2tmp2)
    ! Construct the one-dimensional spline 
    ! from the "reduced" row/column of points

    splint_2d = splint(x1a,yytmp,y2tmp2,x1)
    
  end function splint_2d




  !---------------------------------------------------------------------!
  !----- 3D ------------------------------------------------------------!
  !---------------------------------------------------------------------!
  function splint_3d(x1a,x2a,x3a,ya,x1,x2,x3)
    
    implicit none

    real(p_double), dimension(:), intent(in) :: x1a,x2a,x3a
    real(p_double), dimension(:,:,:), intent(in) :: ya
    real(p_double), intent(in) :: x1,x2,x3
    real(p_double), dimension(size(ya,1),size(ya,2),size(ya,3)) :: y3a
    real(p_double) :: splint_3d

    integer :: i,j,m,n,l
    real(p_double), dimension(size(x1a),size(x2a)) :: yytmp

    m = size(x1a)
    n = size(x2a)
    l = size(x3a)

    ! along first dimension (a "bundle" of 1d lines)
    do i=1,m ! Values 1 x 10^30 signal a natural spline.
      do j=1,n 
        call spline(x3a,ya(i,j,:),1.0e30_p_double,1.0e30_p_double,y3a(i,j,:))
      end do
    end do

    do i=1,m       ! evaluate the "plane" appropriate at x3
      do j=1,n
        yytmp(i,j) = splint(x3a,ya(i,j,:),y3a(i,j,:),x3)
      end do
    end do

    ! remaining 2 dimensions
    splint_3d = splint_2d(x1a,x2a,yytmp,x1,x2)

  end function splint_3d



  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  !-----      Gridding
  !---------------------------------------------------------------------!
  subroutine interp_grid1d(ddim_1d,bmin_1d,bmax_1d,dx,nc_x1,x1_box)

    implicit none

    integer, dimension(1), intent(in) :: ddim_1d
    real(p_double), dimension(1), intent(in) :: bmin_1d, bmax_1d
    real(p_double), dimension(1), intent(in) :: dx
    integer, intent(out) :: nc_x1
    real(p_double), allocatable, dimension(:), intent(out) :: x1_box

    !! Dummies
    real(p_double) :: bmin_idx, bmax_idx
    integer :: i

    ! In cell indices
    bmin_idx = bmin_1d(1) / dx(1)
    bmax_idx = bmax_1d(1) / dx(1)

    ! how many grids in new simulation of the desired box 
    nc_x1 = int(bmax_idx-bmin_idx)
      
    ! construct grid arrays
    !! fix first ends - namely the first cells of the new 
    !! simulation always correspond to the first cells of 
    !! the data
    allocate(x1_box(nc_x1))
    do i=1,nc_x1
      !! starting from 1, the first cell
      x1_box(i) = 1.0_p_double + (i-1)*(real(ddim_1d(1),p_double)-1.0_p_double)/real(nc_x1,p_double)
    end do

  end subroutine interp_grid1d

  !---------------------------------------------------------------------!
  subroutine interp_grid2d(ddim_2d,bmin_2d,bmax_2d,dx,nc_x1,nc_x2,x1_box,x2_box)

    implicit none

    integer, dimension(2), intent(in) :: ddim_2d
    real(p_double), dimension(2), intent(in) :: bmin_2d, bmax_2d
    real(p_double), dimension(2), intent(in) :: dx
    integer, intent(out) :: nc_x1, nc_x2
    real(p_double), allocatable, dimension(:), intent(out) :: x1_box, x2_box

    !! Dummies
    real(p_double), dimension(2) :: bmin_idx, bmax_idx
    integer :: i,j

    ! In cell indices; array operation
    bmin_idx = bmin_2d / dx
    bmax_idx = bmax_2d / dx

    ! how many grids in new simulation of the desired box 
    nc_x1 = int(bmax_idx(1)-bmin_idx(1))
    nc_x2 = int(bmax_idx(2)-bmin_idx(2))
      
    ! construct grid arrays - fix first ends
    !! see 'interp_grid1d' for details
    allocate(x1_box(nc_x1))
    allocate(x2_box(nc_x2))

    do i=1,nc_x1
      x1_box(i) = 1.0_p_double + (i-1)*(real(ddim_2d(1),p_double)-1.0_p_double)/real(nc_x1,p_double)
    end do
    do j=1,nc_x2
      x2_box(j) = 1.0_p_double + (j-1)*(real(ddim_2d(2),p_double)-1.0_p_double)/real(nc_x2,p_double)
    end do

  end subroutine interp_grid2d

  !---------------------------------------------------------------------!
  subroutine interp_grid3d(ddim_3d,bmin_3d,bmax_3d,dx, &
                           nc_x1,nc_x2,nc_x3,x1_box,x2_box,x3_box)

    implicit none

    integer, dimension(3), intent(in) :: ddim_3d
    real(p_double), dimension(3), intent(in) :: bmin_3d, bmax_3d
    real(p_double), dimension(3), intent(in) :: dx
    integer, intent(out) :: nc_x1, nc_x2, nc_x3
    real(p_double), allocatable, dimension(:), intent(out) :: x1_box, x2_box, x3_box

    !! Dummies
    real(p_double), dimension(3) :: bmin_idx, bmax_idx
    integer :: i,j,k

    ! In cell indices; array operation
    bmin_idx = bmin_3d / dx
    bmax_idx = bmax_3d / dx

    ! how many grids in new simulation of the desired box 
    nc_x1 = int(bmax_idx(1)-bmin_idx(1))
    nc_x2 = int(bmax_idx(2)-bmin_idx(2))
    nc_x3 = int(bmax_idx(3)-bmin_idx(3))
      
    ! construct grid arrays - fix first ends
    !! see 'interp_grid1d' for details
    allocate(x1_box(nc_x1))
    allocate(x2_box(nc_x2))
    allocate(x3_box(nc_x3))

    do i=1,nc_x1
      x1_box(i) = 1.0_p_double + (i-1)*(real(ddim_3d(1),p_double)-1.0_p_double)/real(nc_x1,p_double)
    end do
    do j=1,nc_x2
      x2_box(j) = 1.0_p_double + (j-1)*(real(ddim_3d(2),p_double)-1.0_p_double)/real(nc_x2,p_double)
    end do
    do k=1,nc_x3
      x3_box(k) = 1.0_p_double + (k-1)*(real(ddim_3d(3),p_double)-1.0_p_double)/real(nc_x3,p_double)
    end do

  end subroutine interp_grid3d


  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  !-----      Linear Interpolation implementations
  !---------------------------------------------------------------------!
  !----- 1D ------------------------------------------------------------!
  !---------------------------------------------------------------------!

  function lin_1d(xa,ya,x)

    ! This function uses bilinear interpolation to estimate the value
    ! of a function ya at point x
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by xa
    
    implicit none

    real(p_double), dimension(:), intent(in) :: xa,ya
    real(p_double), intent(in) :: x
    real(p_double) :: lin_1d
    
    real(p_double) :: xl,xh
    integer :: xj

    ! indices of x in xa
    xj = locate(xa,x)

    !!
    xl = xa(xj)
    xh = xa(xj+1)

    !
    lin_1d = ( ya(xj)*(xh-x) + ya(xj+1)*(x-xl) ) / (xh-xl)

  end function lin_1d



  !---------------------------------------------------------------------!
  !----- 2D ------------------------------------------------------------!
  !---------------------------------------------------------------------!

  function bilin_2d(x1a,x2a,ya,x1,x2)

    ! This function uses bilinear interpolation to estimate the value
    ! of a function ya at point (x1,x2)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x1a and the grid y values specified by x2a
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation; 
    !            https://github.com/cfinch/
    
    implicit none

    real(p_double), dimension(:), intent(in) :: x1a,x2a
    real(p_double), dimension(:,:), intent(in) :: ya
    real(p_double), intent(in) :: x1,x2

    real(p_double) :: denom,x1l,x1h,x2l,x2h
    integer :: x1j,x2j

    real(p_double) :: bilin_2d

    ! indices of x1,x2 in x1a,x2a
    x1j = locate(x1a,x1)
    x2j = locate(x2a,x2)

    !!
    x1l = x1a(x1j)
    x1h = x1a(x1j+1)

    x2l = x2a(x2j)
    x2h = x2a(x2j+1)

    !
    denom = (x1h - x1l) * (x2h - x2l)

    !
    bilin_2d = ( ya(x1j,x2j)*(x1h-x1)*(x2h-x2) + &
                 ya(x1j+1,x2j)*(x1-x1l)*(x2h-x2) + & 
                 ya(x1j,x2j+1)*(x1h-x1)*(x2-x2l) + &
                 ya(x1j+1,x2j+1)*(x1-x1l)*(x2-x2l) ) / denom
    
  end function bilin_2d

  !---------------------------------------------------------------------!
  !----- 3D ------------------------------------------------------------!
  !---------------------------------------------------------------------!

  

end module m_interp_fun




