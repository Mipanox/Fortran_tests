program test
  use spl_fun
  use nrtype

  implicit none

  integer, parameter :: nx = 6     !! number of points in x
  integer, parameter :: ny = 6     !! number of points in y
  integer, parameter :: nz = 6     !! number of points in z

  real(SP) :: x(nx) ,y(ny) ,z(nz)
  real(SP) :: xt(nx-1),yt(nx-1),zt(nz-1)

  ! function values on given grids / temporary storing spline output
  real(SP) :: fcn_1d(nx)      , tmp_1d(nx)
  real(SP) :: fcn_2d(nx,ny)   , tmp_1d(nx,ny)
  real(SP) :: fcn_3d(nx,ny,nz), tmp_1d(nx,ny,nz)

  ! arrays for values on test grids (with 1 fewer index)
  !  true_o*d stores true function values from f1,f2,f3 defined below
  !  err_*d then stores the difference between true and interpolated values
  real(SP) :: out_1d(nx-1)          ,true_o1d(nx-1)          , err_1d(nx-1)
  real(SP) :: out_2d(nx-1,ny-1)     ,true_o2d(nx-1,ny-1)     , err_2d(nx-1,ny-1)
  real(SP) :: out_3d(nx-1,ny-1,nz-1),true_o3d(nx-1,ny-1,nz-1), err_3d(nx-1,ny-1,nz-1)

  integer :: i,j,k

  ! Setup coordinate grids:
  !  integer values for data grid
  !  half integers for grid to be interpolated
  do i=1,nx-1
     x(i)  = dble(i-1)
     xt(i) = dble(i)-0.5
  end do
  x(nx) = dble(nx-1)

  do j=1,ny-1
     y(j)  = dble(j-1)
     yt(i) = dble(j)-0.5
  end do
  y(ny) = dble(ny-1)

  do k=1,nz-1
     z(k)  = dble(k-1)
     zt(k) = dble(k)-0.5
  end do
  z(nz) = dble(nz-1)

  ! Setup test functions
  do i=1,nx
    fcn_1d(i) = f1(x(i))
    do j=1,ny
      fcn_2d(i,j) = f2(x(i),y(j))
      do k=1,nz
        fcn_3d(i,j,k) = f3(x(i),y(j),z(k))
      end do
    end do
  end do

  ! Evaluate "true" values for half-integral grids
  do i=1,nx-1
    true_o1d(i) = f1(xt(i))
    do j=1,ny-1
      true_o2d(i,j) = f2(xt(i),yt(j))
      do k=1,nz-1
        true_o3d(i,j,k) = f3(xt(i),yt(j),zt(k))
      end do
    end do
  end do

  !---------------------------------------------------------------------
  !-- Test spline
  !---------------------------------------------------------------------
  call spline(x,fcn_1d,1.0e30_sp,1.0e30_sp,tmp_1d)

  do i=1,nx-1
    out_1d(i) = splint(x,fcn_1d,tmp_1d,xt(i))
    err_1d(i) = abs(out_1d(i)-true_o1d(i))

    write (*,*) 'Error at index ', i, ' is ', err_1d(i)
  
  end do



  !---------------------------------------------------------------------

  contains

  function f1(x) !! 1d test function
    implicit none
    real(SP) :: x
    f1 = 0.5_sp * (x*exp(-x) + sin(x) )
  end function f1

  function f2(x,y) !! 2d test function
    implicit none
    real(SP) x,y,piov2
    piov2 = 2.0_sp * atan(1.0_sp)
    f2 = 0.5_sp * (y*exp(-x) + sin(piov2*y) )
  end function f2

end program test
