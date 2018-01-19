!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Testing algorithms with linear interpolations !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_lin
  use m_interp_fun
  use nrtype

  implicit none

  integer, parameter :: nx = 8     !! number of points in x
  integer, parameter :: ny = 8     !! number of points in y
  integer, parameter :: nz = 8     !! number of points in z

  real(dp), dimension(nx) :: x
  real(dp), dimension(ny) :: y
  real(dp), dimension(nz) :: z
  real(dp), dimension(nx-1) :: xt
  real(dp), dimension(ny-1) :: yt
  real(dp), dimension(nz-1) :: zt

  ! function values on given grids / temporary storing spline output
  real(dp), dimension(nx)       :: fcn_1d, tmp_1d
  real(dp), dimension(nx,ny)    :: fcn_2d, tmp_2d
  real(dp), dimension(nx,ny,nz) :: fcn_3d, tmp_3d

  ! arrays for values on test grids (with 1 fewer index)
  !  true_o*d stores true function values from f1,f2,f3 defined below
  !  err_*d then stores the difference between true and interpolated values
  real(dp), dimension(nx-1)           :: out_1d, true_o1d, err_1d
  real(dp), dimension(nx-1,ny-1)      :: out_2d, true_o2d, err_2d
  real(dp), dimension(nx-1,ny-1,nz-1) :: out_3d, true_o3d, err_3d

  integer :: i,j,k

  ! Setup coordinate grids:
  !  integer values for data grid
  !  half integers for grid to be interpolated
  do i=1,nx-1
     x(i)  = dble(i-1)!/dble(nx-1)
     xt(i) = (dble(i)-0.5)!/dble(nx-1)
  end do
  x(nx) = dble(nx-1)!/dble(nx-1)

  do j=1,ny-1
     y(j)  = (dble(j-1))!/dble(ny-1)
     yt(j) = (dble(j)-0.5)!/dble(ny-1)
  end do
  y(ny) = dble(ny-1)!/dble(ny-1)

  do k=1,nz-1
     z(k)  = (dble(k-1))!/dble(nz-1)
     zt(k) = (dble(k)-0.5)!/dble(nz-1)
  end do
  z(nz) = dble(nz-1)!/dble(nz-1)

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
  !-- 2D
  !---------------------------------------------------------------------
  do i=1,nx-1
    do j=1,ny-1
      out_2d(i,j) = linint(x,y,fcn_2d,xt(i),yt(j))
      err_2d(i,j) = abs(out_2d(i,j)-true_o2d(i,j))

      if (mod(i,2)==0 .and. mod(j,2)==0) then
        write (*,*) 'At index (',i,',',j,') coordinates', xt(i),yt(j)
        write (*,*) ' value  is ', true_o2d(i,j), 'and'
        write (*,*) ' output is ', out_2d(i,j), 'and'
        write (*,*) ' error  is ', err_2d(i,j)/true_o2d(i,j) * 100, '%'
      end if
    end do
  end do

  !---------------------------------------------------------------------
  !-- Test spline - 3D
  !---------------------------------------------------------------------
  !do i=1,nx-1
  !  do j=1,ny-1
  !    do k=1,nz-1
  !      out_3d(i,j,k) = splint_3d(x,y,z,fcn_3d,xt(i),yt(j),zt(k))
  !      err_3d(i,j,k) = abs(out_3d(i,j,k)-true_o3d(i,j,k))

  !      if (mod(i,2)==0 .and. mod(j,2)==0 .and. mod(k,2)==0) then
          ! write (*,*) 'At index (',i,',',j,',',k,') coordinates', xt(i),yt(j),zt(k)
          ! write (*,*) ' value is ', true_o3d(i,j,k), ' and '
          ! write (*,*) ' error is ', err_3d(i,j,k)
  !      end if
  !    end do
  !  end do
  !end do

  !! Santiy check; inputing the exact same grid
  !do i=1,nx
  !  do j=1,ny
  !    do k=1,nz
  !      out_3d(i,j,k) = splint_3d(x,y,z,fcn_3d,x(i),y(j),z(k))
  !      err_3d(i,j,k) = abs(out_3d(i,j,k)-fcn_3d(i,j,k))

  !      if (mod(i,2)==0 .and. mod(j,2)==0 .and. mod(k,2)==0) then
          ! write (*,*) 'At coordinates', x(i),y(j),z(k)
          ! write (*,*) ' value is ', fcn_3d(i,j,k), ' and '
  !        ! write (*,*) ' error is ', err_3d(i,j,k)
  !      end if
  !    end do
  !  end do
  !end do





  !---------------------------------------------------------------------

  contains

    function f1(x) !! 1d test function
      implicit none
      real(dp) :: x,f1
      f1 = 0.5_dp * (x*exp(-x) + sin(x) )
    end function f1

    function f2(x,y) !! 2d test function
      implicit none
      real(dp) :: x,y,piov2,f2
      piov2 = 2.0_dp !* atan(1.0_dp)
      f2 = exp(-x/100.0_p_double) * sin(piov2*y/50.0_p_double) + cos((x/50.0_p_double)**0.5)
    end function f2
 
    function f3 (x,y,z) !! 3d test function
      implicit none
      real(dp) :: x,y,z,piov2,f3
      piov2 = 2.0_dp!*atan(1.0_dp)
      f3 = 0.5_dp*( y*exp(-x) + z*sin(piov2*y) )
    end function f3

end program test_lin
