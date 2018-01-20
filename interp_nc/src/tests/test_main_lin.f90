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
  real(dp), dimension(nx)       :: fcn_1d, tmp_1d, out_1d_, err_1d_
  real(dp), dimension(nx,ny)    :: fcn_2d, tmp_2d
  real(dp), dimension(nx,ny,nz) :: fcn_3d, tmp_3d, out_3d_, err_3d_

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
     xt(i) = (dble(i)-0.03)!/dble(nx-1)
  end do
  x(nx) = dble(nx-1)!/dble(nx-1)

  do j=1,ny-1
     y(j)  = (dble(j-1))!/dble(ny-1)
     yt(j) = (dble(j)-0.03)!/dble(ny-1)
  end do
  y(ny) = dble(ny-1)!/dble(ny-1)

  do k=1,nz-1
     z(k)  = (dble(k-1))!/dble(nz-1)
     zt(k) = (dble(k)-0.03)!/dble(nz-1)
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

  
  ! It is not wise to test on the grid points (as done in spline tests)
  !  since the denominator(s) in linear interpolations would go to infinity

  !---------------------------------------------------------------------
  !-- 1D
  !---------------------------------------------------------------------
  do i=1,nx-1
    out_1d(i) = linint(x,fcn_1d,xt(i))
    err_1d(i) = abs(out_1d(i)-true_o1d(i))

    ! write (*,*) 'At index (',i,') coordinate', xt(i)
    ! write (*,*) ' value  is ', true_o1d(i), 'and'
    ! write (*,*) ' output is ', out_1d(i), 'and'
    ! write (*,*) ' error  is ', err_1d(i)/true_o1d(i)*100, '%'

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
  !-- 3D
  !---------------------------------------------------------------------
  do i=1,nx-1
    do j=1,ny-1
      do k=1,nz-1
        out_3d(i,j,k) = linint(x,y,z,fcn_3d,xt(i),yt(j),zt(k))
        err_3d(i,j,k) = abs(out_3d(i,j,k)-true_o3d(i,j,k))

        if (mod(i,2)==0 .and. mod(j,2)==0 .and. mod(k,2)==0) then
          write (*,*) 'At index (',i,',',j,',',k,') coordinates', xt(i),yt(j),zt(k)
          write (*,*) ' value  is ', true_o3d(i,j,k), ' and '
          write (*,*) ' output is ', out_3d(i,j,k), 'and'
          write (*,*) ' error  is ', err_3d(i,j,k)/true_o3d(i,j,k)*100, '%'
        end if
      end do
    end do
  end do


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
      f3 = 0.5_dp*( y/50.0_p_double*exp(-x/50.0_p_double) + z*sin(piov2*y/50.0_p_double) )
    end function f3

end program test_lin
