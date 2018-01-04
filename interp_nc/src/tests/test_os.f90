!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Testing scenarios in reading in RAW data file !!
!-------------------------------------------------!
!! -- see `From box info to grid to be           !!
!!    interpolated on` in `hdf5_ext_int.ipynb`   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_os
  use m_interp_fun
  use nrtype

  implicit none

  ! Variables for setting up test data
  !! datadim
  integer, dimension(1), parameter :: ddim_1d = 80
  integer, dimension(2), parameter :: ddim_2d = (/80, 90/)
  integer, dimension(3), parameter :: ddim_3d = (/80, 90, 70/)

  real(dp), dimension(ddim_1d(1)) :: x1
  real(dp), dimension(ddim_2d(2)) :: x2
  real(dp), dimension(ddim_3d(3)) :: x3

  !! function values on given grids / temporary storing spline output
  real(dp), dimension(ddim_1d(1))                       :: data_1d, tmp_1d
  real(dp), dimension(ddim_2d(1),ddim_2d(2))            :: data_2d, tmp_2d
  real(dp), dimension(ddim_3d(1),ddim_3d(2),ddim_3d(3)) :: data_3d, tmp_3d

  ! Parameters of the test grids
  !! Box extent in new setup
  real(dp), dimension(1), parameter :: bmin_1d = 2.68, &
                                       bmax_1d = 40.
  real(dp), dimension(2), parameter :: bmin_2d = (/2.68,  8.5/), &
                                       bmax_2d = (/40. , 36.5/)
  real(dp), dimension(3), parameter :: bmin_3d = (/2.68,  8.5,  6.1/), &
                                       bmax_3d = (/40. , 36.5, 20.5/)

  !! Sizes of cells in simulation units
  real(dp), dimension(1), parameter :: dx_1d = 0.9
  real(dp), dimension(2), parameter :: dx_2d = (/0.9, 1.6/)
  real(dp), dimension(3), parameter :: dx_3d = (/0.9, 1.6, 0.8/)

  !! Cell grids to be determined by box extent, etc.
  integer :: nc_x1, nc_x2, nc_x3
  real(dp), allocatable, dimension(:) :: x1_box, x2_box, x3_box


  ! Variables for testing accuracy (true values at the new grid points)
  real(dp), allocatable, dimension(:)     :: out_1d, tru_1d, err_1d
  real(dp), allocatable, dimension(:,:)   :: out_2d, tru_2d, err_2d
  real(dp), allocatable, dimension(:,:,:) :: out_3d, tru_3d, err_3d


  integer :: i,j,k

  !! Setup coordinate grids:
  !!  integer values for data grid
  !!  half integers for grid to be interpolated
  do i=1,ddim_1d(1)
     x1(i) = dble(i)
  end do
  do j=1,ddim_2d(2)
     x2(j) = dble(j)
  end do
  do k=1,ddim_3d(3)
     x3(k) = dble(k)
  end do

  !! Fill in test arrays
  do i=1,ddim_1d(1)
    data_1d(i) = f1(x1(i))
    do j=1,ddim_2d(2)
      data_2d(i,j) = f2(x1(i),x2(j))
      do k=1,ddim_3d(3)
        data_3d(i,j,k) = f3(x1(i),x2(j),x3(k))
      end do
    end do
  end do
  

  !---------------------------------------------------------------------
  !-- Test - 1D
  !---------------------------------------------------------------------

  !!-- Grid allocation
  call interp_grid1d(ddim_1d,bmin_1d,bmax_1d,dx_1d,nc_x1,x1_box)

  do i=1,nc_x1
    ! write (*,*) x1_box(i)
  end do
  ! write (*,*) nc_x1


  !!-- Natural Spline
  call spline(x1,data_1d,1.0e30_dp,1.0e30_dp,tmp_1d)

  allocate(tru_1d(nc_x1))
  allocate(out_1d(nc_x1))
  allocate(err_1d(nc_x1))

  !! true values
  do i=1,nc_x1
    tru_1d(i) = f1(x1_box(i))
  end do

  do i=1,nc_x1
    out_1d(i) = splint(x1,data_1d,tmp_1d,x1_box(i))
    err_1d(i) = abs(out_1d(i)-tru_1d(i))

    if (mod(i,5) == 0) then
      ! write (*,*) 'At index ', i, ', coordinate', x1_box(i)
      ! write (*,*) ' value is ', tru_1d(i), ' and '
      ! write (*,*) ' error is ', err_1d(i)
    end if
  end do

  deallocate(tru_1d)
  deallocate(out_1d)
  deallocate(err_1d)

  !---------------------------------------------------------------------
  !-- Test - 2D
  !---------------------------------------------------------------------

  !!-- Grid allocation
  call interp_grid2d(ddim_2d,bmin_2d,bmax_2d,dx_2d,nc_x1,nc_x2,x1_box,x2_box)

  do i=1,nc_x1
    ! print *, x1_box(i)
  end do
  do j=1,nc_x2
    ! print *, x2_box(j)
  end do
  ! print *, nc_x1, nc_x2
  print *, splint_2d(x1,x2,data_2d,x1_box(5),x2_box(10))

  !! true values
  allocate(tru_2d(nc_x1,nc_x2))
  allocate(out_2d(nc_x1,nc_x2))
  allocate(err_2d(nc_x1,nc_x2))

  do i=1,nc_x1
    do j=1,nc_x2
      tru_2d(i,j) = f2(x1_box(i),x2_box(j))
    end do
  end do

  !!-- Natural Spline
  do i=1,nc_x1
    do j=1,nc_x2
      out_2d(i,j) = splint_2d(x1,x2,data_2d,x1_box(i),x2_box(j))
      err_2d(i,j) = abs(out_2d(i,j)-tru_2d(i,j))

      if (mod(i,5)==0 .and. mod(j,5)==0) then
        write (*,*) 'At index (',i,',',j,') coordinates', x1_box(i),x2_box(j)
        write (*,*) ' value  is ', tru_2d(i,j), 'and'
        write (*,*) ' output is ', out_2d(i,j), 'and'
        write (*,*) ' error  is ', err_2d(i,j)
      end if
    end do
  end do

  deallocate(tru_2d)
  deallocate(out_2d)
  deallocate(err_2d)

  !---------------------------------------------------------------------

  contains

    !-------------------------------------------------------------
    function f1(x) !! 1d test function
      implicit none
      real(dp) :: x,f1
      f1 = 0.5_dp * (x*exp(-x) + sin(x) )
    end function f1

    function f2(x,y) !! 2d test function
      implicit none
      real(dp) :: x,y,piov2,f2
      piov2 = 2.0_dp !* atan(1.0_dp)
      f2 = exp(-x/100.0_dp) * sin(piov2*y/50.0_dp) + cos((x/50.0_dp)**0.5)
      !f2 = 0.5_dp * (y/50.0_dp*exp(-x/50.0_dp) + sin(piov2*y/50.0_dp) )
    end function f2
 
    function f3 (x,y,z) !! 3d test function
      implicit none
      real(dp) :: x,y,z,piov2,f3
      piov2 = 2.0_dp!*atan(1.0_dp)
      f3 = 0.5_dp*( y/50.0_dp*exp(-x/100.0_dp) + z/50.0_dp*sin(piov2*y/50.0_dp) )
    end function f3

end program test_os
