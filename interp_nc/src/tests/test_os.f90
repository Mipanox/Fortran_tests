!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Testing scenarios in reading in RAW data file !!
!-------------------------------------------------!
!! -- see `From box info to grid to be           !!
!!    interpolated on` in `hdf5_ext_int.ipynb`   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_os
  use spl_fun
  use nrtype

  implicit none

  ! Variables for setting up test data
  !! datadim
  integer, dimension(1), parameter :: ddim_1d = 80
  integer, dimension(2), parameter :: ddim_2d = (/80, 90/)
  integer, dimension(3), parameter :: ddim_3d = (/80, 90, 70/)

  real(SP), dimension(ddim_1d(1)) :: x1
  real(SP), dimension(ddim_2d(2)) :: x2
  real(SP), dimension(ddim_3d(3)) :: x3

  !! function values on given grids / temporary storing spline output
  real(SP), dimension(ddim_1d(1))                       :: data_1d, tmp_1d
  real(SP), dimension(ddim_2d(1),ddim_2d(2))            :: data_2d, tmp_2d
  real(SP), dimension(ddim_3d(1),ddim_3d(2),ddim_3d(3)) :: data_3d, tmp_3d

  ! Parameters of the test grids
  !! Box extent in new setup
  real(SP), dimension(1), parameter :: bmin_1d = 2.68, &
                                       bmax_1d = 40.
  real(SP), dimension(2), parameter :: bmin_2d = (/2.68,  8.5/), &
                                       bmax_2d = (/40. , 36.5/)
  real(SP), dimension(3), parameter :: bmin_3d = (/2.68,  8.5,  6.1/), &
                                       bmax_3d = (/40. , 36.5, 20.5/)

  !! Sizes of cells in simulation units
  real(SP), dimension(1), parameter :: dx_1d = 0.9
  real(SP), dimension(2), parameter :: dx_2d = (/0.9, 1.6/)
  real(SP), dimension(3), parameter :: dx_3d = (/0.9, 1.6, 0.8/)

  !! Cell grids to be determined by box extent, etc.
  integer :: nc_x1, nc_x2, nc_x3
  real(SP), allocatable, dimension(:) :: x1_box, x2_box, x3_box


  ! Variables for testing accuracy (true values at the new grid points)
  real(SP), allocatable, dimension(:)     :: out_1d, tru_1d, err_1d
  real(SP), allocatable, dimension(:,:)   :: out_2d, tru_2d, err_2d
  real(SP), allocatable, dimension(:,:,:) :: out_3d, tru_3d, err_3d


  integer :: i,j,k

  !! Setup coordinate grids:
  !!  integer values for data grid
  !!  half integers for grid to be interpolated
  do i=1,ddim_1d(1)
     x1(i)  = dble(i)
  end do
  do j=1,ddim_2d(2)
     x2(j)  = (dble(j))
  end do
  do k=1,ddim_3d(3)
     x3(k)  = (dble(k))
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
  call spline(x1,data_1d,1.0e30_sp,1.0e30_sp,tmp_1d)

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
      write (*,*) 'At index ', i, ', coordinate', x1_box(i)
      write (*,*) ' value is ', tru_1d(i), ' and '
      write (*,*) ' error is ', err_1d(i)
    end if
  end do



  !---------------------------------------------------------------------

  contains

    subroutine interp_grid1d(ddim_1d,bmin_1d,bmax_1d,dx,nc_x1,x1_box)
      implicit none
      integer, dimension(1), intent(in) :: ddim_1d
      real(SP), dimension(1), intent(in) :: bmin_1d, bmax_1d
      real(SP), dimension(1), intent(in) :: dx
      integer, intent(out) :: nc_x1
      real(SP), allocatable, dimension(:), intent(out) :: x1_box

      !! Dummies
      real(SP) :: bmin_idx, bmax_idx

      ! In cell indices
      bmin_idx = bmin_1d(1) / dx(1)
      bmax_idx = bmax_1d(1) / dx(1)

      ! how many grids in new simulation of the desired box 
      nc_x1 = int(bmax_idx-bmin_idx)
      
      ! construct grid arrays; using FORTRAN-like style
      !! fix first ends - namely the first cells of the new 
      !! simulation always correspond to the first cells of 
      !! the data
      allocate(x1_box(nc_x1))
      do i=1,nc_x1
        !! starting from 1, the first cell
        x1_box(i) = 1.0_sp + (i-1)*(real(ddim_1d(1),SP)-1.0_sp)/real(nc_x1,SP)
      end do

    end subroutine interp_grid1d



    !-------------------------------------------------------------
    function f1(x) !! 1d test function
      implicit none
      real(SP) :: x,f1
      f1 = 0.5_sp * (x*exp(-x) + sin(x) )
    end function f1

    function f2(x,y) !! 2d test function
      implicit none
      real(SP) :: x,y,piov2,f2
      piov2 = 2.0_sp * atan(1.0_sp)
      f2 = 0.5_sp * (y*exp(-x) + sin(piov2*y) )
    end function f2
 
    function f3 (x,y,z) !! 3d test function
      implicit none
      real(SP) :: x,y,z,piov2,f3
      piov2 = 2.0_sp*atan(1.0_sp)
      f3 = 0.5_sp*( y*exp(-x) + z*sin(piov2*y) )
    end function f3

end program test_os
