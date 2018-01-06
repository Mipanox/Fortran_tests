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
  integer, dimension(1), parameter :: ddim_1d = 50
  integer, dimension(2), parameter :: ddim_2d = (/50, 60/)
  integer, dimension(3), parameter :: ddim_3d = (/50, 60, 1/)

  real(p_double), dimension(ddim_1d(1)) :: x1
  real(p_double), dimension(ddim_2d(2)) :: x2
  real(p_double), dimension(ddim_3d(3)) :: x3

  !! function values on given grids / temporary storing spline output
  real(p_double), dimension(ddim_1d(1))                       :: data_1d, tmp_1d
  real(p_double), dimension(ddim_2d(1),ddim_2d(2))            :: data_2d, tmp_2d
  real(p_double), dimension(ddim_3d(1),ddim_3d(2),ddim_3d(3)) :: data_3d, tmp_3d

  ! Parameters of the test grids
  !! Box extent in new setup
  real(p_double), dimension(1), parameter :: bmin_1d = 2.68, &
                                       bmax_1d = 40.
  real(p_double), dimension(2), parameter :: bmin_2d = (/2.68,  8.5/), &
                                       bmax_2d = (/40. , 36.5/)
  real(p_double), dimension(3), parameter :: bmin_3d = (/2.68,  8.5,  6.1/), &
                                       bmax_3d = (/40. , 36.5, 20.5/)

  !! Sizes of cells in simulation units
  real(p_double), dimension(1), parameter :: dx_1d = 0.9
  real(p_double), dimension(2), parameter :: dx_2d = (/0.9, 1.6/)
  real(p_double), dimension(3), parameter :: dx_3d = (/0.9, 1.6, 0.8/)

  !! Cell grids to be determined by box extent, etc.
  integer :: nc_x1, nc_x2, nc_x3
  real(p_double), allocatable, dimension(:) :: x1_box, x2_box, x3_box


  ! Variables for testing accuracy (true values at the new grid points)
  real(p_double), allocatable, dimension(:)     :: out_1d, tru_1d, err_1d
  real(p_double), allocatable, dimension(:,:)   :: out_2d, tru_2d, err_2d
  real(p_double), allocatable, dimension(:,:,:) :: out_3d, tru_3d, err_3d


  integer :: i,j,k
  integer, dimension(2) :: idx1,idx2,idx3
  real(p_double) :: tmp

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
  call interp_grid(ddim_1d,bmin_1d,bmax_1d,dx_1d,nc_x1,x1_box)

  do i=1,nc_x1
    ! write (*,*) x1_box(i)
  end do
  ! write (*,*) nc_x1


  !!-- Natural Spline
  call spline(x1,data_1d,1.0e30_p_double,1.0e30_p_double,tmp_1d)

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
  call interp_grid(ddim_2d,bmin_2d,bmax_2d,dx_2d,nc_x1,nc_x2,x1_box,x2_box)

  do i=1,nc_x1
    ! print *, x1_box(i)
  end do
  do j=1,nc_x2
    ! print *, x2_box(j)
  end do

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
      ! Only interpolating between nearby points
      idx1 = (/ mod(max(locate(x1,x1_box(i))-10,1),ddim_2d(1)-1), &
                min(locate(x1,x1_box(i))+10,ddim_2d(1))  /)
      idx2 = (/ mod(max(locate(x2,x2_box(j))-10,1),ddim_2d(2)-1), &
                min(locate(x2,x2_box(j))+10,ddim_2d(2))  /)
      
      out_2d(i,j) = splint_2d(x1(idx1(1):idx1(2)) ,x2(idx2(1):idx2(2)), &
                              data_2d(idx1(1):idx1(2),idx2(1):idx2(2)), &
                              x1_box(i),x2_box(j))
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
  !-- Test - 3D
  !---------------------------------------------------------------------

  !!-- Grid allocation
  call interp_grid(ddim_3d,bmin_3d,bmax_3d,dx_3d, &
                   nc_x1,nc_x2,nc_x3,x1_box,x2_box,x3_box)

  do i=1,nc_x1
    ! print *, x1_box(i)
  end do
  do j=1,nc_x2
    ! print *, x2_box(j)
  end do
  do k=1,nc_x3
    ! print *, x3_box(k)
  end do

  !! true values
  allocate(tru_3d(nc_x1,nc_x2,nc_x3))
  allocate(out_3d(nc_x1,nc_x2,nc_x3))
  allocate(err_3d(nc_x1,nc_x2,nc_x3))

  do i=1,nc_x1
    do j=1,nc_x2
      do k=1,nc_x3
        tru_3d(i,j,k) = f3(x1_box(i),x2_box(j),x3_box(k))
      end do
    end do
  end do
  
  !!-- Natural Spline
  !do i=1,nc_x1
  !  do j=1,nc_x2
  !    do k=1,nc_x3
  !      out_3d(i,j,k) = splint_3d(x1,x2,x3,data_3d,x1_box(i),x2_box(j),x3_box(k))
  !      err_3d(i,j,k) = abs(out_3d(i,j,k)-tru_3d(i,j,k))

  !      if (mod(i,5)==0 .and. mod(j,5)==0 .and. mod(k,10)==0) then
          !write (*,*) 'At index (',i,',',j,',',k,') coordinates', x1_box(i),x2_box(j),x3_box(k)
          !write (*,*) ' value  is ', tru_3d(i,j,k), 'and'
          !write (*,*) ' output is ', out_3d(i,j,k), 'and'
          !write (*,*) ' error  is ', err_3d(i,j,k)
  !      end if
  !    end do
  !  end do
  !end do

  deallocate(tru_3d)
  deallocate(out_3d)
  deallocate(err_3d)






  !---------------------------------------------------------------------

  contains

    !-------------------------------------------------------------
    function f1(x) !! 1d test function
      implicit none
      real(p_double) :: x,f1
      f1 = 0.5_p_double * (x*exp(-x) + sin(x) )
    end function f1

    function f2(x,y) !! 2d test function
      implicit none
      real(p_double) :: x,y,piov2,f2
      piov2 = 2.0_p_double !* atan(1.0_p_double)
      f2 = exp(-x/10.0_p_double) * sin(piov2*y/5.0_p_double) + cos((x/50.0_p_double)**0.5)
      !f2 = 0.5_p_double * (y/50.0_p_double*exp(-x/50.0_p_double) + sin(piov2*y/50.0_p_double) )
    end function f2
 
    function f3 (x,y,z) !! 3d test function
      implicit none
      real(p_double) :: x,y,z,piov2,f3
      piov2 = 2.0_p_double!*atan(1.0_p_double)
      f3 = 0.5_p_double*( y/50.0_p_double*exp(-x/100.0_p_double) + z/50.0_p_double*sin(piov2*y/50.0_p_double) )
    end function f3

end program test_os
