! Only a small part of the nrutil described in Appendix C1.2 in the book
!! The file contains a single module named nrutil, which contains specific
!! implementations for all the Numerical Recipes utility functions described
!! in detail in Chapter 23. (details see p.1364 of the book)

module nrutil
  use nrtype

  ! Parameters for crossover from serial to parallel algorithms 
  ! (these are used only within thin nrutil module):
  implicit none
  integer(I4B), parameter :: NPAR_ARTH = 16, NPAR2_ARTH = 8
  integer(I4B), parameter :: NPAR_GEOP = 4 , NPAR2_GEOP = 2
  integer(I4B), parameter :: NPAR_CUMSUM = 16
  integer(I4B), parameter :: NPAR_CUMPROD = 8
  integer(I4B), parameter :: NPAR_POLY = 8
  integer(I4B), parameter :: NPAR_POLYTERM = 8

  ! Next, generic interfaces for routines with overloaded versions. 
  ! Naming conventions for appended codes in the names of overloaded 
  ! routines are as follows: r=real, d=double precision, i=integer,
  ! c=complex, z=double-precision complex, h=character, l=logical.
  ! Any of r,d,i,c,z,h,l may be followed by v=vector or m=matrix
  ! (v,m suffixes are used only when needed to resolve ambguities).

  !! LARGELY OMITTED AT THE MOMENT

  ! Routines for argument checking and error handling (nrerror is not
  ! currently overloaded and so do not have a generic interface here):
  interface assert_eq
    module procedure assert_eq2, assert_eq3, assert_eq4, assert_eqn
  end interface

  contains

  function assert_eq2(n1,n2,string) 
    ! Report and die if integers not all equal (used for size checking)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1, n2
    integer :: assert_eq2

    if (n1 == n2) then
      assert_eq2 = n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
      STOP 'program terminated by assert_eq2'
    end if
  end function assert_eq2

  function assert_eq3(n1,n2,n3,string) 
    ! Report and die if integers not all equal (used for size checking)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1, n2, n3
    integer :: assert_eq3

    if (n1 == n2 .and. n2 == n3) then
      assert_eq3 = n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
      STOP 'program terminated by assert_eq3'
    end if
  end function assert_eq3

  function assert_eq4(n1,n2,n3,n4,string) 
    ! Report and die if integers not all equal (used for size checking)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1, n2, n3, n4
    integer :: assert_eq4

    if (n1 == n2 .and. n2 == n3 .and. n3==n4) then
      assert_eq4 = n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
      STOP 'program terminated by assert_eq4'
    end if
  end function assert_eq4

  function assert_eqn(nn,string) 
    ! Report and die if integers not all equal (used for size checking)
    character(len=*), intent(in) :: string
    integer, dimension(:), intent(in) :: nn
    integer :: assert_eqn

    if ( all(nn(2:) == nn(1)) ) then
      assert_eqn = nn(1)
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
      STOP 'program terminated by assert_eqn'
    end if
  end function assert_eqn

  subroutine nrerror(string)
    ! Report a message, then die
    character(len=*), intent(in) :: string
    write (*,*) 'nrerror: ', string
    STOP 'program terminated by nrerror'
  end subroutine nrerror

end module nrutil
