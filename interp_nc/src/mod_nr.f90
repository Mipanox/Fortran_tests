! Only a small part of the nrutil described in Appendix C2 in the book
!! The file contains explicit interfaces for all the Numerical Recipes
!! routines (except those already in the module nrutil). The interfaces
!! are in alphabetical order, by the generic interface name, if one 
!! exists, or by the specific routine name if there is no generic name.

module nr
  ! On a purely serial machine, for greater efficiency, remove the 
  ! generic name tridag from the following interface, and put it on 
  ! the next one after that.
  interface tridag
    recursive subroutine tridag_par(a,b,c,r,u)
      use nrtype
      real(SP), dimension(:), intent(in) :: a,b,c,r
      real(SP), dimension(:), intent(out) :: u
    end subroutine tridag_par
  end interface

  interface 
    subroutine tridag_ser(a,b,c,r,u)
      use nrtype
      real(SP), dimension(:), intent(in) :: a,b,c,r
      real(SP), dimension(:), intent(out) :: u
    end subroutine tridag_ser
  end interface

  interface
    function locate(xx,x)
      use nrtype
      real(SP), dimension(:), intent(in) :: xx
      real(SP), intent(in) :: x
      integer(I4B) :: locate
    end function locate
  end interface

  interface 
    subroutine spline(x,y,yp1,ypn,y2)
      use nrtype
      real(SP), dimension(:), intent(in) :: x,y
      real(SP), intent(in) :: yp1, ypn
      real(SP), dimension(:), intent(out) :: y2
    end subroutine spline
  end interface

  interface
    function splint(xa,ya,y2a,x)
      use nrtype
      real(SP), dimension(:), intent(in) :: xa,ya,y2a
      real(SP), intent(in) :: x
      real(SP) :: splint
    end function splint
  end interface

end module nr
