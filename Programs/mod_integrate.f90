! This module contains the main integration methods

module mod_integrate
implicit none
integer,parameter::NquadMomentum =100
integer,parameter::NquadPosition =5000
real(8),parameter::Xmax = 1000.d0
real(8),allocatable::daP(:),dbP(:),dp(:),dwP(:),eP(:)
real(8),allocatable::daX(:),dbX(:),dx(:),dwX(:),eX(:)
integer,parameter:: ipolyMomentum=2 !  Legendre polynomial from 0 to 1 
integer,parameter:: ipolyPosition=2 !  Legendre polynomial from 0 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe

contains

! Initialize the Qudrature routines
subroutine initquad()
! Allocate the Gauss-Legendre quadrature routine
    allocate(daP(NquadMomentum),dbP(NquadMomentum),dp(NquadMomentum),dwP(NquadMomentum),eP(NquadMomentum))
    allocate(daX(NquadPosition),dbX(NquadPosition),dx(NquadPosition),dwX(NquadPosition),eX(NquadPosition))
    
    al = 0.d0
    
    call drecur(NquadMomentum,ipolyMomentum,al,dbe,daP,dbP,iderr)
    call dgauss(NquadMomentum,daP,dbP,depsma,dp,dwP,ierr,eP)
    
    call drecur(NquadPosition,ipolyPosition,al,dbe,daX,dbX,iderr)
    call dgauss(NquadPosition,daX,dbX,depsma,dx,dwX,ierr,eX)
end subroutine

! This integrates our function for us
! Computes the overlap between two functions
real(8) function integrate(n1,alpha1,n2,alpha2,func)
implicit none
integer::n1,n2
real(8)::s,alpha1,alpha2,xi,f
real(8),external::func
integer::i

s=0.d0
do i=1,NquadPosition
xi = dx(i)*Xmax
f = func(n1,alpha1,xi)*func(n2,alpha2,xi)*(xi**(0.5d0*(alpha1+alpha2)))
s = s +f*dwX(i)*Xmax
end do 

integrate = s

end function



end module
