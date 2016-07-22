! This is the module that generates the Laguerre Polynomials
module mod_Laguerre
implicit none

contains




! This function evaluates the associated Laguerre polynomial with an exponential at a given point.
! This function behaves quite well even for large values of x. The exponential helps keep the polynomial from
! overflow problems.
! LagExp = Exp(-1/2*x)*Lag[n,alpha,x]
real*8 function LagExp(n,alpha,x)
implicit none
integer,intent(in):: n
real*8,intent(in):: alpha
real*8,intent(in):: x
real*8::L,Lm,Lp
integer:: i 

! Initialize the values
Lm = 1.d0*dexp(-0.5d0*x)
L = (1.d0+ alpha -x)*dexp(-0.5d0*x)
Lp = 0.d0


    if(n.eq.0) then
        Lp = Lm
    else if (n.eq.1) then
        Lp = L
    else
    
        do i=1,n-1
        Lp = ((2.d0*dfloat(i)+1.d0+alpha-x)/(dfloat(i)+1.d0)*L - ((dfloat(i)+alpha)/(dfloat(i)+1.d0))*Lm)
        Lm = L
        L = Lp
        end do
    
    end if
    
    LagExp = Norm(n,alpha)*Lp

end function

! N = Dsqrt[Gamma[n+1]/Gamma[n+alpha+1]]
real*8 function Norm(n,alpha)
implicit none
integer::n
real(8)::alpha,logNorm

LogNorm = 0.5d0*(DLGAMA(dfloat(n+1)) - DLGAMA(dfloat(n+1)+alpha))
Norm = dexp(LogNorm)


end function Norm

! This function evaluates the associated Laguerre polynomial with an exponential at a given point.
! This function behaves quite well even for large values of x. The exponential helps keep the polynomial from
! overflow problems.
! LagExp = Exp(-1/2*x)*Lag[n,alpha,x]
real*8 function Orthog_LagExp(n,alpha,x)
implicit none
integer,intent(in):: n
real*8,intent(in):: alpha
real*8,intent(in):: x
real*8::L,Lm,Lp
integer:: i 

! Initialize the values
!Lm = 1.d0*dexp(-0.5d0*x)
Lm = 1.d0/dsqrt(Dgamma(alpha+1.d0))*dexp(-0.5d0*x)
L = (1.d0+ alpha -x)*dexp(-0.5d0*x)*(1.d0/dsqrt(dGamma(alpha+2.d0)))
Lp = 0.d0


    if(n.eq.0) then
        Lp = Lm
    else if (n.eq.1) then
        Lp = L
    else
    
        do i=1,n-1
        Lp = (((2.d0*dfloat(i)+1.d0+alpha-x)/(dfloat(i)+1.d0))*dsqrt((dfloat(i)+1)/(dfloat(i)+1.d0+alpha))*L &
        & - ((dfloat(i)+alpha)/(dfloat(i)+1.d0))*(dsqrt((dfloat(i)*dfloat(i+1))/((dfloat(i)+alpha)*(dfloat(i+1)+alpha))))*Lm)
        Lm = L
        L = Lp
        end do
    
    end if
    
    Orthog_LagExp = Lp

end function


end module
