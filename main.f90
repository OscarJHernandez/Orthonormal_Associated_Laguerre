! This is the main program that tests out the Associted Laguerre polynomials
program main
use mod_Laguerre
use mod_integrate
implicit none
integer::n1,n2
real(8)::alpha1,alpha2

alpha1 = 0.5d0
alpha2 = 0.5d0

call initquad()

do n1=500,501
do n2=500,501

print *, n1,n2,integrate(n1,alpha1,n2,alpha2,LagExp),integrate(n1,alpha1,n2,alpha2,Orthog_LagExp)
!print *, n1,n2,LagExp(n1,alpha1,1.d0)

end do
end do



print *, 'Program Terminated'
end program
