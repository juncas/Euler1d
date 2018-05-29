subroutine post_entropy
use flow_1d
implicit none
integer :: i
double precision :: rho,u,p,t,rho0,u0,p0,t0
double precision :: l2_rho,l2_u,l2_p,l2,linf

l2_rho=0.d0
l2_u=0.d0
l2_p=0.d0
l2=0.d0
linf=0.d0
do i=bfsize+1,il+bfsize
	rho0=1.d0 + 0.2d0*sin(pi*(x(i)-time))
	u0=1.d0
	p0=1.d0

	rho = q(i,1)
	u   = q(i,2)/rho
	p   = (q(i,3)-0.5d0*rho*u**2)*(gamma-1.d0)

	l2_rho  = l2_rho + (rho-rho0)**2
	l2_u    = l2_u   + (u-u0)**2
	l2_p    = l2_p   + (p-p0)**2


enddo
l2 = sqrt((l2_rho+l2_u+l2_p)/il)

open(10,file="entropy_error.dat",ACCESS='APPEND')
write(10,*),il,l2,sqrt(l2_rho/il),sqrt(l2_u/il),sqrt(l2_p/il)
close(10)
end 

