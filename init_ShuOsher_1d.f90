subroutine init_ShuOsher_1d
	use flow_1d
	implicit none
	integer :: i
	double precision :: rho,u,p,t
	xl=10.d0
	tl=1.8d0
	hx=xl/(il-1.d0)
	do i=1,il+2*bfsize
		x(i)=hx*(i-bfsize-1.d0)-xl/2.d0
	enddo
	do i=1,il+2*bfsize
		if(x(i).lt.-4.d0) then
			rho=27.d0/7.d0
			u=4.d0*sqrt(35.d0)/9.d0
			p=31.d0/3.d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		else
			rho=1.d0+0.2*sin(5.d0*x(i))
			u=0.d0
			p=1.d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		endif
	enddo
end subroutine
