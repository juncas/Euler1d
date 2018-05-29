subroutine init_Lax_1d
	use flow_1d
	implicit none
	integer :: i
	double precision :: rho,u,p,t
	xl=1.d0
	tl=0.13d0
	hx=xl/(il-1.d0)
	do i=1,il+2*bfsize
		x(i)=hx*(i-bfsize-1.d0)-xl/2.d0
	enddo
	do i=1,il+2*bfsize
		if(x(i).lt.0) then
			rho=0.445d0
			u=0.698d0
			p=3.528d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		else
			rho=0.5d0
			u=0.d0
			p=0.571d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		endif
	enddo
end subroutine
