subroutine init_contact_1d
	use flow_1d
	implicit none
	integer :: i
	double precision :: rho,u,p,t
	xl=10.d0
	tl=1.5d0
	hx=xl/(il-1.d0)
	do i=1,il+2*bfsize
		x(i)=hx*(i-bfsize-1.d0)-5.d0
	enddo
	do i=1,il+2*bfsize
		if(x(i).lt.0) then
			rho=0.125d0
			u=0.112
			p=1.d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		else
			rho=10.d0
			u=0.112d0
			p=1.d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		endif
	enddo
	close(10)

end subroutine
