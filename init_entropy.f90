subroutine init_entropy_1d
	use flow_1d
	implicit none
	integer :: i
	double precision :: rho,u,p,t
	xl=2.d0
	tl=2.0d0
	hx=xl/(il-1.d0)
	do i=1,il+2*bfsize
		x(i)=hx*(i-bfsize-1.d0)
	enddo
	do i=1,il+2*bfsize
			rho=1.d0 + 0.2d0*sin(pi*x(i))
			u=1.d0
			p=1.d0
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5d0*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
	enddo
end 
subroutine bndy_entropy_1d
	use flow_1d
	do i=1,bfsize
		q(i,1)=q(il+i-1,1)
		q(i,2)=q(il+i-1,2)
		q(i,3)=q(il+i-1,3)           

		q(il+bfsize+i,1)=q(bfsize+i+1,1)
		q(il+bfsize+i,2)=q(bfsize+i+1,2)
		q(il+bfsize+i,3)=q(bfsize+i+1,3)           
	end do
end
