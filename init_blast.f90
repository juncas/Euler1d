subroutine init_blast_1d
	use flow_1d
	implicit none
	integer :: i
	double precision :: rho,u,p,t
	xl=1.d0
	tl=0.038d0
	hx=xl/(il-1.d0)
	do i=1,il+2*bfsize
		x(i)=hx*(i-bfsize-1.d0)
	enddo
	do i=1,il+2*bfsize
		if(x(i).lt.0.1) then
			rho=1
			u=0
			p=1000
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		elseif(x(i).gt.0.1.and.x(i).lt.0.9) then
			rho=1
			u=0
			p=0.1
			q(i,1)=rho
			q(i,2)=rho*u
			q(i,3)=p/(gamma-1.d0)+0.5*rho*u**2
			flux(i,1)=rho*u
			flux(i,2)=rho*u*u+p
			flux(i,3)=u*q(i,3)+u*p
		elseif(x(i).gt.0.9) then
			rho=1
			u=0
			p=100
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

subroutine bndy_blast_1d
	use flow_1d
	do i=1,bfsize
		q(i,1)=q(2*bfsize+1-i,1)
		q(i,2)=-q(2*bfsize+1-i,2)
		q(i,3)=q(2*bfsize+1-i,3)           

		q(il+bfsize+i,1)=q(il+bfsize+1-i,1)
		q(il+bfsize+i,2)=-q(il+bfsize+1-i,2)
		q(il+bfsize+i,3)=q(il+bfsize+1-i,3)           
	end do
end