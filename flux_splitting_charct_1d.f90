subroutine flux_splitting_chrt_1d
use flow_1d
implicit none
double precision rho,u,p,c0,h0
integer :: i,j,k,m,info,n
double precision,dimension(nl,nl,il+2*bfsize) :: l,r,omega,a,b,c
double precision,dimension(nl,il+2*bfsize) :: d,fluxl,fluxr,alpha,alpha_mean
double precision :: rhor,rhol,ur,ul,pr,pl
double precision :: w(nl,2),temp1(nl),temp2(nl),temp(nl)
double precision :: l1,l2,alpha0

exmax=0.d0
do i=1,il+2*bfsize
	rho=q(i,1)
	u=q(i,2)/rho
	p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*(u**2))
	c0=sqrt(gamma*p/rho)
	exmax=max(exmax,abs(u)+c0)
	flux(i,1)=rho*u
	flux(i,2)=rho*u*u+p
	flux(i,3)=u*q(i,3)+u*p
enddo

do i=1,il+2*bfsize !LF splitting
	do n=1,nl
		fluxp(i,n)=0.5*(flux(i,n)+exmax*q(i,n))
		fluxn(i,n)=0.5*(flux(i,n)-exmax*q(i,n))
	enddo
enddo

do i=bfsize,il+bfsize+1
	!construction L and R with Roe mean
	rhor=q(i+1,1)
	ur=q(i+1,2)/rhor
	pr=(gamma-1.d0)*(q(i+1,3)-0.5d0*rhor*(ur**2))        

	rhol=q(i,1)        
	ul=q(i,2)/rhol
	pl=(gamma-1.d0)*(q(i,3)-0.5d0*rhol*(ul**2))

	c0=sqrt(gamma*pl/rhol)
	alpha(1,i)=ul
	alpha(2,i)=ul-c0
	alpha(3,i)=ul+c0

	w(1,1)=sqrt(rhor)
	w(1,2)=sqrt(rhol)
	w(2,1)=w(1,1)*ur
	w(2,2)=w(1,2)*ul
	w(3,1)=(q(i+1,3)+pr)/w(1,1)
	w(3,2)=(q(i,3)+pl)/w(1,2)

	u=(w(2,1)+w(2,2))/(w(1,1)+w(1,2))
	h0=(w(3,1)+w(3,2))/(w(1,1)+w(1,2))
	c0=sqrt((gamma-1)*(h0-0.5d0*u**2))

	alpha_mean(1,i)=u
	alpha_mean(2,i)=u-c0
	alpha_mean(3,i)=u+c0

	r(1,1,i)=1
	r(1,2,i)=1
	r(1,3,i)=1
	r(2,1,i)=u
	r(2,2,i)=(u-c0)
	r(2,3,i)=(u+c0)
	r(3,1,i)=0.5*u**2
	r(3,2,i)=(h0-u*c0)
	r(3,3,i)=(h0+u*c0)

	l1=(gamma-1)/c0**2
	l2=0.5*u**2

	l(1,1,i)=1.d0-l2*l1
	l(1,2,i)=u*l1
	l(1,3,i)=-l1
	l(2,1,i)=0.5d0*(u+l1*c0*l2)/c0
	l(2,2,i)=-0.5d0*(1.d0/c0+l1*u)
	l(2,3,i)=0.5d0*l1
	l(3,1,i)=0.5d0*(-u/c0+l1*l2)
	l(3,2,i)=0.5d0*(1.d0/c0-l1*u)
	l(3,3,i)=0.5d0*l1
	! this part can be build in one single subroutine
	!--------------------------------------------------------------------------------------------------
	!construction omega with L_1/2 and q
enddo
!construct D,E,F,b
select case(ischeme_inv)
case(2)
	call scheme_weno7th_1d_character_l(l,r,fluxp,fluxl,nl,il,bfsize)
	call scheme_weno7th_1d_character_r(l,r,fluxn,fluxr,nl,il,bfsize)
	do i=bfsize,il+bfsize
		do j=1,nl
			h(i,j)=fluxl(j,i)+fluxr(j,i)
		enddo
	enddo
case(1)
	call scheme_weno5th_1d_character_l(l,r,fluxp,fluxl,nl,il,bfsize)
	call scheme_weno5th_1d_character_r(l,r,fluxn,fluxr,nl,il,bfsize)
	do i=bfsize,il+bfsize
		do j=1,nl
			h(i,j)=fluxl(j,i)+fluxr(j,i)
		enddo
	enddo
end select
end subroutine

