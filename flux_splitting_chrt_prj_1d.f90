subroutine flux_splitting_chrt_prj_1d
use flow_1d
implicit none
double precision,dimension(nl,il+2*bfsize) :: prim
double precision,dimension(nl,il+2*bfsize) :: d,fluxl,fluxr

double precision :: rho,u,p,alpha0
integer :: i,j,n

exmax=0.d0

do i=1,il+2*bfsize

	rho=q(i,1)
	u=q(i,2)/rho
	p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*(u**2))
	alpha0=sqrt(gamma*p/rho)

	exmax=max(exmax,abs(u)+alpha0)

	flux(i,1)=rho*u
	flux(i,2)=rho*u*u+p
	flux(i,3)=u*q(i,3)+u*p

	prim(1,i)=rho
	prim(2,i)=u
	prim(3,i)=p
enddo


do i=1,il+2*bfsize
	do n=1,nl
		fluxp(i,n)=0.5d0*(flux(i,n)+exmax*q(i,n))
		fluxn(i,n)=0.5d0*(flux(i,n)-exmax*q(i,n))
	enddo
enddo

call get_tot_l(fluxp,flag_l,theta_l,il,bfsize)
call scheme_weno5th_tot_1d_character_l(theta_l,flag_l,prim,fluxp,fluxl,nl,il,bfsize)

call get_tot_r(fluxn,flag_r,theta_r,il,bfsize)
call scheme_weno5th_tot_1d_character_r(theta_r,flag_r,prim,fluxn,fluxr,nl,il,bfsize)

do i=bfsize,il+bfsize
	do j=1,nl
		h(i,j)=fluxl(j,i)+fluxr(j,i)
	enddo
enddo
end subroutine

module weno5_coef
implicit none
double precision,parameter :: a300 = 1.d0/3.d0
double precision,parameter :: a301 =-7.d0/6.d0
double precision,parameter :: a302 = 11.d0/6.d0
double precision,parameter :: a310 =-1.d0/6.d0
double precision,parameter :: a311 = 5.d0/6.d0
double precision,parameter :: a312 = 1.d0/3.d0
double precision,parameter :: a320 = 1.d0/3.d0
double precision,parameter :: a321 = 5.d0/6.d0
double precision,parameter :: a322 =-1.d0/6.d0
double precision,parameter :: c30  = 0.1d0
double precision,parameter :: c31  = 0.6d0
double precision,parameter :: c32  = 0.3d0
end module

subroutine get_tot_l(flux,flag,tot,il,bfsize)
use flow_1d,only:ss,u2,hx,switch_l
use weno5_coef
implicit none
integer,intent(in) :: il,bfsize
double precision,intent(in),dimension(il+2*bfsize,3) :: flux
double precision,intent(out),dimension(3,il+2*bfsize) :: tot
logical,intent(out),dimension(il+2*bfsize) :: flag


double precision :: is0,is1,is2,tau5
double precision :: aa0,aa1,aa2
double precision :: w0,w1,w2
double precision :: tot0
double precision :: lf(-2:2)
integer :: i,m,n

u2=0.0
do n=1,3
	do i=1,il+2*bfsize
		u2(n) = u2(n) + flux(i,n)**2*hx
	enddo
enddo 
do n=1,3
	u2(n) = sqrt(u2(n))
	if(u2(n).lt.1.0e-20) u2(n)=1.0
enddo 

tot = 0.0
do n=1,3
	do i=bfsize,il+bfsize
	        do m=-2,2
	            lf(m)=flux(i+m,n)
	        enddo
	        is0= 13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0 &
	            +(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
	        is1= 13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0 &
	        	+(lf(-1)-lf(1))**2/4.d0
	        is2= 13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0 &
	        	+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

	        tot(1,i)=tot(1,i) + is0/u2(n)/3.0
	        tot(2,i)=tot(2,i) + is1/u2(n)/3.0
	        tot(3,i)=tot(3,i) + is2/u2(n)/3.0
	enddo
enddo

do i=bfsize,il+bfsize
	is0 = tot(1,i)
	is1 = tot(2,i)
	is2 = tot(3,i)
	switch_l(i) = is0 + is1 + is2
    if(switch_l(i).gt.1.d0) then
    	flag(i) = .true.
    else
    	flag(i) = .false.
    endif

	tau5 = abs(is2-is1)
    aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
    aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
    aa2=c32*(1.0d0+(tau5/(ss+is2))**2)

	! aa0=c30/(ss+is0)**2
	! aa1=c31/(ss+is1)**2
	! aa2=c32/(ss+is2)**2
	w0=aa0/(aa0+aa1+aa2)
	w1=aa1/(aa0+aa1+aa2)
	w2=aa2/(aa0+aa1+aa2)    

	tot(1,i)=w0
	tot(2,i)=w1
	tot(3,i)=w2

enddo
end subroutine

subroutine get_tot_r(flux,flag,tot,il,bfsize)
use flow_1d,only:ss,u2,hx,switch_r
use weno5_coef
implicit none
integer,intent(in) :: il,bfsize
double precision,intent(in),dimension(il+2*bfsize,3) :: flux
double precision,intent(out),dimension(3,il+2*bfsize) :: tot
logical,intent(out),dimension(il+2*bfsize) :: flag

double precision :: is0,is1,is2,tau5
double precision :: aa0,aa1,aa2
double precision :: w0,w1,w2
double precision :: tot0
double precision :: lf(-2:2)
integer :: i,m,n

u2=0.0
do n=1,3
	do i=1,il+2*bfsize
		u2(n) = u2(n) + flux(i,n)**2*hx
	enddo
enddo 
do n=1,3
	u2(n) = sqrt(u2(n))
	if(u2(n).lt.1.0e-20) u2(n)=1.0
enddo 

tot = 0.0
do n=1,3
	do i=bfsize,il+bfsize
	        do m=-2,2
	            lf(m)=flux(i+1-m,n)
	        enddo
	        is0= 13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0 &
	            +(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
	        is1= 13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0 &
	        	+(lf(-1)-lf(1))**2/4.d0
	        is2= 13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0 &
	        	+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

	        tot(1,i)=tot(1,i) + is0/u2(n)/3.0
	        tot(2,i)=tot(2,i) + is1/u2(n)/3.0
	        tot(3,i)=tot(3,i) + is2/u2(n)/3.0
	enddo
enddo

do i=bfsize,il+bfsize
	is0 = tot(1,i)
	is1 = tot(2,i)
	is2 = tot(3,i)
	switch_r(i) = is0 + is1 + is2
    if(switch_r(i).gt.1.d0) then
    	flag(i) = .true.
    else
    	flag(i) = .false.
    endif

	tau5 = abs(is2-is1)
    aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
    aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
    aa2=c32*(1.0d0+(tau5/(ss+is2))**2)

	! aa0=c30/(ss+is0)**2
	! aa1=c31/(ss+is1)**2
	! aa2=c32/(ss+is2)**2

	w0=aa0/(aa0+aa1+aa2)
	w1=aa1/(aa0+aa1+aa2)
	w2=aa2/(aa0+aa1+aa2)    

	tot(1,i)=w0
	tot(2,i)=w1
	tot(3,i)=w2
enddo
end subroutine


subroutine scheme_weno5th_tot_1d_character_l(tot,flag,prim,fluxp,fluxl,nl,il,bfsize)
use flow_1d,only:ss,gamma
use weno5_coef
implicit none
double precision,intent(in),dimension(nl,il+2*bfsize) :: prim
double precision,intent(in),dimension(il+2*bfsize,nl) :: fluxp
double precision,intent(in),dimension(3,il+2*bfsize) :: tot
logical,intent(in),dimension(il+2*bfsize) :: flag

double precision,intent(inout),dimension(nl,il+2*bfsize) :: fluxl

double precision :: is0,is1,is2
double precision :: tau5
double precision :: aa0,aa1,aa2
double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision,dimension(-2:2) :: lf
double precision ::char_var(nl)
double precision ::r(3,3),l(3,3)

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,alpha0
double precision :: l1,l2

integer :: i,j,k,n,m,il,nl,bfsize


do i=bfsize,il+bfsize
	if(flag(i)) then
		rhol=prim(1,i)
		ul=prim(2,i)
		pl=prim(3,i)   
		hl=pl/(gamma-1.d0)/rhol+0.5d0*ul*ul+pl/rhol

		rhor=prim(1,i+1)
		ur=prim(2,i+1)
		pr=prim(3,i+1)    
		hr=pr/(gamma-1.d0)/rhor+0.5d0*ur*ur+pr/rhor 

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1.d0)*(h0-0.5d0*u0**2))

		r(1,1)=1.d0
		r(1,2)=1.d0
		r(1,3)=1.d0
		r(2,1)=u0
		r(2,2)=(u0-alpha0)
		r(2,3)=(u0+alpha0)
		r(3,1)=0.5d0*u0**2
		r(3,2)=(h0-u0*alpha0)
		r(3,3)=(h0+u0*alpha0)

		l1=(gamma-1.d0)/alpha0**2
		l2=0.5d0*u0**2

		l(1,1)=1.d0-l2*l1
		l(1,2)=u0*l1
		l(1,3)=-l1
		l(2,1)=0.5d0*(u0+l1*alpha0*l2)/alpha0
		l(2,2)=-0.5d0*(1.d0/alpha0+l1*u0)
		l(2,3)=0.5d0*l1
		l(3,1)=0.5d0*(-u0/alpha0+l1*l2)
		l(3,2)=0.5d0*(1.d0/alpha0-l1*u0)
		l(3,3)=0.5d0*l1

		do j=1,nl
			do m=-2,2        
				lf(m)=(l(j,1)*fluxp(i+m,1)+l(j,2)*fluxp(i+m,2)+l(j,3)*fluxp(i+m,3))                
			enddo
			is0=13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0+(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
			is1=13.d0*(lf(-1)-2.d0*lf(0) +lf(1))**2/12.d0+(lf(-1)-lf(1))**2/4.d0
			is2=13.d0*(lf(0) -2.d0*lf(1) +lf(2))**2/12.d0+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

			q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
			q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
			q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

			! aa0=c30/(ss+is0)**2
			! aa1=c31/(ss+is1)**2
			! aa2=c32/(ss+is2)**2

            tau5 = abs(is2-is1)
            aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
            aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
            aa2=c32*(1.0d0+(tau5/(ss+is2))**2)

			w0=aa0/(aa0+aa1+aa2)
			w1=aa1/(aa0+aa1+aa2)
			w2=aa2/(aa0+aa1+aa2)      

			char_var(j)=w0*q30+w1*q31+w2*q32
		enddo
		do j=1,nl
			fluxl(j,i)=(r(j,1)*char_var(1)+r(j,2)*char_var(2)+r(j,3)*char_var(3)) 
		enddo
	else
        w0=tot(1,i)
        w1=tot(2,i)
        w2=tot(3,i)
		do j=1,nl
	        do m=-2,2
	            lf(m)=fluxp(i+m,j)
	        enddo
	        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
	        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
	        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

	        fluxl(j,i)=w0*q30+w1*q31+w2*q32
	    enddo
	endif
enddo

end subroutine

subroutine scheme_weno5th_tot_1d_character_r(tot,flag,prim,fluxn,fluxr,nl,il,bfsize)
    use flow_1d,only:ss,gamma
use weno5_coef
implicit none
double precision,intent(in),dimension(nl,il+2*bfsize) :: prim
double precision,intent(in),dimension(il+2*bfsize,nl) :: fluxn
double precision,intent(in),dimension(3,il+2*bfsize) :: tot
logical,intent(in),dimension(il+2*bfsize) :: flag

double precision,intent(inout),dimension(nl,il+2*bfsize) :: fluxr

double precision :: is0,is1,is2
double precision :: tau5
double precision :: aa0,aa1,aa2
double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision,dimension(-2:2) :: lf
double precision ::char_var(nl)
double precision ::r(3,3),l(3,3)

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,alpha0
double precision :: l1,l2

integer :: i,j,k,n,m,il,nl,bfsize


do i=bfsize,il+bfsize
	if(flag(i)) then
		rhol=prim(1,i)
		ul=prim(2,i)
		pl=prim(3,i)   
		hl=pl/(gamma-1)/rhol+0.5*ul*ul+pl/rhol

		rhor=prim(1,i+1)
		ur=prim(2,i+1)
		pr=prim(3,i+1)    
		hr=pr/(gamma-1.d0)/rhor+0.5*ur*ur+pr/rhor 

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1.d0)*(h0-0.5d0*u0**2))

		r(1,1)=1.d0
		r(1,2)=1.d0
		r(1,3)=1.d0
		r(2,1)=u0
		r(2,2)=(u0-alpha0)
		r(2,3)=(u0+alpha0)
		r(3,1)=0.5d0*u0**2
		r(3,2)=(h0-u0*alpha0)
		r(3,3)=(h0+u0*alpha0)

		l1=(gamma-1.d0)/alpha0**2
		l2=0.5d0*u0**2

		l(1,1)= 1.d0-l2*l1
		l(1,2)= u0*l1
		l(1,3)=-l1
		l(2,1)= 0.5d0*(u0+l1*alpha0*l2)/alpha0
		l(2,2)=-0.5d0*(1.d0/alpha0+l1*u0)
		l(2,3)= 0.5d0*l1
		l(3,1)= 0.5d0*(-u0/alpha0+l1*l2)
		l(3,2)= 0.5d0*(1.d0/alpha0-l1*u0)
		l(3,3)= 0.5d0*l1

		do j=1,nl
			do m=-2,2                
				lf(m)=(l(j,1)*fluxn(i+1-m,1)+l(j,2)*fluxn(i+1-m,2)+l(j,3)*fluxn(i+1-m,3))                
			enddo	
			is0=13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0+(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
			is1=13.d0*(lf(-1)-2.d0*lf(0) +lf(1))**2/12.d0+(lf(-1)-lf(1))**2/4.d0
			is2=13.d0*(lf(0) -2.d0*lf(1) +lf(2))**2/12.d0+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

			q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
			q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
			q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

            tau5 = abs(is2-is1)
            aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
            aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
            aa2=c32*(1.0d0+(tau5/(ss+is2))**2)
            
			! aa0=c30/(ss+is0)**2
			! aa1=c31/(ss+is1)**2
			! aa2=c32/(ss+is2)**2

			w0=aa0/(aa0+aa1+aa2)
			w1=aa1/(aa0+aa1+aa2)
			w2=aa2/(aa0+aa1+aa2)      

			char_var(j)=w0*q30+w1*q31+w2*q32
		enddo
		do j=1,nl
			fluxr(j,i)=(r(j,1)*char_var(1)+r(j,2)*char_var(2)+r(j,3)*char_var(3)) 
		enddo
	else
        w0=tot(1,i)
        w1=tot(2,i)
        w2=tot(3,i)
		do j=1,nl
	        do m=-2,2
	            lf(m)=fluxn(i+1-m,j)
	        enddo
	        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
	        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
	        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)	        

	        fluxr(j,i)=w0*q30+w1*q31+w2*q32
    	enddo
	endif
enddo
end subroutine
