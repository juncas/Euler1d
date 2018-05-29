subroutine flux_splitting_AUSM
use flow_1d
implicit none
double precision,dimension(il+2*bfsize) :: rho,rho_l,rho_r
double precision,dimension(il+2*bfsize) :: u,u_l,u_r
double precision,dimension(il+2*bfsize) :: p,p_l,p_r
double precision,dimension(il+2*bfsize) :: c,c_l,c_r

double precision :: Mach,Mach_l,Mach_r,Mp,Mn,M0
double precision :: Pp,Pn,P0
double precision :: alpha0
integer :: i,n

exmax = 0.0

do i=1,il+2*bfsize   
    rho(i)=q(i,1)
    u(i)=q(i,2)/rho(i)
    p(i)=(gamma-1.d0)*(q(i,3)-0.5d0*rho(i)*(u(i)**2))
    c(i)=sqrt(gamma*p(i)/rho(i))
    exmax=MAX(exmax,abs(u(i))+c(i))   
enddo

    call scheme_WENO5th_l(rho,rho_l,il,bfsize)
    call scheme_WENO5th_l(u,u_l,il,bfsize)
    call scheme_WENO5th_l(p,p_l,il,bfsize)
    call scheme_WENO5th_l(c,c_l,il,bfsize)

    call scheme_WENO5th_r(rho,rho_r,il,bfsize)
    call scheme_WENO5th_r(u,u_r,il,bfsize)
    call scheme_WENO5th_r(p,p_r,il,bfsize)
    call scheme_WENO5th_r(c,c_r,il,bfsize)

do i=bfsize,il+bfsize  
    alpha0=(c_l(i)*rho_l(i)+c_r(i)*rho_r(i))/(rho_l(i)+rho_r(i))

    Mach_l=u_l(i)/alpha0
    Mach_r=u_r(i)/alpha0

    call Get_M_p(Mp,Mach_l)
    call Get_M_n(Mn,Mach_r)
    M0=Mp+Mn

    call Get_P_p(Pp,Mach_l)
    call Get_P_n(Pn,Mach_r)
    P0=Pp*p_l(i)+Pn*p_r(i)

    h(i,1)= alpha0*0.5d0*(M0+abs(M0))*rho_l(i)+alpha0*0.5d0*(M0-abs(M0))*rho_r(i+1)
    h(i,2)= alpha0*0.5d0*(M0+abs(M0))*rho_l(i)*u_l(i)+alpha0*0.5d0*(M0-abs(M0))*rho_r(i+1)*u_r(i+1)+P0
    h(i,3)= alpha0*0.5d0*(M0+abs(M0))*(0.5d0*rho_l(i)*u_l(i)*u_l(i)+p_l(i)/(gamma-1.d0)+p_l(i)) &
           +alpha0*0.5d0*(M0-abs(M0))*(0.5d0*rho_r(i)*u_r(i)*u_r(i)+p_r(i)/(gamma-1.d0)+p_r(i))
enddo

end subroutine

subroutine Get_M_p(Machp,Mach)
implicit none
double precision,intent(in) :: Mach
double precision,intent(inout) :: Machp

if(abs(Mach).le.1.d0) then
    Machp=0.25d0*(Mach+1.d0)**2+0.125d0*(Mach**2-1.d00)**2
else
    Machp=0.5d0*(Mach+abs(Mach))
endif

end subroutine

subroutine Get_M_n(Machn,Mach)
implicit none
double precision,intent(in) :: Mach
double precision,intent(inout) :: Machn

if(abs(Mach).le.1.d0) then
    Machn=-0.25d0*(Mach-1.d0)**2-0.125d0*(Mach**2-1.0d0)**2
else
    Machn=0.5d0*(Mach-abs(Mach))
endif

end subroutine

subroutine Get_P_p(p_p,Mach)
implicit none
double precision,intent(in) :: Mach
double precision,intent(inout) :: p_p

if(abs(Mach).le.1.d0) then
    p_p=0.25d0*(Mach+1.d0)**2*(2.d0-Mach)+0.1875d0*Mach*(Mach**2-1.0d0)**2
else
    p_p=0.5d0*(Mach+abs(Mach))/Mach
endif

end subroutine

subroutine Get_P_n(p_n,Mach)
implicit none
double precision,intent(in) :: Mach
double precision,intent(inout) :: p_n

if(abs(Mach).le.1.d0) then
    p_n=0.25d0*(Mach-1.d0)**2*(2.d0+Mach)-0.1875d0*Mach*(Mach**2-1.0)**2
else
    p_n=0.5d0*(Mach-abs(Mach))/Mach
endif

end subroutine