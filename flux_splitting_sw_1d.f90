subroutine flux_splitting_sw_1d
use flow_1d
implicit none    
double precision :: rho,u,v,p,c,a0,wii,epsl
double precision,dimension(3,2) :: lambda
integer :: i,j,n
epsl=0.125
do i=1,il+2*bfsize !LF splitting  
    rho=q(i,1)
    u=q(i,2)/rho
    p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*(u**2))
    c=sqrt(gamma*p/rho)
    a0=rho/(2.d0*gamma)
    !-----------------------------x direction---------------------------------------------------
    lambda(1,1)=(u+abs(sqrt(u**2+epsl)))/2.d0
    lambda(2,1)=(u-c+abs(sqrt((u-c)**2+epsl)))/2.d0
    lambda(3,1)=(u+c+abs(sqrt((u+c)**2+epsl)))/2.d0
    wii=(3.d0-gamma)*(lambda(2,1)+lambda(3,1))*c**2/(2.d0*(gamma-1.d0))
    fluxp(i,1)=2.d0*(gamma-1.d0)*lambda(1,1)+lambda(2,1)+lambda(3,1)
    fluxp(i,2)=2.d0*(gamma-1.d0)*lambda(1,1)*u+lambda(2,1)*(u-c)+lambda(3,1)*(u+c)
    fluxp(i,3)=(gamma-1.d0)*lambda(1,1)*(u**2)+0.5d0*lambda(2,1)*((u-c)**2)+0.5d0*lambda(3,1)*((u+c)**2)+wii

    lambda(1,1)=(u-abs(sqrt(u**2+epsl)))/2.d0
    lambda(2,1)=(u-c-abs(sqrt((u-c)**2+epsl)))/2.d0
    lambda(3,1)=(u+c-abs(sqrt((u+c)**2+epsl)))/2.d0
    wii=(3.d0-gamma)*(lambda(2,1)+lambda(3,1))*c**2/(2.d0*(gamma-1.d0))
    fluxn(i,1)=2.d0*(gamma-1.d0)*lambda(1,1)+lambda(2,1)+lambda(3,1)
    fluxn(i,2)=2.d0*(gamma-1.d0)*lambda(1,1)*u+lambda(2,1)*(u-c)+lambda(3,1)*(u+c)
    fluxn(i,3)=(gamma-1.d0)*lambda(1,1)*(u**2)+0.5d0*lambda(2,1)*((u-c)**2)+0.5d0*lambda(3,1)*((u+c)**2)+wii
    exmax=MAX(exmax,abs(u-c),abs(u+c))
    do n=1,nl
        fluxn(i,n)=fluxn(i,n)*a0
        fluxp(i,n)=fluxp(i,n)*a0               
    enddo        
enddo
call scheme_invis_1d
end subroutine