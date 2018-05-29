subroutine flux_splitting_LLF_1d
use flow_1d
implicit none    
double precision :: rho,u,v,p,c
integer :: i,n

exmax = 0.0
do i=1,il+2*bfsize !LF splitting
    rho=q(i,1)
    u=q(i,2)/rho
    p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*(u**2))
    c=sqrt(gamma*p/rho)
    
    exmax=MAX(exmax,abs(u)+c)   
    flux(i,1)=q(i,1)*u
    flux(i,2)=q(i,2)*u+p
    flux(i,3)=q(i,3)*u+u*p

    do n=1,nl
        fluxp(i,n)=0.5d0*(flux(i,n)+(abs(u)+c)*q(i,n))
        fluxn(i,n)=0.5d0*(flux(i,n)-(abs(u)+c)*q(i,n))                
    enddo        
enddo
    call scheme_invis_1d
end subroutine

subroutine flux_splitting_GLF_1d
use flow_1d
implicit none    
double precision :: rho,u,v,p,c,a0,wii,epsl
double precision,dimension(3,2) :: lambda
integer :: i,j,n

exmax = 0.0
do i=1,il+2*bfsize !LF splitting
    rho=q(i,1)
    u=q(i,2)/rho
    p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*(u**2))
    c=sqrt(gamma*p/rho)

    exmax=MAX(exmax,abs(u)+c)    
    flux(i,1)=q(i,1)*u
    flux(i,2)=q(i,2)*u+p
    flux(i,3)=q(i,3)*u+u*p 
enddo

do i=1,il+2*bfsize !LF splitting
    do n=1,nl
        fluxp(i,n)=0.5d0*(flux(i,n)+exmax*q(i,n))
        fluxn(i,n)=0.5d0*(flux(i,n)-exmax*q(i,n))        
    enddo        
enddo
    call scheme_invis_1d
end subroutine