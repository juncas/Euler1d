subroutine flux_splitting_HLLC
use flow_1d
implicit none
double precision,dimension(il+2*bfsize) :: rho,rho_l,rho_r
double precision,dimension(il+2*bfsize) :: u,u_l,u_r
double precision,dimension(il+2*bfsize) :: p,p_l,p_r
double precision,dimension(il+2*bfsize) :: c

double precision :: S_l,S_r,S_star,c_l,c_r
double precision :: u_bar,d_bar,eta
double precision :: U_star(3)
double precision :: p_star

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

call scheme_WENO5th_r(rho,rho_r,il,bfsize)
call scheme_WENO5th_r(u,u_r,il,bfsize)
call scheme_WENO5th_r(p,p_r,il,bfsize)

do i=bfsize,il+bfsize  
    c_l=sqrt(gamma*p_l(i)/rho_l(i))
    c_r=sqrt(gamma*p_r(i)/rho_r(i))

    eta = 0.5*sqrt(rho_l(i))*sqrt(rho_r(i))/(sqrt(rho_l(i))+sqrt(rho_r(i)))**2
    u_bar = (sqrt(rho_l(i))*u_l(i)+sqrt(rho_r(i))*u_r(i))/(sqrt(rho_l(i))+sqrt(rho_r(i)))
    d_bar = sqrt((sqrt(rho_l(i))*c_l**2+sqrt(rho_r(i))*c_r**2)/(sqrt(rho_l(i))+sqrt(rho_r(i)))+eta*(u_l(i)-u_r(i))**2)

    S_l = u_bar-d_bar
    S_r = u_bar+d_bar
    S_star = (p_l(i)-p_r(i)-rho_l(i)*u_l(i)*(S_l-u_l(i))+rho_r(i)*u_r(i)*(S_r-u_r(i)))/(rho_l(i)*u_l(i)-rho_r(i)*u_r(i)+S_r*rho_r(i)-S_l*rho_l(i))


    if (0.0.le.S_l) then
        h(i,1)= u_l(i)*rho_l(i)
        h(i,2)= u_l(i)*rho_l(i)*u_l(i)+p_l(i)
        h(i,3)= u_l(i)*(0.5*rho_l(i)*u_l(i)*u_l(i)+p_l(i)/(gamma-1))+u_l(i)*p_l(i)
    elseif((S_l.lt.0.0).and.(0.0.le.S_star)) then    
        U_star(1)=rho_l(i)*(S_l-u_l(i))/(S_l-S_star)
        U_star(2)=rho_l(i)*(S_l-u_l(i))/(S_l-S_star)*S_star
        U_star(3)=rho_l(i)*(S_l-u_l(i))/(S_l-S_star)*(0.5*u_l(i)*u_l(i)+p_l(i)/(gamma-1)/rho_l(i)+(S_star-u_l(i))*(S_star+p_l(i)/rho_l(i)/(S_l-u_l(i))))         
        p_star = p_l(i)+rho_l(i)*(S_l-u_l(i))*(S_star-u_l(i))

        h(i,1)= S_star*U_star(1)
        h(i,2)= S_star*U_star(2)+p_star
        h(i,3)= S_star*U_star(3)+S_star*p_star
    elseif((S_star.lt.0.0).and.(0.0.le.S_r)) then  
        U_star(1)=rho_r(i)*(S_r-u_r(i))/(S_r-S_star)
        U_star(2)=rho_r(i)*(S_r-u_r(i))/(S_r-S_star)*S_star
        U_star(3)=rho_r(i)*(S_r-u_r(i))/(S_r-S_star)*(0.5*u_r(i)*u_r(i)+p_r(i)/(gamma-1)/rho_r(i)+(S_star-u_r(i))*(S_star+p_r(i)/rho_r(i)/(S_r-u_r(i))))         
        p_star = p_r(i)+rho_r(i)*(S_r-u_r(i))*(S_star-u_r(i))

        h(i,1)= S_star*U_star(1)
        h(i,2)= S_star*U_star(2)+p_star
        h(i,3)= S_star*U_star(3)+S_star*p_star
    elseif((0.0.gt.S_r)) then                
        h(i,1)= u_r(i)*rho_r(i)
        h(i,2)= u_r(i)*rho_r(i)*u_r(i)+p_r(i)
        h(i,3)= u_r(i)*(0.5*rho_r(i)*u_r(i)*u_r(i)+p_r(i)/(gamma-1))+u_r(i)*p_r(i)
    endif    
enddo

end subroutine
