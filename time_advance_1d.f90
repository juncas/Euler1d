    subroutine time_rk3_1d
    use flow_1d
    implicit none
    integer :: i,n
    double precision ::dtxi,dteta
    select case (irk)
    case(1)
        if(cfl.lt.1.0) then
            dt=cfl*hx/exmax
        else 
            dt=0.05d0*hx**2
        endif
        q0=q        
        do n=1,nl
            do i=bfsize+1,il+bfsize                
                q1(i,n)=q0(i,n) &
                -dt/hx*(h(i,n)-h(i-1,n))             
                q(i,n)=q1(i,n)
            end do
        end do
    case(2)
        do n=1,nl
            do i=bfsize+1,il+bfsize               
                q2(i,n)=3.d0/4.d0*q0(i,n)+1.d0/4.d0*q1(i,n) &
                -dt/hx*(h(i,n)-h(i-1,n))/4.d0                         
                q(i,n)=q2(i,n)                
            end do
        end do
    case(3)
        do n=1,nl
            do i=bfsize+1,il+bfsize                
                q3(i,n)=1.d0/3.d0*q0(i,n)+2.d0/3.d0*q2(i,n) &
                -dt/hx*(h(i,n)-h(i-1,n))*2.d0/3.d0                        
                q(i,n)=q3(i,n)               
            end do
        end do
    end select
    end subroutine
