subroutine postprocess
	use flow_1d
	implicit none

    select case(iproblem)
    case(1)
        !call init_Lax_1d
    case(2)
        !call init_ShuOsher_1d
    case(3)
        !call init_Sod_1d
    case(4)
        !call init_blast_1d
    case(5)
        !call init_contact_1d
    case(6)
        call post_entropy
    end select

end

subroutine postprocess_fly
    use flow_1d
    implicit none
    double precision :: count_ch,count_tot
    integer :: i

    if(iflux_splitting.eq.3) then
        open(11,file="SC_percentage.dat",ACCESS='APPEND')
        count_ch=0.0
        do i=bfsize+1,il+bfsize
            if(flag_l(i).or.flag_r(i)) then
                count_ch = count_ch+1.0
            endif
        end do 
        write(11,*) time,count_ch/il*100.d0
        close(11)
    endif

    if(iflux_splitting.eq.8) then
        open(11,file="TOT_percentage.dat",ACCESS='APPEND')
        count_tot=0.0
        do i=bfsize+1,il+bfsize
            if(flag_l(i).or.flag_r(i)) then
                count_tot = count_tot+1.0
            endif
        end do 
        write(11,*) time,count_tot/il*100.d0
        close(11)
    endif

end