    subroutine flux_splitting_1d
    use flow_1d
    ! print*,"splitting Method=",iflux_splitting
    select case (iflux_splitting)
    case(1)
        call flux_splitting_AUSM
    case(2)
        call flux_splitting_chrt_1d
    case(3)
        call flux_splitting_chrt_semi_1d
    case(4)
        call flux_splitting_LLF_1d  
    case(5)
        call flux_splitting_GLF_1d
    case(6)
        call flux_splitting_sw_1d
    case(7)
        call flux_splitting_HLLC
    case(8)
        call flux_splitting_chrt_prj_1d
    end select

    end subroutine