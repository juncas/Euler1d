    subroutine init_1d
    use flow_1d
    select case(iproblem)
    case(1)
        call init_Lax_1d
    case(2)
        call init_ShuOsher_1d
    case(3)
        call init_Sod_1d
    case(4)
        call init_blast_1d
    case(5)
        call init_contact_1d
    case(6)
        call init_entropy_1d
    end select
    end subroutine