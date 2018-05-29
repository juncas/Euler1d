subroutine boundary_1d
	use flow_1d
    select case(iproblem)
    case(1)
        !call bndy_Lax_1d
    case(2)
        !call bndy_ShuOsher_1d
    case(3)
        !call bndy_Sod_1d
    case(4)
        call bndy_blast_1d
    case(5)
        !call bndy_contact_1d
    case(6)
        call bndy_entropy_1d
    end select
end subroutine
