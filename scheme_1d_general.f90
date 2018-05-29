subroutine scheme_invis_1d
	use flow_1d
	select case (ischeme_inv)
	case(1)
		call scheme_WENO5th
	case(2)
		call scheme_WENO7th
	end select
end subroutine
