subroutine output_1d
	use flow_1d
	implicit none
	double precision :: rho,u,p,t,theta,count=0.0
	integer :: i,j
	character(len=40) :: filename,temp

	write(temp,'(i8)') it
	select case(iflux_splitting)
	case(1)
		filename='flow_1d_AUSM_'
	case(2)
		filename='flow_1d_CH_'
	case(3)
		filename='flow_1d_SC_'
	case(4)
		filename='flow_1d_LLF_'
	case(5)
		filename='flow_1d_GLF_'
	case(6)
		filename='flow_1d_SW_'
	case(7)
		filename='flow_1d_HLLC_'
	case(8)
		filename='flow_1d_TOT_'
	end select

	select case(ischeme_inv)
	case(1)
		filename=trim(adjustl(filename))//'WENO5.dat'
	case(2)
		filename=trim(adjustl(filename))//'WENO7.dat'
	end select

!-'//trim(adjustl(temp))//'
	if(last.eq.0) then
	open(10,file=filename)
	write(10,*)'title="one dimensional"'
	write(10,*)'variables="x","rho","u","p","t","alpha","theta","switch_l","switch_r"'
	write(10,*)'zone i=',il
	do i=bfsize+1,il+bfsize
		rho=q(i,1)
		u=q(i,2)/rho
		p=(gamma-1.d0)*(q(i,3)-0.5d0*rho*u**2)
		t=p/rho
		if(flag_l(i).or.flag_r(i)) then
			theta=1.0
		else
			theta=0.0
		endif
		write(10,*) x(i),rho,u,p,t,sqrt(gamma*p/rho),theta,switch_l(i),switch_r(i)
	end do 
	close(10)
endif


end subroutine
 