    subroutine scheme_WENO7th
    use flow_1d
    implicit none
    integer :: i

    call scheme_weno7_p
    call scheme_weno7_n
    h=hp+hn
    end  subroutine

    
    subroutine scheme_weno7_p
	use flow_1d
	implicit none
	double precision :: q40,q41,q42,q43
	double precision :: aa0,aa1,aa2,aa3
	double precision :: w0,w1,w2,w3
	double precision :: is0,is1,is2,is3
	double precision :: tau7
	double precision :: a400,a401,a402,a403
	double precision :: a410,a411,a412,a413
	double precision :: a420,a421,a422,a423
	double precision :: a430,a431,a432,a433
	double precision :: c40,c41,c42,c43
	double precision,dimension(-3:3) :: lf
	integer :: i,j,k,n,m

	a400=-1.d0/4.d0
	a401=13.d0/12.d0
	a402=-23.d0/12.d0
	a403=25.d0/12.d0

	a410=1.d0/12.d0
	a411=-5.d0/12.d0
	a412=13.d0/12.d0
	a413=3.d0/12.d0

	a420=-1.d0/12.d0
	a421=7.d0/12.d0
	a422=7.d0/12.d0
	a423=-1.d0/12.d0

	a430=1.d0/4.d0
	a431=13.d0/12.d0
	a432=-5.d0/12.d0
	a433=1.d0/12.d0

	c40=1.d0/35.d0
	c41=12.d0/35.d0
	c42=18.d0/35.d0
	c43=4.d0/35.d0

	do i=bfsize,il+bfsize
		do j=1,nl
			do m=-3,3                
				lf(m)=fluxp(i+m,j)                
			enddo	
			is0=lf(-3)*(547.d0*lf(-3)-3882.d0*lf(-2)+4642.d0*lf(-1)-1854.d0*lf(0)) &
				+lf(-2)*(7043.d0*lf(-2)-17246.d0*lf(-1)+7042.d0*lf(0)) &
				+lf(-1)*(11003.d0*lf(-1)-9402.d0*lf(0)) &
				+2107.d0*lf(0)**2.d0 

			is1=lf(-2)*(267.d0*lf(-2)-1642.d0*lf(-1)+1602.d0*lf(0)-494.d0*lf(1)) &
				+lf(-1)*(2843.d0*lf(-1)-5966.d0*lf(0)+1922.d0*lf(1)) &
				+lf(0)*(3443.d0*lf(0)-2522.d0*lf(1)) &
				+547.d0*lf(1)**2.d0 


			is2=lf(-1)*(547.d0*lf(-1)-2522.d0*lf(0)+1922.d0*lf(1)-494.d0*lf(2)) &
				+lf(0)*(3443.d0*lf(0)-5966.d0*lf(1)+1602.d0*lf(2)) &
				+lf(1)*(2843.d0*lf(1)-1642.d0*lf(2)) &
				+267.d0*lf(2)**2.d0 

			is3=lf(0)*(2107.d0*lf(0)-9402.d0*lf(1)+7042.d0*lf(2)-1854.d0*lf(3)) &
				+lf(1)*(11003.d0*lf(1)-17246.d0*lf(2)+4642.d0*lf(3)) &
				+lf(2)*(7043.d0*lf(2)-3882.d0*lf(3)) &
				+547.d0*lf(3)**2.d0 

			q40=a400*lf(-3)+a401*lf(-2)+a402*lf(-1)+a403*lf(0)  
			q41=a410*lf(-2)+a411*lf(-1)+a412*lf(0)+a413*lf(1)
			q42=a420*lf(-1)+a421*lf(0)+a422*lf(1)+a423*lf(2)
			q43=a430*lf(0)+a431*lf(1)+a432*lf(2)+a433*lf(3)

			tau7=abs(is0-is3)

			aa0=c40*(1.d0+(tau7/(ss+is0))**2)
			aa1=c41*(1.d0+(tau7/(ss+is1))**2)
			aa2=c42*(1.d0+(tau7/(ss+is2))**2)
			aa3=c43*(1.d0+(tau7/(ss+is3))**2)

			w0=aa0/(aa0+aa1+aa2+aa3)
			w1=aa1/(aa0+aa1+aa2+aa3)
			w2=aa2/(aa0+aa1+aa2+aa3)
			w3=aa3/(aa0+aa1+aa2+aa3)

			hp(i,j)=w0*q40+w1*q41+w2*q42+w3*q43
		enddo
	enddo
end subroutine

subroutine scheme_weno7_n
	use flow_1d
	implicit none
	double precision :: q40,q41,q42,q43
	double precision :: aa0,aa1,aa2,aa3
	double precision :: w0,w1,w2,w3
	double precision :: is0,is1,is2,is3
	double precision :: tau7
	double precision :: a400,a401,a402,a403
	double precision :: a410,a411,a412,a413
	double precision :: a420,a421,a422,a423
	double precision :: a430,a431,a432,a433
	double precision :: c40,c41,c42,c43
	double precision,dimension(-3:3) :: lf
	integer :: i,j,k,n,m

	a400=-1.d0/4.d0
	a401=13.d0/12.d0
	a402=-23.d0/12.d0
	a403=25.d0/12.d0

	a410=1.d0/12.d0
	a411=-5.d0/12.d0
	a412=13.d0/12.d0
	a413=3.d0/12.d0

	a420=-1.d0/12.d0
	a421=7.d0/12.d0
	a422=7.d0/12.d0
	a423=-1.d0/12.d0

	a430=1.d0/4.d0
	a431=13.d0/12.d0
	a432=-5.d0/12.d0
	a433=1.d0/12.d0

	c40=1.d0/35.d0
	c41=12.d0/35.d0
	c42=18.d0/35.d0
	c43=4.d0/35.d0

	do i=bfsize,il+bfsize
		do j=1,nl
			do m=-3,3                
				lf(m)=fluxn(i+1-m,j)    
			enddo	
			is0=lf(-3)*(547.d0*lf(-3)-3882.d0*lf(-2)+4642.d0*lf(-1)-1854.d0*lf(0)) &
				+lf(-2)*(7043.d0*lf(-2)-17246.d0*lf(-1)+7042.d0*lf(0)) &
				+lf(-1)*(11003.d0*lf(-1)-9402.d0*lf(0)) &
				+2107.d0*lf(0)**2.d0 

			is1=lf(-2)*(267.d0*lf(-2)-1642.d0*lf(-1)+1602.d0*lf(0)-494.d0*lf(1)) &
				+lf(-1)*(2843.d0*lf(-1)-5966.d0*lf(0)+1922.d0*lf(1)) &
				+lf(0)*(3443.d0*lf(0)-2522.d0*lf(1)) &
				+547.d0*lf(1)**2.d0 


			is2=lf(-1)*(547.d0*lf(-1)-2522.d0*lf(0)+1922.d0*lf(1)-494.d0*lf(2)) &
				+lf(0)*(3443.d0*lf(0)-5966.d0*lf(1)+1602.d0*lf(2)) &
				+lf(1)*(2843.d0*lf(1)-1642.d0*lf(2)) &
				+267.d0*lf(2)**2.d0 

			is3=lf(0)*(2107.d0*lf(0)-9402.d0*lf(1)+7042.d0*lf(2)-1854.d0*lf(3)) &
				+lf(1)*(11003.d0*lf(1)-17246.d0*lf(2)+4642.d0*lf(3)) &
				+lf(2)*(7043.d0*lf(2)-3882.d0*lf(3)) &
				+547.d0*lf(3)**2.d0 

			q40=a400*lf(-3)+a401*lf(-2)+a402*lf(-1)+a403*lf(0)  
			q41=a410*lf(-2)+a411*lf(-1)+a412*lf(0)+a413*lf(1)
			q42=a420*lf(-1)+a421*lf(0)+a422*lf(1)+a423*lf(2)
			q43=a430*lf(0)+a431*lf(1)+a432*lf(2)+a433*lf(3)

			tau7=abs(is0-is3)

			aa0=c40*(1.d0+(tau7/(ss+is0))**2)
			aa1=c41*(1.d0+(tau7/(ss+is1))**2)
			aa2=c42*(1.d0+(tau7/(ss+is2))**2)
			aa3=c43*(1.d0+(tau7/(ss+is3))**2)

			w0=aa0/(aa0+aa1+aa2+aa3)
			w1=aa1/(aa0+aa1+aa2+aa3)
			w2=aa2/(aa0+aa1+aa2+aa3)
			w3=aa3/(aa0+aa1+aa2+aa3)

			hn(i,j)=w0*q40+w1*q41+w2*q42+w3*q43
		enddo
	enddo
end subroutine
