subroutine scheme_weno5th_1d_character_l(l,r,fluxp,fluxl,nl,il,bfsize)
use flow_1d,only:ss
use weno5_coef
implicit none
double precision,intent(in),dimension(nl,nl,il+2*bfsize) :: l,r
double precision,intent(in),dimension(il+2*bfsize,nl) :: fluxp

double precision,intent(inout),dimension(nl,il+2*bfsize) :: fluxl

double precision :: is0,is1,is2
double precision :: w0,w1,w2
double precision :: w00,w01,w10,w11
double precision :: tau5,tau4,tau3
double precision :: aa0,aa1,aa2
double precision :: a0,a1,aa00,aa01,aa10,aa11
double precision :: fi0,fi1,fi00,fi01,fi10,fi11
double precision :: g0,g1,g00,g01,g10,g11
double precision :: c0,c1
double precision :: q30,q31,q32,temp1(nl),temp2(nl)
double precision,dimension(-2:2) :: lf
integer :: i,j,k,n,m,il,nl,bfsize


do i=bfsize,il+bfsize
	do j=1,nl
		do m=-2,2        
			lf(m)=(l(j,1,i)*fluxp(i+m,1)+l(j,2,i)*fluxp(i+m,2)+l(j,3,i)*fluxp(i+m,3))                
		enddo
		is0=13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0+(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
		is1=13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0+(lf(-1)-lf(1))**2/4.d0
		is2=13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

		q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
		q31=a310*lf(-1)+a311*lf(0)  +a312*lf(1)
		q32=a320*lf(0)  +a321*lf(1)+a322*lf(2)

        tau5 = abs(is2-is1)
        aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
        aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
        aa2=c32*(1.0d0+(tau5/(ss+is2))**2)

			! aa0=c30/(ss+is0)**2
			! aa1=c31/(ss+is1)**2
			! aa2=c32/(ss+is2)**2

		w0=aa0/(aa0+aa1+aa2)
		w1=aa1/(aa0+aa1+aa2)
		w2=aa2/(aa0+aa1+aa2)
		temp1(j)=w0*q30+w1*q31+w2*q32
	enddo
	do j=1,nl
		fluxl(j,i)=(r(j,1,i)*temp1(1)+r(j,2,i)*temp1(2)+r(j,3,i)*temp1(3)) 
	enddo
enddo


end subroutine

subroutine scheme_weno5th_1d_character_r(l,r,fluxn,fluxr,nl,il,bfsize)
use flow_1d,only:ss
use weno5_coef
implicit none
	double precision,intent(in),dimension(nl,nl,il+2*bfsize) :: l,r
	double precision,intent(in),dimension(il+2*bfsize,nl) :: fluxn

	double precision,intent(inout),dimension(nl,il+2*bfsize) :: fluxr

	double precision :: is0,is1,is2
	double precision :: w0,w1,w2
	double precision :: w00,w01,w10,w11
	double precision :: tau5,tau4,tau3
	double precision :: aa0,aa1,aa2
	double precision :: a0,a1,aa00,aa01,aa10,aa11
	double precision :: fi0,fi1,fi00,fi01,fi10,fi11
	double precision :: g0,g1,g00,g01,g10,g11
	double precision :: c0,c1
	double precision :: q30,q31,q32,temp1(nl),temp2(nl)
	double precision,dimension(-2:2) :: lf
	integer :: i,j,k,n,m,il,nl,bfsize

	do i=bfsize,il+bfsize
		do j=1,nl
			do m=-2,2                
				lf(m)=(l(j,1,i)*fluxn(i+1-m,1)+l(j,2,i)*fluxn(i+1-m,2)+l(j,3,i)*fluxn(i+1-m,3))                
			enddo	
			is0=13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0+(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0
			is1=13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0+(lf(-1)-lf(1))**2/4.d0
			is2=13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

			q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
			q31=a310*lf(-1)+a311*lf(0)  +a312*lf(1)
			q32=a320*lf(0)  +a321*lf(1)+a322*lf(2)

            tau5 = abs(is2-is1)
            aa0=c30*(1.0d0+(tau5/(ss+is0))**2)
            aa1=c31*(1.0d0+(tau5/(ss+is1))**2)
            aa2=c32*(1.0d0+(tau5/(ss+is2))**2)


			! aa0=c30/(ss+is0)**2
			! aa1=c31/(ss+is1)**2
			! aa2=c32/(ss+is2)**2

			w0=aa0/(aa0+aa1+aa2)
			w1=aa1/(aa0+aa1+aa2)
			w2=aa2/(aa0+aa1+aa2)            
			temp1(j)=w0*q30+w1*q31+w2*q32
		enddo
		do j=1,nl
			fluxr(j,i)=(r(j,1,i)*temp1(1)+r(j,2,i)*temp1(2)+r(j,3,i)*temp1(3)) 
		enddo
	enddo
end subroutine
