   subroutine scheme_WENO5th
    use flow_1d
    implicit none
    integer :: i

    call scheme_WENO5th_p
    call scheme_WENO5th_n
    h=hp+hn
    end  subroutine

    subroutine scheme_WENO5th_p
    use flow_1d
    implicit none
    double precision :: is0,is1,is2
    double precision :: w0,w1,w2
    double precision :: tau5
    double precision :: aa0,aa1,aa2
    double precision :: a300,a301,a302,a310,a311,a312,a320,a321,a322
    double precision :: c30,c31,c32
    double precision :: q30,q31,q32
    integer :: i,n,j,m
    double precision :: lf(-2:2)
    a300=1.d0/3.d0
    a301=-7.d0/6.d0
    a302=11.d0/6.d0

    a310=-1.d0/6.d0
    a311=5.d0/6.d0
    a312=1.d0/3.d0

    a320=1.d0/3.d0
    a321=5.d0/6.d0
    a322=-1.d0/6.d0

    c30=0.1d0
    c31=0.6d0
    c32=0.3d0
    do i=bfsize,il+bfsize
        do j=1,nl
            do m=-2,2
                lf(m)=fluxp(i+m,j)
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

            w0=aa0/(aa0+aa1+aa2)
            w1=aa1/(aa0+aa1+aa2)
            w2=aa2/(aa0+aa1+aa2)
            hp(i,j)=w0*q30+w1*q31+w2*q32
        enddo
    enddo
    end subroutine

    subroutine scheme_WENO5th_n
    use flow_1d
    implicit none
    double precision :: is0,is1,is2
    double precision :: w0,w1,w2
    double precision :: tau5
    double precision :: aa0,aa1,aa2
    double precision :: a300,a301,a302,a310,a311,a312,a320,a321,a322
    double precision :: c30,c31,c32
    double precision :: q30,q31,q32
    integer :: i,n,j,m
    double precision :: lf(-2:2)
    a300=1.d0/3.d0
    a301=-7.d0/6.d0
    a302=11.d0/6.d0

    a310=-1.d0/6.d0
    a311=5.d0/6.d0
    a312=1.d0/3.d0

    a320=1.d0/3.d0
    a321=5.d0/6.d0
    a322=-1.d0/6.d0

    c30=0.1d0
    c31=0.6d0
    c32=0.3d0
    do i=bfsize,il+bfsize
        do j=1,nl
            do m=-2,2
                lf(m)=fluxn(i+1-m,j)
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

            w0=aa0/(aa0+aa1+aa2)
            w1=aa1/(aa0+aa1+aa2)
            w2=aa2/(aa0+aa1+aa2)
            hn(i,j)=w0*q30+w1*q31+w2*q32
        enddo
    enddo
    end subroutine


    subroutine scheme_WENO5th_l(prim,recon_prim,il,bfsize)
    use flow_1d,only:ss
    implicit none

    integer,intent(in) :: il,bfsize
    double precision,intent(in),dimension(il+2*bfsize) :: prim
    double precision,intent(out),dimension(il+2*bfsize) :: recon_prim

    double precision :: is0,is1,is2
    double precision :: w0,w1,w2
    double precision :: tau5
    double precision :: aa0,aa1,aa2
    double precision :: a300,a301,a302,a310,a311,a312,a320,a321,a322
    double precision :: c30,c31,c32
    double precision :: q30,q31,q32
    integer :: i,n,j,m
    double precision :: lf(-2:2)

    a300=1.d0/3.d0
    a301=-7.d0/6.d0
    a302=11.d0/6.d0

    a310=-1.d0/6.d0
    a311=5.d0/6.d0
    a312=1.d0/3.d0

    a320=1.d0/3.d0
    a321=5.d0/6.d0
    a322=-1.d0/6.d0

    c30=0.1d0
    c31=0.6d0
    c32=0.3d0
    do i=bfsize,il+bfsize
            do m=-2,2
                lf(m)=prim(i+m)
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

            w0=aa0/(aa0+aa1+aa2)
            w1=aa1/(aa0+aa1+aa2)
            w2=aa2/(aa0+aa1+aa2)
            recon_prim(i)=w0*q30+w1*q31+w2*q32
    enddo
    end subroutine

    subroutine scheme_WENO5th_r(prim,recon_prim,il,bfsize)
    use flow_1d,only:ss
    implicit none

    integer,intent(in) :: il,bfsize
    double precision,intent(in),dimension(il+2*bfsize) :: prim
    double precision,intent(out),dimension(il+2*bfsize) :: recon_prim

    double precision :: is0,is1,is2
    double precision :: w0,w1,w2
    double precision :: tau5
    double precision :: aa0,aa1,aa2
    double precision :: a300,a301,a302,a310,a311,a312,a320,a321,a322
    double precision :: c30,c31,c32
    double precision :: q30,q31,q32
    integer :: i,n,j,m
    double precision :: lf(-2:2)
    a300=1.d0/3.d0
    a301=-7.d0/6.d0
    a302=11.d0/6.d0

    a310=-1.d0/6.d0
    a311=5.d0/6.d0
    a312=1.d0/3.d0

    a320=1.d0/3.d0
    a321=5.d0/6.d0
    a322=-1.d0/6.d0

    c30=0.1d0
    c31=0.6d0
    c32=0.3d0
    do i=bfsize,il+bfsize
            do m=-2,2
                lf(m)=prim(i+1-m)
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

            w0=aa0/(aa0+aa1+aa2)
            w1=aa1/(aa0+aa1+aa2)
            w2=aa2/(aa0+aa1+aa2)
            recon_prim(i)=w0*q30+w1*q31+w2*q32
    enddo
    end subroutine

