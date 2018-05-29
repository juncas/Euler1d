    module flow_1d
    implicit none
    integer :: il,jl,nl,bfsize
    double precision,pointer :: q(:,:)
    double precision,pointer :: q0(:,:),q1(:,:),q2(:,:),q3(:,:)
    double precision,pointer :: fluxp(:,:),fluxn(:,:),flux(:,:)
    double precision,pointer :: hp(:,:),hn(:,:)
    double precision,pointer :: h(:,:)
    double precision,pointer :: x(:)    

    double precision,pointer :: theta_l(:,:),theta_r(:,:)
    logical,pointer :: flag_l(:),flag_r(:)
    logical :: if_post
    double precision,pointer :: switch_l(:),switch_r(:)

    double precision :: xl,tl
    double precision :: hx,dt,time
    double precision :: cfl,ss,exmax
    double precision :: u2(3)
    integer :: isave,ischeme_inv,iproblem,iflux_splitting
    integer :: qz,z,it
    integer :: last = 1,irk
    double precision,parameter :: gamma=1.4d0
    double precision,parameter :: pi=4.0d0*atan(1.0d0)
    contains
    subroutine allct_data(il,nl)
    implicit none
    integer,intent(in) :: il,nl
    allocate(q(il,nl))
    allocate(x(il))
    allocate(q0(il,nl),q1(il,nl),q2(il,nl),q3(il,nl))
    allocate(fluxp(il,nl),fluxn(il,nl),flux(il,nl))
    allocate(hp(il,nl),hn(il,nl))
    allocate(h(il,nl))
    allocate(theta_l(3,il))
    allocate(theta_r(3,il))
    allocate(flag_l(il),flag_r(il))
    allocate(switch_l(il),switch_r(il))
    theta_l=1.0
    theta_r=1.0
    flag_l = .false.
    flag_r = .false.
    switch_l = 0.0
    switch_r = 0.0
    hp=0.d0
    hn=0.d0
    h=0.d0
    fluxp=0.d0
    fluxn=0.d0
    flux=0.d0
    q0=0.d0
    q1=0.d0
    q2=0.d0
    q3=0.d0
    end subroutine
    end module flow_1d
