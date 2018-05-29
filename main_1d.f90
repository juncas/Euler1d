program Euler1d
use flow_1d
double precision :: dtx,dty
double precision :: cput0,cput1
call read_input_1d
call allct_data(il+2*bfsize,nl)
call init_1d


it=1
time=0.d0
call CPU_TIME(cput0)
call postprocess_fly
do while(time.le.tl)
    do irk=1,3
        call boundary_1d
        call flux_splitting_1d
        call time_rk3_1d
    enddo
    write(*,*) it,'time= ',time,'dt=', dt 
    if(MOD(it,isave).eq.0) then
        print*,"Saving..."
        print*,"Done"
        call output_1d
    endif
    it=it+1
    time=time+dt
    if(if_post) call postprocess_fly
enddo 
last = 0
call CPU_TIME(cput1)
print*,'Finished, cpu time is ',cput1-cput0

open(10,file="cpuTime.dat",ACCESS='APPEND')
write(10,*) il,cput1-cput0
close(10)
call output_1d
call postprocess
end program
