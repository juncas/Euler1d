fc=gfortran -ffree-line-length-200
include= 
lib= 

mod=flow_mod_1d.f90

utility = main_1d.f90\
		  boundary_1d.f90\
		  read_input_1d.f90\
		  output_1d.f90\
		  init_1d.f90\
		  postprocess.f90\
		  post_entropy.f90\
		  tridiagsolve.f90

init    = init_ShuOsher_1d.f90\
          init_Lax_1d.f90\
          init_Sod_1d.f90\
          init_entropy.f90\
          init_blast.f90\
          init_contact_1d.f90

flux    = flux_splitting_1d.f90\
		  flux_splitting_AUSM_1d.f90\
	      flux_splitting_charct_1d.f90\
	      flux_splitting_charct_semi_1d.f90\
	      flux_splitting_HLLC_1d.f90\
	      flux_splitting_LF_1d.f90\
	      flux_splitting_sw_1d.f90\
	      flux_splitting_chrt_prj_1d.f90

scheme  = scheme_1d_general.f90\
		  scheme_weno5th.f90\
		  scheme_weno5th_1d_character.f90\
		  scheme_weno7th.f90\
		  scheme_weno7th_1d_character.f90

TimeAdvance = time_advance_1d.f90

srcs=$(mod) $(utility) $(init) $(flux) $(scheme) $(TimeAdvance)

OBJS=$(srcs:.f90=.o)

default:$(OBJS)
	$(fc) -O2 -o Euler-1d.out $^ $(lib) 
$(OBJS):%.o:%.f90
	$(fc) -O2 -c $(include) $<
clean:
	rm -f *.out *.o *.mod *.dat
