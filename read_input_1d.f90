subroutine read_input_1d
	use flow_1d
	implicit none

	open(10,file='input.in',form='formatted')
	write(*,*) 'Reading...'
	read(10,*)
	read(10,*) il,nl,bfsize
	write(*,*) 'Xn=',il,'buff size=',bfsize
	write(*,*) 'Dimension=',nl
	read(10,*)
	read(10,*) cfl,ss,qz,z,if_post
	read(10,*)
	read(10,*) xl,tl
	read(10,*)
	read(10,*) isave,dt
	read(10,*)
	read(10,*) iflux_splitting
	read(10,*)
	read(10,*) ischeme_inv
	read(10,*)
	read(10,*) iproblem
endsubroutine

