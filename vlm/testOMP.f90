PROGRAM testOMP

  USE OMP_LIB

  IMPLICIT NONE

  WRITE(*,*) "Num procs", omp_get_num_procs()

END PROGRAM testOMP
