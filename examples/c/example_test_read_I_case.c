#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

int main(int argc, char* argv[])
{
  int my_rank;
  int ntasks;
  int format = PIO_IOTYPE_ADIOS;
  int niotasks;
  const int ioproc_start = 0;
  int iosysid;
  int ncid1;
  int ncid2;
  int varid;
  double buffer[30];
  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  int ioproc_stride = 1;
  niotasks = ntasks / ioproc_stride;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR

  char filename1[] = "I1850GSWCNPRDCTCBC_f19_g16_ADIOS.elm.r.0001-01-02-00000.nc";
  ret = PIOc_openfile(iosysid, &ncid1, &format, filename1, PIO_NOWRITE); ERR

  char filename2[] = "I1850GSWCNPRDCTCBC_f19_g16_ADIOS.elm.rh0.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid2, &format, filename2, PIO_NOWRITE); ERR


  varid = 151; /* double budg_fluxG(budg_flux = 30); */
  ret = PIOc_get_var_double(ncid1, varid, buffer); ERR

  ret = PIOc_closefile(ncid2); ERR
  ret = PIOc_closefile(ncid1); ERR


  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}

