#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#include <assert.h>
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
  int ncid_read;
  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  int ioproc_stride = 2;
  if (ntasks == 1)
    ioproc_stride = 1;

  niotasks = ntasks / ioproc_stride;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR

  char filename_eam_r[] = "F2010_ne4_oQU240_ADIOS.eam.r.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_r, PIO_NOWRITE); ERR

  /* There is no global attribute named NLON in the test file so the inquiry below should return PIO_ENOTATT */
  char att_val[PIO_MAX_NAME];
  ret = PIOc_get_att(ncid_read, PIO_GLOBAL, "NLON", att_val);
  assert(ret == PIO_ENOTATT);

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
