#include <stdio.h>
#include <mpi.h>
#include <pio.h>
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
  int varid_ANRAIN;
  PIO_Offset att_len = 0;
  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  if (ntasks != 16) {
    if (my_rank == 0)
      printf("This example must be run with 16 MPI tasks!\n");

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return -1;
  }

  int ioproc_stride = 4;
  niotasks = ntasks / ioproc_stride;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR

  char filename_eam_h0[] = "F2010_ne4_oQU240_ADIOS.eam.h0.0001-01-01-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_h0, PIO_NOWRITE); ERR

  varid_ANRAIN = -1;
  ret = PIOc_inq_varid(ncid_read, "ANRAIN", &varid_ANRAIN); ERR

  att_len = 0;
  ret = PIOc_inq_attlen(ncid_read, varid_ANRAIN, "mdims", &att_len); ERR
  printf("rank = %d, att_len = %lld\n", my_rank, att_len);

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
