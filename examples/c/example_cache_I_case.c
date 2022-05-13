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
  int varid_spinup_state;
  int var_val;
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

  char filename[] = "I1850GSWCNPRDCTCBC_f19_g16_ADIOS.elm.r.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename, PIO_NOWRITE); ERR

  varid_spinup_state = -1;
  ret = PIOc_inq_varid(ncid_read, "spinup_state", &varid_spinup_state); ERR

  var_val = -1;
  ret = PIOc_get_var_int(ncid_read, varid_spinup_state, &var_val); ERR
  if (my_rank == 0)
    printf("1st PIOc_get_var_int call, var_val = %d\n", var_val);

  var_val = -1;
  ret = PIOc_get_var_int(ncid_read, varid_spinup_state, &var_val); ERR
  if (my_rank == 0)
    printf("2nd PIOc_get_var_int call, var_val = %d\n", var_val);

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}

