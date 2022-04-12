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
  /* Zero-based rank of processor. */
  int my_rank;

  /* Number of processors involved in current execution. */
  int ntasks;

  /* This example only tests two formats (IO types): PnetCDF and ADIOS */
  int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};

  /* Number of processors that will do IO. In this example we
     will do IO from all processors. */
  int niotasks;

  /* Stride in the MPI rank between IO tasks. Always 1 in this example. */
  const int ioproc_stride = 1;

  /* Zero based rank of first processor to be used for I/O. */
  const int ioproc_start = 0;

  /* The ID for the parallel I/O system. It is set by
     PIOc_Init_Intracomm(). It references an internal structure
     containing the general IO subsystem data and MPI
     structure. It is passed to PIOc_finalize() to free
     associated resources, after all I/O, but before
     MPI_Finalize is called. */
  int iosysid;

  /* The ncid of the netCDF file created and read in this example. */
  int ncid_write;
  int ncid_read;

  /* The IDs of the netCDF variables in the example file. */
  int varid_dummy_scalar_var_int_write = -1;
  int varid_dummy_scalar_var_int_read = -1;

  /* Test filename. */
  char filename[PIO_MAX_NAME];

  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /* Keep things simple - 1 IO task per MPI process */
  niotasks = ntasks;

  /* Initialize the PIO IO system. This specifies how
     many and which processors are involved in I/O. */
  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                      PIO_REARR_SUBSET, &iosysid); ERR

  for (int fmt = 0; fmt < 2; fmt++) {
    /* Create a filename to write. */
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    /* Define some variables for PIOc_write_darray. */
    ret = PIOc_def_var(ncid_write, "dummy_scalar_var_int", PIO_INT, 0, NULL, &varid_dummy_scalar_var_int_write); ERR
    if (my_rank == 0)
      printf("varid_dummy_scalar_var_int_write = %d\n", varid_dummy_scalar_var_int_write);

    ret = PIOc_enddef(ncid_write); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    ret = PIOc_inq_varid(ncid_read, "dummy_scalar_var_int", &varid_dummy_scalar_var_int_read); ERR
    if (my_rank == 0)
      printf("varid_dummy_scalar_var_int_read = %d\n", varid_dummy_scalar_var_int_read);

    ret = PIOc_closefile(ncid_read); ERR
  }

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
