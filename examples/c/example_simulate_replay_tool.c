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
  int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};
  int niotasks;
  const int ioproc_start = 0;
  int dimids_time_var[2];
  int iosysid;
  int ncid_write;
  int ncid_read;
  int varid_varr0001;
  int ioid;
  float write_buffer[56];
  float read_buffer[56];
  int ndims = 0;
  int *gdimlen = NULL;
  PIO_Offset fmaplen = 0;
  PIO_Offset *compmap = NULL;

  /* Test filename. */
  char filename[PIO_MAX_NAME];

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

  int ioproc_stride = 1;
  niotasks = ntasks;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR

  /* Read the decomp file for testing */
  ret = PIOc_readmap("decomp_1D_866_16p.dat", &ndims, &gdimlen, &fmaplen, &compmap, MPI_COMM_WORLD);
  assert(ndims == 1);
  assert(gdimlen != NULL && gdimlen[0] == 866);
  assert(fmaplen == 56);
  assert(compmap != NULL);

  ret = PIOc_InitDecomp(iosysid, PIO_FLOAT, 1, gdimlen, fmaplen, compmap, &ioid, NULL, NULL, NULL); ERR

  for (int i = 0; i < fmaplen; i++) {
    write_buffer[i] = my_rank * 0.1;
    read_buffer[i] = 0.0;
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    ret = PIOc_def_dim(ncid_write, "time", NC_UNLIMITED, &dimids_time_var[0]); ERR
    ret = PIOc_def_dim(ncid_write, "dim000001", (PIO_Offset)gdimlen[0], &dimids_time_var[1]); ERR
    ret = PIOc_def_var(ncid_write, "varr0001", PIO_FLOAT, 2, dimids_time_var, &varid_varr0001); ERR

    ret = PIOc_enddef(ncid_write); ERR

    PIOc_setframe(ncid_write, varid_varr0001, 0);
    ret = PIOc_write_darray(ncid_write, varid_varr0001, ioid, fmaplen, write_buffer, NULL); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    varid_varr0001 = -1;
    ret = PIOc_inq_varid(ncid_read, "varr0001", &varid_varr0001); ERR

    PIOc_setframe(ncid_read, varid_varr0001, 0);
    ret = PIOc_read_darray(ncid_read, varid_varr0001, ioid, fmaplen, read_buffer); ERR
    for (int i = 0; i < fmaplen; i++) {
      /* We should not compare values correspond to holes in the decomposition map */
      if (compmap[i] > 0) {
        if (fabs(read_buffer[i] - write_buffer[i]) > 1E-5) {
          printf("rank = %d, fmt = %d, read wrong data for float type variable varr0001 at index %d\n", my_rank, fmt, i);
          break;
        }
      }
    }

    ret = PIOc_closefile(ncid_read); ERR
  }

  free(compmap);
  free(gdimlen);
  ret = PIOc_freedecomp(iosysid, ioid); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
