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
  int dimids[2];
  int iosysid;
  int ncid_write;
  int ncid_read;
  int varid_fractions_ax_ifrac;
  int gdimlen[2] = {1, 4};
  int ioid;
  PIO_Offset *compmap = NULL;
  int element_per_pe;
  double write_buffer[1];
  double read_buffer[1];
  int start_index;

  /* Test filename. */
  char filename[PIO_MAX_NAME];

  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  if (ntasks != 4) {
    if (my_rank == 0)
      printf("This example must be run with 4 MPI tasks!\n");

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return -1;
  }

  int ioproc_stride = 1;
  niotasks = ntasks;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR

  element_per_pe = 1;
  start_index = my_rank + 1;

  compmap = malloc(element_per_pe * sizeof(PIO_Offset));
  for (int i = 0; i < element_per_pe; i++)
    compmap[i] = start_index + i;

  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 2, gdimlen, element_per_pe, compmap, &ioid, NULL, NULL, NULL); ERR

  for (int i = 0; i < 1; i++)
    write_buffer[i] = my_rank + 0.5;

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    ret = PIOc_def_dim(ncid_write, "fractions_ax_ny", (PIO_Offset)gdimlen[0], &dimids[0]); ERR
    ret = PIOc_def_dim(ncid_write, "fractions_ax_nx", (PIO_Offset)gdimlen[1], &dimids[1]); ERR
    ret = PIOc_def_var(ncid_write, "fractions_ax_ifrac", PIO_DOUBLE, 2, dimids, &varid_fractions_ax_ifrac); ERR

    ret = PIOc_enddef(ncid_write); ERR

    ret = PIOc_write_darray(ncid_write, varid_fractions_ax_ifrac, ioid, element_per_pe, write_buffer, NULL); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    for (int i = 0; i < 1; i++)
      read_buffer[i] = -1.0;

    sprintf(filename, "example1_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    varid_fractions_ax_ifrac = -1;
    ret = PIOc_inq_varid(ncid_read, "fractions_ax_ifrac", &varid_fractions_ax_ifrac); ERR

    ret = PIOc_read_darray(ncid_read, varid_fractions_ax_ifrac, ioid, element_per_pe, read_buffer); ERR
    for (int i = 0; i < element_per_pe; i++) {
      if (fabs(read_buffer[i] - write_buffer[i]) > 1E-5) {
        printf("rank = %d, fmt = %d, read wrong data at index %d, expected value: %lf, actual value: %lf\n", my_rank, fmt, i, write_buffer[i], read_buffer[i]);
        break;
      }
    }

    ret = PIOc_closefile(ncid_read); ERR
  }

  free(compmap);
  ret = PIOc_freedecomp(iosysid, ioid); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}

