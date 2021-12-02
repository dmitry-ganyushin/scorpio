#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

#define MAX_NFLDS 40
#define MAX_FNAME_LEN 16

int main(int argc, char* argv[])
{
  int my_rank;
  int ntasks;

  int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};

  int niotasks;

  const int ioproc_start = 0;

  int dimid_max_nflds;
  int dimid_max_fname_len;

  int dimids[2];

  int iosysid;

  int ncid_write;
  int ncid_read;

  int varid_names;

  PIO_Offset start_2D[2];
  PIO_Offset count_2D[2];

  /* Buffers for sample put/get var data. */
  char put_all_strings_buffer[MAX_NFLDS][MAX_FNAME_LEN];
  char get_single_string_buffer[MAX_FNAME_LEN];

  char filename[PIO_MAX_NAME];

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

  /* Create name strings */
  for (int fld = 0; fld < MAX_NFLDS; fld++)
    snprintf(put_all_strings_buffer[fld], MAX_FNAME_LEN, "FIELD_%d", fld); 

  for (int fmt = 0; fmt < 2; fmt++) {
    /* Create a filename to write. */
    sprintf(filename, "test_put_get_strings_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    /* Define dimensions. */
    ret = PIOc_def_dim(ncid_write, "max_nflds", MAX_NFLDS, &dimid_max_nflds); ERR
    ret = PIOc_def_dim(ncid_write, "max_fname_len", MAX_FNAME_LEN, &dimid_max_fname_len); ERR

    /* Define a text variable for PIOc_put_vars. */
    dimids[0] = dimid_max_nflds;
    dimids[1] = dimid_max_fname_len;
    ret = PIOc_def_var(ncid_write, "field_names", PIO_CHAR, 2, dimids, &varid_names); ERR

    ret = PIOc_enddef(ncid_write); ERR

    /* Put all strings to text variable field_names */
    ret = PIOc_put_var_text(ncid_write, varid_names, put_all_strings_buffer); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "test_put_get_strings_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    /* Get text variable */
    varid_names = -1;
    ret = PIOc_inq_varid(ncid_read, "field_names", &varid_names); ERR

    /* Always get full string */
    start_2D[1] = 0;
    count_2D[1] = MAX_FNAME_LEN;

    for (int fld = 0; fld < MAX_NFLDS; fld++) {
      start_2D[0] = fld; /* start for ith string */
      count_2D[0] = 1; /* Get only one string each time */
      memset(get_single_string_buffer, 0, MAX_FNAME_LEN);
      ret = PIOc_get_vara_text(ncid_read, varid_names, start_2D, count_2D, get_single_string_buffer); ERR

      //printf("rank = %d, field = %d, string read: %s\n", my_rank, fld, get_single_string_buffer);
      if (strncmp(get_single_string_buffer, put_all_strings_buffer[fld], MAX_FNAME_LEN))
        printf("rank = %d, field = %d, read wrong name string (expected: %s, actual: %s)\n", my_rank, fld, put_all_strings_buffer[fld], get_single_string_buffer);
    }

    ret = PIOc_closefile(ncid_read); ERR
  }

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
