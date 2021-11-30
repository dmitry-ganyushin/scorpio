#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

#define COLUMN_LEN 4096
#define LEVSNO_LEN 5
#define BUDG_FLUX_LEN 30
#define NTAPES 2
#define MAX_CHARS 256

int main(int argc, char* argv[])
{
  int my_rank;
  int ntasks;

  int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};

  int niotasks;

  const int ioproc_start = 0;

  int dimid_column;
  int dimid_levsno;
  int dimid_budg_flux;
  int dimid_ntapes;
  int dimid_max_chars;

  int dimids[2];

  PIO_Offset dimlen;

  int total_dims;

  int gdimlen[2];

  int iosysid;

  int ncid_write;
  int ncid_read;

  int varid_timemgr_rst_type;
  int varid_FSD24_PERIOD;
  int varid_cols1d_wtxy;
  int varid_DZSNO;
  int varid_locfnh;
  int varid_locfnhr;
  int varid_budg_fluxG;

  /* Sample data for scalar variables. */
  int put_scalar_int_data = -1;
  int get_scalar_int_data = -1;

  PIO_Offset start_1D[1];
  PIO_Offset count_1D[1];

  PIO_Offset start_2D[2];
  PIO_Offset count_2D[2];

  int ioid_1D;
  int ioid_2D;

  /* Buffers for sample write/read darray data. */
  double write_darray_buffer_1D[COLUMN_LEN];
  double read_darray_buffer_1D[COLUMN_LEN];
  double write_darray_buffer_2D[COLUMN_LEN * LEVSNO_LEN];
  double read_darray_buffer_2D[COLUMN_LEN * LEVSNO_LEN];

  /* Buffers for sample put/get var data. */
  double put_var_buffer_1D[BUDG_FLUX_LEN];
  double get_var_buffer_1D[BUDG_FLUX_LEN];
  char put_string_array[NTAPES][MAX_CHARS];
  char get_string_array[NTAPES][MAX_CHARS];

  PIO_Offset *compmap_1D;
  PIO_Offset *compmap_2D;

  int element_per_pe_1D;
  int element_per_pe_2D;

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

  element_per_pe_1D = COLUMN_LEN / ntasks;
  compmap_1D = malloc(element_per_pe_1D * sizeof(PIO_Offset));
  for (int i = 0; i < element_per_pe_1D; i++)
    compmap_1D[i] = my_rank * element_per_pe_1D + i + 1;

  gdimlen[0] = COLUMN_LEN;
  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 1, gdimlen, element_per_pe_1D, compmap_1D, &ioid_1D, NULL, NULL, NULL); ERR
  free(compmap_1D);

  element_per_pe_2D = (COLUMN_LEN * LEVSNO_LEN) / ntasks;
  compmap_2D = malloc(element_per_pe_2D * sizeof(PIO_Offset));
  for (int i = 0; i < element_per_pe_2D; i++)
    compmap_2D[i] = my_rank * element_per_pe_2D + i + 1;

  gdimlen[0] = COLUMN_LEN;
  gdimlen[1] = LEVSNO_LEN;
  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 2, gdimlen, element_per_pe_2D, compmap_2D, &ioid_2D, NULL, NULL, NULL); ERR
  free(compmap_2D);

  /* Prepare sample data for write buffers and initialize read buffers. */
  for (int i = 0; i < element_per_pe_1D; i++) {
    write_darray_buffer_1D[i] = my_rank;
    read_darray_buffer_1D[i] = 0;
  }

  for (int i = 0; i < element_per_pe_2D; i++) {
    write_darray_buffer_2D[i] = my_rank;
    read_darray_buffer_2D[i] = 0;
  }

  for (int i = 0; i < BUDG_FLUX_LEN; i++) {
    put_var_buffer_1D[i] = (i + 1) * 0.1;
    get_var_buffer_1D[i] = 0.0;
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    /* Create a filename to write. */
    sprintf(filename, "test_elm.r_%d.nc", fmt);

    ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER); ERR

    /* Define two scalar variables for PIOc_put_vars. */
    ret = PIOc_def_var(ncid_write, "timemgr_rst_type", PIO_INT, 0, NULL, &varid_timemgr_rst_type); ERR
    ret = PIOc_def_var(ncid_write, "FSD24_PERIOD", PIO_INT, 0, NULL, &varid_FSD24_PERIOD); ERR

    /* Define dimensions. */
    ret = PIOc_def_dim(ncid_write, "column", COLUMN_LEN, &dimid_column); ERR
    ret = PIOc_def_dim(ncid_write, "levsno", LEVSNO_LEN, &dimid_levsno); ERR
    ret = PIOc_def_dim(ncid_write, "budg_flux", BUDG_FLUX_LEN, &dimid_budg_flux); ERR
    ret = PIOc_def_dim(ncid_write, "ntapes", NTAPES, &dimid_ntapes); ERR
    ret = PIOc_def_dim(ncid_write, "max_chars", MAX_CHARS, &dimid_max_chars); ERR

    /* Define two double variables for PIOc_write_darray. */
    dimids[0] = dimid_column;
    dimids[1] = dimid_levsno;
    ret = PIOc_def_var(ncid_write, "cols1d_wtxy", PIO_DOUBLE, 1, dimids, &varid_cols1d_wtxy); ERR
    ret = PIOc_def_var(ncid_write, "DZSNO", PIO_DOUBLE, 2, dimids, &varid_DZSNO); ERR

    /* Define two text variables for PIOc_put_vars. */
    dimids[0] = dimid_ntapes;
    dimids[1] = dimid_max_chars;
    ret = PIOc_def_var(ncid_write, "locfnh", PIO_CHAR, 2, dimids, &varid_locfnh); ERR
    ret = PIOc_def_var(ncid_write, "locfnhr", PIO_CHAR, 2, dimids, &varid_locfnhr); ERR

    /* Define one double variables for PIOc_put_vars. */
    dimids[0] = dimid_budg_flux;
    ret = PIOc_def_var(ncid_write, "budg_fluxG", PIO_DOUBLE, 1, dimids, &varid_budg_fluxG); ERR

    ret = PIOc_enddef(ncid_write); ERR

    /* Put one scalar variable. */
    put_scalar_int_data = 1;
    ret = PIOc_put_var_int(ncid_write, varid_timemgr_rst_type, &put_scalar_int_data); ERR

    /* Write one 1D double variable. */
    ret = PIOc_write_darray(ncid_write, varid_cols1d_wtxy, ioid_1D, element_per_pe_1D, write_darray_buffer_1D, NULL); ERR

    /* Write one 2D double variable. */
    ret = PIOc_write_darray(ncid_write, varid_DZSNO, ioid_2D, element_per_pe_2D, write_darray_buffer_2D, NULL); ERR

    /* Put one scalar variable. */
    put_scalar_int_data = 12;
    ret = PIOc_put_var_int(ncid_write, varid_FSD24_PERIOD, &put_scalar_int_data); ERR

    /* Put two text variables. */
    start_2D[0] = 0;
    count_2D[0] = NTAPES;
    start_2D[1] = 0;
    count_2D[1] = MAX_CHARS;

    /* Empty strings for locfnh */
    memset(put_string_array, 0, sizeof(put_string_array));
    ret = PIOc_put_vara_text(ncid_write, varid_locfnh, start_2D, count_2D, put_string_array[0]); ERR

    /* Non-empty strings for locfnhr */
    for (int t = 0; t < NTAPES; t++)
      snprintf(put_string_array[t], MAX_CHARS, "tape_%d.nc", t);
    ret = PIOc_put_vara_text(ncid_write, varid_locfnhr, start_2D, count_2D, put_string_array[0]); ERR
 
    /* Put one 1D double variable. */
    ret = PIOc_put_var_double(ncid_write, varid_budg_fluxG, put_var_buffer_1D); ERR

    ret = PIOc_closefile(ncid_write); ERR
  }

  for (int fmt = 0; fmt < 2; fmt++) {
    sprintf(filename, "test_elm.r_%d.nc", fmt);

    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    /* Get one scalar variable. */
    varid_timemgr_rst_type = -1;
    ret = PIOc_inq_varid(ncid_read, "timemgr_rst_type", &varid_timemgr_rst_type); ERR

    get_scalar_int_data = -1;
    ret = PIOc_get_var_int(ncid_read, varid_timemgr_rst_type, &get_scalar_int_data); ERR
    if (get_scalar_int_data != 1)
      printf("rank = %d, read wrong data for timemgr_rst_type\n", my_rank);

    /* Close file and open it later. */
    ret = PIOc_closefile(ncid_read); ERR

    /* Reopen file. */
    ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE); ERR

    total_dims = 0;
    ret = PIOc_inq_ndims(ncid_read, &total_dims); ERR

    dimid_column = -1;
    ret = PIOc_inq_dimid(ncid_read, "column", &dimid_column); ERR

    if (dimid_column < 0 || dimid_column > total_dims)
        printf("rank = %d, read wrong ID for dimension column\n", my_rank);

    dimlen = -1;
    ret = PIOc_inq_dimlen(ncid_read, dimid_column, &dimlen); ERR
    if (dimlen != COLUMN_LEN)
        printf("rank = %d, read wrong length for dimension column\n", my_rank);

    varid_DZSNO = -1;
    ret = PIOc_inq_varid(ncid_read, "DZSNO", &varid_DZSNO); ERR

    dimids[0] = -1;
    dimids[1] = -1;
    ret = PIOc_inq_vardimid(ncid_read, varid_DZSNO, dimids); ERR

    if (dimids[0] < 0 || dimids[0] > total_dims)
        printf("rank = %d, read wrong ID for dimension column\n", my_rank);

    dimlen = -1;
    ret = PIOc_inq_dimlen(ncid_read, dimids[0], &dimlen); ERR
    if (dimlen != COLUMN_LEN)
          printf("rank = %d, read wrong length for dimension column\n", my_rank);

    if (dimids[1] < 0 || dimids[1] > total_dims)
        printf("rank = %d, read wrong ID for dimension levsno\n", my_rank);

    dimlen = -1;
    ret = PIOc_inq_dimlen(ncid_read, dimids[1], &dimlen); ERR
    if (dimlen != LEVSNO_LEN)
          printf("rank = %d, read wrong length for dimension levsno\n", my_rank);

    /* Read one 1D double variable. */
    for (int i = 0; i < element_per_pe_1D; i++)
      read_darray_buffer_1D[i] = 0.0;

    ret = PIOc_inq_varid(ncid_read, "cols1d_wtxy", &varid_cols1d_wtxy); ERR
    ret = PIOc_read_darray(ncid_read, varid_cols1d_wtxy, ioid_1D, read_darray_buffer_1D, read_darray_buffer_1D); ERR

    for (int i = 0; i < element_per_pe_1D; i++) {
      double diff = read_darray_buffer_1D[i] - write_darray_buffer_1D[i];
      if (fabs(diff) > 1E-5) {
          printf("rank = %d, read wrong data for cols1d_wtxy at index %d\n", my_rank, i);
          break;
      }
    }

    /* Read one 2D double variable. */
    for (int i = 0; i < element_per_pe_2D; i++)
      read_darray_buffer_2D[i] = 0.0;

    ret = PIOc_inq_varid(ncid_read, "DZSNO", &varid_DZSNO); ERR
    ret = PIOc_read_darray(ncid_read, varid_DZSNO, ioid_2D, read_darray_buffer_2D, read_darray_buffer_2D); ERR

    for (int i = 0; i < element_per_pe_2D; i++) {
      double diff = read_darray_buffer_2D[i] - write_darray_buffer_2D[i];
      if (fabs(diff) > 1E-5) {
          printf("rank = %d, read wrong data for DZSNO at index %d\n", my_rank, i);
          break;
      }
    }

    /* Get one scalar variable. */
    varid_FSD24_PERIOD = -1;
    ret = PIOc_inq_varid(ncid_read, "FSD24_PERIOD", &varid_FSD24_PERIOD); ERR

    get_scalar_int_data = -1;
    ret = PIOc_get_var_int(ncid_read, varid_FSD24_PERIOD, &get_scalar_int_data); ERR
    if (get_scalar_int_data != 12)
      printf("rank = %d, read wrong data for FSD24_PERIOD\n", my_rank);

    /* Get two text variables. */
    varid_locfnh = -1;
    ret = PIOc_inq_varid(ncid_read, "locfnh", &varid_locfnh); ERR

    memset(get_string_array, 0, sizeof(get_string_array));
    ret = PIOc_get_var_text(ncid_read, varid_locfnh, get_string_array[0]); ERR

    for (int t = 0; t < NTAPES; t++) {
      if (strlen(get_string_array[t]) > 0)
        printf("rank = %d, tape = %d, read wrong text (expected: empty, actual: %s) for locfnh\n", my_rank, t, get_string_array[t]);
    }

    memset(get_string_array, 0, sizeof(get_string_array));
    ret = PIOc_get_var_text(ncid_read, varid_locfnhr, get_string_array[0]); ERR

    for (int t = 0; t < NTAPES; t++) {
      if (strncmp(get_string_array[t], put_string_array[t], MAX_CHARS))
        printf("rank = %d, tape = %d, read wrong text (expected: %s, actual: %s) for locfnhr\n", my_rank, t, put_string_array[t], get_string_array[t]);
    }

    /* Get one 1D double variable. */
    varid_budg_fluxG = -1;

    ret = PIOc_inq_varid(ncid_read, "budg_fluxG", &varid_budg_fluxG); ERR
    for (int i = 0; i < BUDG_FLUX_LEN; i++) get_var_buffer_1D[i] = -1;
    ret = PIOc_get_var_double(ncid_read, varid_budg_fluxG, get_var_buffer_1D); ERR
    for (int i = 0; i < BUDG_FLUX_LEN; i++) {
      double diff = get_var_buffer_1D[i] - put_var_buffer_1D[i];
      if (fabs(diff) > 1E-5) {
          printf("rank = %d, get wrong data for budg_fluxG at index %d\n", my_rank, i);
          break;
      }
    }

    ret = PIOc_closefile(ncid_read); ERR
  }

  ret = PIOc_freedecomp(iosysid, ioid_1D); ERR
  ret = PIOc_freedecomp(iosysid, ioid_2D); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
