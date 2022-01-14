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
  int varid_static_ener_ac;
  int varid_tracer_cnst_curr_fname;
  int varid_fincl;
  int ndims = 0;
  int *gdimlen = NULL;
  PIO_Offset fmaplen = 0;
  PIO_Offset *compmap = NULL;
  int ioid;

  double read_darray_buffer[56];
  double get_att_float_data;
  char get_text_data[12 * 1000 * 26] = "\0";

  int varid_mdimnames;
  int dimids[2];
  PIO_Offset dimlen;
  PIO_Offset start_2D[2];
  PIO_Offset count_2D[2];
  char get_mdimnames_1st_str[16] = "\0";
  char get_mdimnames_2nd_str[16] = "\0";
  char get_mdimnames_both_strs[2 * 16] = "\0";

  int varid_mdims;
  int get_mdims[12 * 386 * 1] = {0};

  int varid_decomp_type;
  int get_decomp_type[12 * 386] = {0};

  int varid_field_name;
  char get_field_name[12][386][27];

  int ret = PIO_NOERR;

#ifdef TIMING
  GPTLinitialize();
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#if 1
  if (ntasks != 16) {
    if (my_rank == 0)
      printf("This example must be run with 16 MPI tasks!\n");

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return -1;
  }
#endif
  int ioproc_stride = 4;
  niotasks = ntasks / ioproc_stride;

  ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid); ERR
#if 1
  /* Read the decomp file for testing */
  ret = PIOc_readmap("decomp_1D_866_16p.dat", &ndims, &gdimlen, &fmaplen, &compmap, MPI_COMM_WORLD);
  assert(ndims == 1);
  assert(gdimlen != NULL && gdimlen[0] == 866);
  assert(fmaplen == 56);
  assert(compmap != NULL);

  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 1, gdimlen, fmaplen, compmap, &ioid, NULL, NULL, NULL); ERR
  free(compmap);
  free(gdimlen);
#endif
  char filename_eam_r[] = "F2010_ne4_oQU240_ADIOS.eam.r.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_r, PIO_NOWRITE); ERR
#if 1
  varid_static_ener_ac = -1;
  ret = PIOc_inq_varid(ncid_read, "static_ener_ac", &varid_static_ener_ac); ERR

  ret = PIOc_read_darray(ncid_read, varid_static_ener_ac, ioid, fmaplen, read_darray_buffer); ERR
#endif
  varid_tracer_cnst_curr_fname = -1;
  ret = PIOc_inq_varid(ncid_read, "tracer_cnst_curr_fname", &varid_tracer_cnst_curr_fname); ERR

  get_att_float_data = -1.0;
  ret = PIOc_get_att(ncid_read, varid_tracer_cnst_curr_fname, "offset_time", &get_att_float_data); ERR

  varid_fincl = -1;
  ret = PIOc_inq_varid(ncid_read, "fincl", &varid_fincl); ERR

  ret = PIOc_get_var_text(ncid_read, varid_fincl, get_text_data); ERR

  varid_mdimnames = -1;
  ret = PIOc_inq_varid(ncid_read, "mdimnames", &varid_mdimnames); ERR

  dimids[0] = -1;
  dimids[1] = -1;
  ret = PIOc_inq_vardimid(ncid_read, varid_mdimnames, dimids); ERR

  dimlen = -1;
  ret = PIOc_inq_dimlen(ncid_read, dimids[0], &dimlen); ERR
  assert(dimlen == 2);

  dimlen = -1;
  ret = PIOc_inq_dimlen(ncid_read, dimids[1], &dimlen); ERR
  assert(dimlen == 16);

  start_2D[0] = 0;
  count_2D[0] = 1;
  start_2D[1] = 0;
  count_2D[1] = 16;
#if 1
  /* Get only 1st string, no heap-buffer-overflow, but the returned string is incorrect (expected: lev, actual: ilev) */
  ret = PIOc_get_vara_text(ncid_read, varid_mdimnames, start_2D, count_2D, get_mdimnames_1st_str); ERR
  if (my_rank == 0)
    printf("get_mdimnames_1st_str = %s\n", get_mdimnames_1st_str);
#endif
#if 1
  start_2D[0] = 1;
  /* Get only 2nd string failed, AddressSanitizer: heap-buffer-overflow scorpio/src/clib/pio_getput_int.c:1625 in PIOc_get_vars_tc */
  ret = PIOc_get_vara_text(ncid_read, varid_mdimnames, start_2D, count_2D, get_mdimnames_2nd_str); ERR
  if (my_rank == 0)
    printf("get_mdimnames_2nd_str = %s\n", get_mdimnames_2nd_str);
#endif

#if 1
  /* Get both strings together failed, AddressSanitizer: heap-buffer-overflow scorpio/src/clib/pio_getput_int.c:1632 in PIOc_get_vars_tc */
  ret = PIOc_get_var_text(ncid_read, varid_mdimnames, get_mdimnames_both_strs); ERR
  if (my_rank == 0)
    printf("get_mdimnames_both_strs = %s\n", get_mdimnames_both_strs);
#endif

  /* Get one 3D int variable. */
  varid_mdims = -1;
  ret = PIOc_inq_varid(ncid_read, "mdims", &varid_mdims); ERR

  for (int i = 0; i < 12 * 386 * 1; i++)
    get_mdims[i] = -2;

  ret = PIOc_get_var_int(ncid_read, varid_mdims, get_mdims); ERR
  if (my_rank == 0) {
    for (int i = 0; i < 12 * 386 * 1; i++)
      printf("get_mdims[%d] = %d\n", i, get_mdims[i]);
  }

  /* Get one 2D int variable. */
  varid_decomp_type = -1;
  ret = PIOc_inq_varid(ncid_read, "decomp_type", &varid_decomp_type); ERR

  for (int i = 0; i < 12 * 386; i++)
    get_decomp_type[i] = -2;

  ret = PIOc_get_var_int(ncid_read, varid_decomp_type, get_decomp_type); ERR
  if (my_rank == 0) {
    for (int i = 0; i < 12 * 386; i++)
      printf("get_decomp_type[%d] = %d\n", i, get_decomp_type[i]);
  }

  varid_field_name = -1;
  ret = PIOc_inq_varid(ncid_read, "field_name", &varid_field_name); ERR

  memset(get_field_name, 0, sizeof(get_field_name));
  ret = PIOc_get_var_text(ncid_read, varid_field_name, get_field_name); ERR
  if (my_rank == 0) {
    for (int f = 0; f < 20; f++)
      printf("get_field_name[0][%d] = %s\n", f, get_field_name[0][f]);
  }

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_freedecomp(iosysid, ioid); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
