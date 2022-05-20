#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <assert.h>
#include <time.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

/* #define PTAPES 12 */
#define PTAPES 2
#define MAXNFLDS 386
#define MAX_CHARS 256

int main(int argc, char* argv[])
{
  int my_rank;
  int ntasks;
  int format = PIO_IOTYPE_ADIOS;
  int niotasks;
  const int ioproc_start = 0;
  int iosysid;
  int ncid_read;

  int varid_fillvalue;
  int varid_meridional_complement;
  int varid_zonal_complement;
  int varid_avgflag;

  int varid_long_name;
  int varid_standard_name;
  int varid_units;
  int varid_sampling_seq;

  PIO_Offset start_2D[2];
  PIO_Offset count_2D[2];

  PIO_Offset start_3D[3];
  PIO_Offset count_3D[3];

  double fillvalue[PTAPES][MAXNFLDS];
  int meridional_complement[PTAPES][MAXNFLDS];
  int zonal_complement[PTAPES][MAXNFLDS];
  char avgflag[PTAPES][MAXNFLDS];

  char long_name[PTAPES][MAXNFLDS][MAX_CHARS];
  char standard_name[PTAPES][MAXNFLDS][MAX_CHARS];
  char units[PTAPES][MAXNFLDS][MAX_CHARS];
  char sampling_seq[PTAPES][MAXNFLDS][MAX_CHARS];

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
  char *filename_eam_r = (char *)malloc(1024);
  if (argc == 1)
      strcpy(filename_eam_r, "F2010_ne4_oQU240_ADIOS.eam.r.0001-01-02-00000.nc");
  else{
      strcpy(filename_eam_r, argv[1]);
  }

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_r, PIO_NOWRITE); ERR

  varid_fillvalue = -1;
  ret = PIOc_inq_varid(ncid_read, "fillvalue", &varid_fillvalue); ERR

  varid_meridional_complement = -1;
  ret = PIOc_inq_varid(ncid_read, "meridional_complement", &varid_meridional_complement); ERR

  varid_zonal_complement = -1;
  ret = PIOc_inq_varid(ncid_read, "zonal_complement", &varid_zonal_complement); ERR

  varid_avgflag = -1;
  ret = PIOc_inq_varid(ncid_read, "avgflag", &varid_avgflag); ERR

  varid_long_name = -1;
  ret = PIOc_inq_varid(ncid_read, "long_name", &varid_long_name); ERR

  varid_standard_name = -1;
  ret = PIOc_inq_varid(ncid_read, "standard_name", &varid_standard_name); ERR

  varid_units = -1;
  ret = PIOc_inq_varid(ncid_read, "units", &varid_units); ERR

  varid_sampling_seq = -1;
  ret = PIOc_inq_varid(ncid_read, "sampling_seq", &varid_sampling_seq); ERR

  for (int t = 0; t < PTAPES; t++) {
    for (int f = 0; f < MAXNFLDS; f++) {
      if (my_rank == 0)
        printf("Calling PIOc_get_var on 8 variables, t = %d, f = %d\n", t, f);

      start_2D[0] = t;
      count_2D[0] = 1;
      start_2D[1] = f;
      count_2D[1] = 1;
      int t0 = clock();
      ret = PIOc_get_vara_double(ncid_read, varid_fillvalue, start_2D, count_2D, &fillvalue[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      ret = PIOc_get_vara_int(ncid_read, varid_meridional_complement, start_2D, count_2D, &meridional_complement[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC);t0 = clock();
      ret = PIOc_get_vara_int(ncid_read, varid_zonal_complement, start_2D, count_2D, &zonal_complement[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      ret = PIOc_get_vara_schar(ncid_read, varid_avgflag, start_2D, count_2D, &avgflag[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      start_3D[0] = t;
      count_3D[0] = 1;
      start_3D[1] = f;
      count_3D[1] = 1;
      start_3D[2] = 0;
      count_3D[2] = MAX_CHARS;

      ret = PIOc_get_vara_text(ncid_read, varid_long_name, start_3D, count_3D, &long_name[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      ret = PIOc_get_vara_text(ncid_read, varid_standard_name, start_3D, count_3D, &standard_name[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      ret = PIOc_get_vara_text(ncid_read, varid_units, start_3D, count_3D, &units[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
      ret = PIOc_get_vara_text(ncid_read, varid_sampling_seq, start_3D, count_3D, &sampling_seq[t][f]); ERR
        if (my_rank == 0) printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
    }
  }
  free(filename_eam_r);
  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
