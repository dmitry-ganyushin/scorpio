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
  int my_rank;
  int ntasks;
  int format = PIO_IOTYPE_ADIOS;
  int niotasks;
  const int ioproc_start = 0;
  int iosysid;
  int ncid_read;
  int ncid_read1;
  int ncid_read2;
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

  char filename_elm_r[] = "F2010_ne4_oQU240_ADIOS.elm.r.0001-01-02-00000.nc";
  char filename_elm_rh0[] = "F2010_ne4_oQU240_ADIOS.elm.rh0.0001-01-02-00000.nc";

  char filename_eam_r[] = "F2010_ne4_oQU240_ADIOS-gnu.eam.r.0001-01-02-00000.nc";

  char filename_cpl_r[] = "F2010_ne4_oQU240_ADIOS-gnu.cpl.r.0001-01-02-00000.nc";

  char filename_h0[] = "F2010_ne4_oQU240_ADIOS-gnu.eam.h0.0001-01-01-00000.nc";

  // Max current_var_cnt = 3242
//  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_elm_r, PIO_NOWRITE); ERR
//  ret = PIOc_closefile(ncid_read); ERR

  // Max current_var_cnt = 1678
  // heap-buffer-overflow: attr_data buffer size (256 bytes) is too small to read this file
  // [Workaround] In pioc_support.c, line 3186: char *attr_data = calloc(sizeof(char), 64 * PIO_MAX_NAME);
//  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_elm_rh0, PIO_NOWRITE); ERR
//  ret = PIOc_closefile(ncid_read); ERR

  // Max current_var_cnt = 72
//  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_r, PIO_NOWRITE); ERR
//  ret = PIOc_closefile(ncid_read); ERR
  // Max current_var_cnt = 18214: larger than 8192 which can corrupt struct file_desc_t
  // [Workaround] In pio.h, line 962: struct adios_att_desc_t adios_attrs[4 * PIO_MAX_ATTRS];
  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_cpl_r, PIO_NOWRITE); ERR
  ret = PIOc_closefile(ncid_read); ERR

    ret = PIOc_openfile(iosysid, &ncid_read1, &format, filename_h0, PIO_NOWRITE); ERR
    ret = PIOc_closefile(ncid_read1); ERR

    ret = PIOc_openfile(iosysid, &ncid_read1, &format, filename_h0, PIO_NOWRITE); ERR
    ret = PIOc_closefile(ncid_read1); ERR

    ret = PIOc_openfile(iosysid, &ncid_read2, &format, filename_eam_r, PIO_NOWRITE); ERR
    ret = PIOc_closefile(ncid_read2); ERR

    ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_cpl_r, PIO_NOWRITE); ERR
    ret = PIOc_closefile(ncid_read); ERR

    ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_cpl_r, PIO_NOWRITE); ERR
    ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
