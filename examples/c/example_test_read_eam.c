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
  int ndims = 0;
  int *gdimlen = NULL;
  PIO_Offset fmaplen = 0;
  PIO_Offset *compmap = NULL;
  int ioid;

  double read_darray_buffer[56];
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

  /* Read the decomp file for testing */
  ret = PIOc_readmap("decomp_1D_866_16p.dat", &ndims, &gdimlen, &fmaplen, &compmap, MPI_COMM_WORLD);
  assert(ndims == 1);
  assert(gdimlen != NULL && gdimlen[0] == 866);
  assert(fmaplen == 56);
  assert(compmap != NULL);

  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 1, gdimlen, fmaplen, compmap, &ioid, NULL, NULL, NULL); ERR
  free(compmap);
  free(gdimlen);

  char filename_eam_r[] = "F2010_ne4_oQU240_ADIOS.eam.r.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_eam_r, PIO_NOWRITE); ERR

  varid_static_ener_ac = -1;
  ret = PIOc_inq_varid(ncid_read, "static_ener_ac", &varid_static_ener_ac); ERR

  ret = PIOc_read_darray(ncid_read, varid_static_ener_ac, ioid, fmaplen, read_darray_buffer); ERR

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_freedecomp(iosysid, ioid); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
