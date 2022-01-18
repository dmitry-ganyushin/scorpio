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
  int varid_indexToCellID;
  int ndims = 0;
  int *gdimlen = NULL;
  PIO_Offset fmaplen = 0;
  PIO_Offset *compmap = NULL;
  int ioid;

  int read_indexToCellID[448];

  int get_indexToCellID[7153];

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
  ret = PIOc_readmap("decomp_1D_7153_16p.dat", &ndims, &gdimlen, &fmaplen, &compmap, MPI_COMM_WORLD); ERR
  assert(ndims == 1);
  assert(gdimlen != NULL && gdimlen[0] == 7153);
  assert(fmaplen == 447 || fmaplen == 448);
  assert(compmap != NULL);

  ret = PIOc_InitDecomp(iosysid, PIO_INT, 1, gdimlen, fmaplen, compmap, &ioid, NULL, NULL, NULL); ERR

  char filename_mpassi[] = "F2010_ne4_oQU240_ADIOS.mpassi.rst.0001-01-02_00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_mpassi, PIO_NOWRITE); ERR

  varid_indexToCellID = -1;
  ret = PIOc_inq_varid(ncid_read, "indexToCellID", &varid_indexToCellID); ERR

  ret = PIOc_read_darray(ncid_read, varid_indexToCellID, ioid, fmaplen, read_indexToCellID); ERR
  for (int i = 0; i < fmaplen; i++) {
    if (read_indexToCellID[i] != compmap[i]) {
      printf("rank = %d, read wrong data for indexToCellID at index %d\n", my_rank, i);
      break;
    }
  }

  for (int i = 0; i < 7153; i++)
    get_indexToCellID[i] = -1;

  ret = PIOc_get_var_int(ncid_read, varid_indexToCellID, get_indexToCellID); ERR
  for (int i = 0; i < 7153; i++) {
    if (get_indexToCellID[i] != (i + 1)) {
      printf("rank = %d, get wrong data for indexToCellID at index %d\n", my_rank, i);
      break;
    }
  }

  ret = PIOc_closefile(ncid_read); ERR

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
