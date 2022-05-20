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
  int varid_a2x_ax_Sa_topo;
  int gdimlen[2] = {1, 866};
  int ioid;
  PIO_Offset *compmap;
  int element_per_pe;
  double read_darray_buffer[55];
  int start_index;
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

  if (my_rank == 0 || my_rank == 8)
    element_per_pe = 55;
  else
    element_per_pe = 54;

  if (my_rank == 0)
    start_index = 1;
  else if (my_rank == 1)
    start_index = 56;
  else if (my_rank == 2)
    start_index = 110;
  else if (my_rank == 3)
    start_index = 164;
  else if (my_rank == 4)
    start_index = 218;
  else if (my_rank == 5)
    start_index = 272;
  else if (my_rank == 6)
    start_index = 326;
  else if (my_rank == 7)
    start_index = 380;
  else if (my_rank == 8)
    start_index = 434;
  else if (my_rank == 9)
    start_index = 489;
  else if (my_rank == 10)
    start_index = 543;
  else if (my_rank == 11)
    start_index = 597;
  else if (my_rank == 12)
    start_index = 651;
  else if (my_rank == 13)
    start_index = 705;
  else if (my_rank == 14)
    start_index = 759;
  else if (my_rank == 15)
    start_index = 813;

  compmap = malloc(element_per_pe * sizeof(PIO_Offset));
  for (int i = 0; i < element_per_pe; i++)
    compmap[i] = start_index + i;

  ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, 2, gdimlen, element_per_pe, compmap, &ioid, NULL, NULL, NULL); ERR
  free(compmap);

  char filename_cpl_r[] = "F2010_ne4_oQU240_ADIOS.cpl.r.0001-01-02-00000.nc";

  ret = PIOc_openfile(iosysid, &ncid_read, &format, filename_cpl_r, PIO_NOWRITE); ERR

  varid_a2x_ax_Sa_topo = -1;
  ret = PIOc_inq_varid(ncid_read, "a2x_ax_Sa_topo", &varid_a2x_ax_Sa_topo); ERR
  ret = PIOc_read_darray(ncid_read, varid_a2x_ax_Sa_topo, ioid, element_per_pe, read_darray_buffer); ERR

  ret = PIOc_closefile(ncid_read); ERR

  ret = PIOc_freedecomp(iosysid, ioid); ERR

  ret = PIOc_finalize(iosysid); ERR

  MPI_Finalize();

#ifdef TIMING
  GPTLfinalize();
#endif

  return 0;
}
