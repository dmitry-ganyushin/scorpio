#include <stdio.h>
#include <mpi.h>
#include <pio.h>

#ifdef TIMING
#include <gptl.h>
#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

int main(int argc, char *argv[]) {
    int my_rank;
    int ntasks;
    int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};
    int niotasks;
    const int ioproc_start = 0;
    int iosysid;
    int ncid_write1;
    int ncid_write2;
    int ncid_read1;

    char filename1[PIO_MAX_NAME];
    char filename2[PIO_MAX_NAME];

    int ret = PIO_NOERR;

#ifdef TIMING
    GPTLinitialize();
#endif

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    int ioproc_stride = 1;
    niotasks = ntasks;

    ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start, PIO_REARR_SUBSET, &iosysid);
    ERR

    for (int fmt = 0; fmt < 2; fmt++) {
        sprintf(filename1, "testfile1_%d.nc", fmt);
        ret = PIOc_createfile(iosysid, &ncid_write1, &(formats[fmt]), filename1, PIO_CLOBBER);
        ERR

        sprintf(filename2, "testfile2_%d.nc", fmt);
        ret = PIOc_createfile(iosysid, &ncid_write2, &(formats[fmt]), filename2, PIO_CLOBBER);
        ERR

        ret = PIOc_closefile(ncid_write1);
        ERR

        ret = PIOc_openfile(iosysid, &ncid_read1, &(formats[fmt]), filename1, PIO_NOWRITE);
        ERR

        ret = PIOc_closefile(ncid_write2);
        ERR

        ret = PIOc_closefile(ncid_read1);
        ERR
    }

    ret = PIOc_finalize(iosysid);
    ERR

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}

