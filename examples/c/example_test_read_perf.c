#include <stdio.h>
#include <mpi.h>
#include <pio.h>
#include <math.h>
#include <time.h>

#ifdef TIMING

#include <gptl.h>

#endif

#define ERR { if (ret != PIO_NOERR) printf("rank = %d, error at line = %d\n", my_rank, __LINE__); }

/* Number of elements of the local darray data that
   will be handled by each processor. */
#define ELEMENTS_PER_PE 1024 * 256

/* Length of variables to be put/get */
#define PUT_GET_VAR_LEN 10

int main(int argc, char *argv[]) {
    /* Zero-based rank of processor. */
    int my_rank;

    /* Number of processors involved in current execution. */
    int ntasks;

    /* This example only tests two formats (IO types): PnetCDF and ADIOS */
    int formats[2] = {PIO_IOTYPE_PNETCDF, PIO_IOTYPE_ADIOS};

    /* Number of processors that will do IO. In this example we
       will do IO from all processors. */
    int niotasks;

    /* Zero based rank of first processor to be used for I/O. */
    const int ioproc_start = 0;

    /* The dimension IDs. */
    int dimid_text_var_len;
    int dimid_darray_var_len;
    int dimid_put_get_var_len;

    /* This simple example uses one-dimensional data. */
    const int NDIMS = 1;

    /* Lengths of the global dimensions */
    int gdimlen[NDIMS];

    /* The ID for the parallel I/O system. It is set by
       PIOc_Init_Intracomm(). It references an internal structure
       containing the general IO subsystem data and MPI
       structure. It is passed to PIOc_finalize() to free
       associated resources, after all I/O, but before
       MPI_Finalize is called. */
    int iosysid;

    /* The ncid of the netCDF file created and read in this example. */
    int ncid_write = 0;
    int ncid_read = 0;

    /* The IDs of the netCDF variables in the example file. */
    int varid_dummy_scalar_var_int;
    int varid_dummy_scalar_var_float;

    int varid_dummy_text_var;

    int varid_dummy_darray_var_int;
    int varid_dummy_darray_var_float;
    int varid_dummy_darray_var_float1;
    int varid_dummy_darray_var_float2;
    int varid_dummy_darray_var_float3;
    int varid_dummy_darray_var_float4;
    int varid_dummy_darray_var_float5;
    int varid_dummy_darray_var_float6;
    int varid_dummy_darray_var_float7;
    int varid_dummy_darray_var_float8;
    int varid_dummy_darray_var_float9;
    int varid_dummy_put_get_var_int;
    int varid_dummy_put_get_var_float;

    /* start/count arrays for get/put var */
    PIO_Offset start[NDIMS];
    PIO_Offset count[NDIMS];

    /* The I/O description IDs as passed back by PIOc_InitDecomp()
       and freed in PIOc_freedecomp(). */
    int ioid_int;
    int ioid_double;

    /* Sample data for global attributes. */
    int put_global_att_int_data = -10;
    float put_global_att_float_data = -9.9f;
    char put_global_att_text_data[PIO_MAX_NAME] = "Dummy global attribute string";

    /* Sample data for variable attributes. */
    int put_att_int_data = -100;
    float put_att_float_data = -99.99f;
    char put_att_text_data[PIO_MAX_NAME] = "Dummy variable attribute string";

    /* Sample data for scalar variables. */
    int put_scalar_int_data = -1000;
    float put_scalar_float_data = -999.999f;

    /* Sample string for text variable. */
    char put_text_data[PIO_MAX_NAME] = "Dummy text variable string";

    /* Buffers for sample write/read darray data. */
    double write_darray_buffer_double[ELEMENTS_PER_PE];
    double read_darray_buffer_double[ELEMENTS_PER_PE];

    /* A 1-D array which holds the decomposition mapping for this example. */
    PIO_Offset *compmap;

    /* Test filename. */
    char filename[PIO_MAX_NAME];


    int ret = PIO_NOERR;

#ifdef TIMING
    GPTLinitialize();
#endif

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /* Set lengths of the global dimensions (1D in this example). */
    gdimlen[0] = ELEMENTS_PER_PE * ntasks;

    int ioproc_stride = 4;
    if (ntasks == 1) ioproc_stride = 1;
    niotasks = ntasks / ioproc_stride;

    /* Initialize the PIO IO system. This specifies how
       many and which processors are involved in I/O. */
    ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride, ioproc_start,
                              PIO_REARR_SUBSET, &iosysid);
    ERR

    /* Describe the decomposition. This is a 1-based array. */
    compmap = malloc(ELEMENTS_PER_PE * sizeof(PIO_Offset));
    for (int i = 0; i < ELEMENTS_PER_PE; i++)
        compmap[i] = my_rank * ELEMENTS_PER_PE + i + 1;

    /* Create the decomposition for this example. */
    ret = PIOc_InitDecomp(iosysid, PIO_INT, NDIMS, gdimlen, ELEMENTS_PER_PE, compmap, &ioid_int, NULL, NULL, NULL);
    ERR
    ret = PIOc_InitDecomp(iosysid, PIO_DOUBLE, NDIMS, gdimlen, ELEMENTS_PER_PE, compmap, &ioid_double, NULL, NULL,
                          NULL);
    ERR
    free(compmap);

    /* Prepare sample data for write buffers and initialize read buffers. */
    for (int i = 0; i < ELEMENTS_PER_PE; i++) {
        write_darray_buffer_double[i] = my_rank * 0.1;
        read_darray_buffer_double[i] = 0.0;
    }

#if 1
    for (int fmt = 0; fmt < 2; fmt++) {
        int t0 = clock();
        /* Create a filename to write. */
        sprintf(filename, "example1_%d.nc", fmt);

        ret = PIOc_createfile(iosysid, &ncid_write, &(formats[fmt]), filename, PIO_CLOBBER);
        ERR

        /* Put some global attributes. */
        ret = PIOc_put_att(ncid_write, PIO_GLOBAL, "dummy_global_att_int", PIO_INT, 1, &put_global_att_int_data);
        ERR
        ret = PIOc_put_att(ncid_write, PIO_GLOBAL, "dummy_global_att_float", PIO_FLOAT, 1, &put_global_att_float_data);
        ERR
        ret = PIOc_put_att_text(ncid_write, PIO_GLOBAL, "dummy_global_att_text", strlen(put_global_att_text_data),
                                put_global_att_text_data);
        ERR

        /* Define some scalar variables. */
        ret = PIOc_def_var(ncid_write, "dummy_scalar_var_int", PIO_INT, 0, NULL, &varid_dummy_scalar_var_int);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_scalar_var_float", PIO_FLOAT, 0, NULL, &varid_dummy_scalar_var_float);
        ERR

        /* Define a text variable. */
        ret = PIOc_def_dim(ncid_write, "text_var_len", PIO_MAX_NAME, &dimid_text_var_len);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_text_var", PIO_CHAR, 1, &dimid_text_var_len, &varid_dummy_text_var);
        ERR

        /* Define some variables for PIOc_write_darray. */
        ret = PIOc_def_dim(ncid_write, "darray_var_len", (PIO_Offset) gdimlen[0], &dimid_darray_var_len);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_int", PIO_INT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_int);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float1", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float1);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float2", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float2);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float3", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float3);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float4", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float4);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float5", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float5);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float6", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float6);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float7", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float7);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float8", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float8);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_darray_var_float9", PIO_FLOAT, NDIMS, &dimid_darray_var_len,
                           &varid_dummy_darray_var_float9);
        ERR


        /* Put some local attributes for variable dummy_darray_var_int. */
        ret = PIOc_put_att(ncid_write, varid_dummy_darray_var_int, "dummy_att_float", PIO_FLOAT, 1,
                           &put_att_float_data);
        ERR
        ret = PIOc_put_att(ncid_write, varid_dummy_darray_var_int, "dummy_att_int", PIO_INT, 1, &put_att_int_data);
        ERR
        ret = PIOc_put_att_text(ncid_write, varid_dummy_darray_var_int, "dummy_att_text", strlen(put_att_text_data),
                                put_att_text_data);
        ERR

        /* Define some variables for PIOc_put_vars. */
        ret = PIOc_def_dim(ncid_write, "put_get_var_len", PUT_GET_VAR_LEN, &dimid_put_get_var_len);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_put_get_var_int", PIO_INT, NDIMS, &dimid_put_get_var_len,
                           &varid_dummy_put_get_var_int);
        ERR
        ret = PIOc_def_var(ncid_write, "dummy_put_get_var_float", PIO_FLOAT, NDIMS, &dimid_put_get_var_len,
                           &varid_dummy_put_get_var_float);
        ERR

        ret = PIOc_enddef(ncid_write);
        ERR

        /* Put some scalar variables. */
        ret = PIOc_put_var_int(ncid_write, varid_dummy_scalar_var_int, &put_scalar_int_data);
        ERR
        ret = PIOc_put_var_float(ncid_write, varid_dummy_scalar_var_float, &put_scalar_float_data);
        ERR

        /* Put a text variable. */
        ret = PIOc_put_var_text(ncid_write, varid_dummy_text_var, put_text_data);
        ERR

        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float1, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float2, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float3, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float4, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR

        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float5, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float6, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float7, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float8, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR
        /* Write to float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_write_darray(ncid_write, varid_dummy_darray_var_float9, ioid_double, ELEMENTS_PER_PE,
                                write_darray_buffer_double, NULL);
        ERR

        ret = PIOc_closefile(ncid_write);
        ERR
        if (my_rank == 0) printf("write done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
    }
#endif
    /* Read support is not implemented for ADIOS type yet: change to "fmt < 2" for testing this new feature later.
       Currently, ADIOS type in SCORPIO simply changes actual type to PnetCDF or NetCDF for reading. */
    for (int fmt = 0; fmt < 2; fmt++) {
        int t0 = clock();
        /* Note: for ADIOS type, the actual file name on disk is a BP directory named example1_1.nc.bp.dir,
           client code like E3SM still assumes example1_1.nc as the file name to be read.  */
        sprintf(filename, "example1_%d.nc", fmt);

        ret = PIOc_openfile(iosysid, &ncid_read, &(formats[fmt]), filename, PIO_NOWRITE);
        ERR

        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float", &varid_dummy_darray_var_float);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float1", &varid_dummy_darray_var_float1);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float2", &varid_dummy_darray_var_float2);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float3", &varid_dummy_darray_var_float3);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float4", &varid_dummy_darray_var_float4);
        ERR

        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float5", &varid_dummy_darray_var_float5);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float6", &varid_dummy_darray_var_float6);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float7", &varid_dummy_darray_var_float7);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float8", &varid_dummy_darray_var_float8);
        ERR
        /* Read float type variable with double type decomposition, type conversions will be performed. */
        ret = PIOc_inq_varid(ncid_read, "dummy_darray_var_float9", &varid_dummy_darray_var_float9);
        ERR

        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float1, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float2, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float3, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float4, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float5, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float6, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float7, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float8, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ERR
        ret = PIOc_read_darray(ncid_read, varid_dummy_darray_var_float9, ioid_double, ELEMENTS_PER_PE,
                               read_darray_buffer_double);
        ERR
        ret = PIOc_closefile(ncid_read);
        ERR
        if (my_rank == 0) printf("reading done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
    }

    ret = PIOc_freedecomp(iosysid, ioid_int);
    ERR
    ret = PIOc_freedecomp(iosysid, ioid_double);
    ERR

    ret = PIOc_finalize(iosysid);
    ERR

    MPI_Finalize();

#ifdef TIMING
    GPTLfinalize();
#endif

    return 0;
}
