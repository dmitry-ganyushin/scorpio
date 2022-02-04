//
// Created by ganyush on 1/31/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <adios2_c.h>
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    adios2_adios *adios = adios2_init(MPI_COMM_WORLD, adios2_debug_mode_on);
    adios2_io *io = adios2_declare_io(adios, "Read");
    adios2_error adiosErr = adios2_set_engine(io, "FileStream");
    adios2_engine *e = adios2_open(io, "localArray.bp", adios2_mode_read);
    size_t step=0;
    adios2_steps(&step, e);
    adios2_step_status status;
    while (adios2_begin_step(e, adios2_step_mode_read, -1.,
                             &status) == adios2_error_none) {
        if (step == 0 || status == adios2_step_status_end_of_stream) {
            break;
        } else {
            adios2_end_step(e);
            step++;
            continue;
        }
    }
    size_t current_step;
    adios2_error err = adios2_current_step(&current_step, e);

    printf("step = %lu\n", step);
    printf("current step = %lu\n", current_step);
    adios2_variable *v0 =  adios2_inquire_variable(io, "v0");
    adios2_varinfo *data_blocks = adios2_inquire_blockinfo(e, v0,
                                                           current_step);
    int32_t number_of_data_blocks = data_blocks->nblocks;
    /* free memeory */
    for (size_t i = 0; i < data_blocks->nblocks; ++i) {
        free(data_blocks->BlocksInfo[i].Start);
        free(data_blocks->BlocksInfo[i].Count);
    }
    free(data_blocks->BlocksInfo);
    free(data_blocks);
    size_t var_size;
    for (size_t i = 0; i < number_of_data_blocks; i++) {
        adios2_set_block_selection(v0, i);
        adiosErr = adios2_selection_size(&var_size, v0);
        char *mem_buffer = (char *) calloc(var_size, 1);
        int t0 = clock();
        adios2_get(e, v0, mem_buffer, adios2_mode_sync);
        printf("done in %f\n", 1.0 *(clock() - t0) / CLOCKS_PER_SEC); t0 = clock();
        free(mem_buffer);
    }

    adios2_end_step(e);

    MPI_Finalize();
    return 0;
}

