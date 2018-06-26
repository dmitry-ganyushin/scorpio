/**
 * @file
 * PIO File Handling
 */
#include <config.h>
#include <pio.h>
#include <pio_internal.h>
#include "../../tools/adios2pio-nm/adios2pio-nm-lib-c.h"

/**
 * Open an existing file using PIO library.
 *
 * If the open fails, try again as netCDF serial before giving
 * up. Input parameters are read on comp task 0 and ignored elsewhere.
 *
 * Note that the file is opened with default fill mode, NOFILL for
 * pnetcdf, and FILL for netCDF classic and netCDF-4 files.
 *
 * @param iosysid : A defined pio system descriptor (input)
 * @param ncidp : A pio file descriptor (output)
 * @param iotype : A pio output format (input)
 * @param filename : The filename to open
 * @param mode : The netcdf mode for the open operation
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile
 * @author Jim Edwards, Ed Hartnett
 */
int PIOc_openfile(int iosysid, int *ncidp, int *iotype, const char *filename,
                  int mode)
{
    return openfile_int(iosysid, ncidp, iotype, filename, mode, 1);
}

/**
 * Open an existing file using PIO library.
 *
 * This is like PIOc_openfile(), but if the open fails, this function
 * will not try to open again as netCDF serial before giving
 * up. Input parameters are read on comp task 0 and ignored elsewhere.
 *
 * Note that the file is opened with default fill mode, NOFILL for
 * pnetcdf, and FILL for netCDF classic and netCDF-4 files.
 *
 * @param iosysid : A defined pio system descriptor (input)
 * @param ncidp : A pio file descriptor (output)
 * @param iotype : A pio output format (input)
 * @param filename : The filename to open
 * @param mode : The netcdf mode for the open operation
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile
 * @author Ed Hartnett
 */
int PIOc_openfile2(int iosysid, int *ncidp, int *iotype, const char *filename,
                   int mode)
{
    return openfile_int(iosysid, ncidp, iotype, filename, mode, 0);
}

/**
 * Open an existing file using PIO library.
 *
 * Input parameters are read on comp task 0 and ignored elsewhere.
 *
 * @param iosysid A defined pio system descriptor
 * @param path The filename to open
 * @param mode The netcdf mode for the open operation
 * @param ncidp pointer to int where ncid will go
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_openfile
 * @author Ed Hartnett
 */
int PIOc_open(int iosysid, const char *path, int mode, int *ncidp)
{
    int iotype;

    LOG((1, "PIOc_open iosysid = %d path = %s mode = %x", iosysid, path, mode));

    /* Set the default iotype. */
#ifdef _NETCDF
    iotype = PIO_IOTYPE_NETCDF;
#else /* Assume that _PNETCDF is defined. */
    iotype = PIO_IOTYPE_PNETCDF;
#endif

    /* Figure out the iotype. */
    if (mode & NC_NETCDF4)
    {
#ifdef _NETCDF4
        if (mode & NC_MPIIO || mode & NC_MPIPOSIX)
            iotype = PIO_IOTYPE_NETCDF4P;
        else
            iotype = PIO_IOTYPE_NETCDF4C;
#endif
    }
    else
    {
#ifdef _PNETCDF
        if (mode & NC_PNETCDF || mode & NC_MPIIO)
            iotype = PIO_IOTYPE_PNETCDF;
#endif
    }

    /* Open the file. If the open fails, do not retry as serial
     * netCDF. Just return the error code. */
    return PIOc_openfile_retry(iosysid, ncidp, &iotype, path, mode, 0);
}

/**
 * Create a new file using pio. Input parameters are read on comp task
 * 0 and ignored elsewhere. NOFILL mode will be turned on in all
 * cases.
 *
 * @param iosysid A defined pio system ID, obtained from
 * PIOc_InitIntercomm() or PIOc_InitAsync().
 * @param ncidp A pointer that gets the ncid of the newly created
 * file.
 * @param iotype A pointer to a pio output format. Must be one of
 * PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_IOTYPE_NETCDF4C, or
 * PIO_IOTYPE_NETCDF4P.
 * @param filename The filename to create.
 * @param mode The netcdf mode for the create operation.
 * @returns 0 for success, error code otherwise.
 * @ingroup PIO_createfile
 * @author Jim Edwards, Ed Hartnett
 */
int PIOc_createfile(int iosysid, int *ncidp, int *iotype, const char *filename,
                    int mode)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    int ret;               /* Return code from function calls. */

#ifdef TIMING
    GPTLstart("PIO:PIOc_createfile");
#endif

    /* Get the IO system info from the id. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
    {
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__,
                        "Unable to create file (%s, mode = %d, iotype=%s). Invalid arguments provided, invalid iosystem id (iosysid = %d)", (filename) ? filename : "NULL", mode, (!iotype) ? "UNKNOWN" : pio_iotype_to_string(*iotype), iosysid);
    }

    /* Create the file. */
    if ((ret = PIOc_createfile_int(iosysid, ncidp, iotype, filename, mode)))
    {
#ifdef TIMING
        GPTLstop("PIO:PIOc_createfile");
#endif
        return pio_err(ios, NULL, ret, __FILE__, __LINE__,
                        "Unable to create file (%s, mode = %d, iotype=%s) on iosystem (iosystem id = %d). Internal error creating the file", (filename) ? filename : "NULL", mode, (!iotype) ? "UNKNOWN" : pio_iotype_to_string(*iotype), iosysid);
    }

    /* Run this on all tasks if async is not in use, but only on
     * non-IO tasks if async is in use. (Because otherwise, in async
     * mode, set_fill would be called twice by each IO task, since
     * PIOc_createfile() will already be called on each IO task.) */
    if (!ios->async || !ios->ioproc)
    {
        /* Set the fill mode to NOFILL. */
        if ((ret = PIOc_set_fill(*ncidp, NC_NOFILL, NULL)))
        {
            return pio_err(ios, NULL, ret, __FILE__, __LINE__,
                            "Unable to create file (%s, mode = %d, iotype=%s) on iosystem (iosystem id = %d). Setting fill mode to NOFILL failed.", (filename) ? filename : "NULL", mode, (!iotype) ? "UNKNOWN" : pio_iotype_to_string(*iotype), iosysid);
        }
    }

#ifdef TIMING
    GPTLstop("PIO:PIOc_createfile");
#endif
    return ret;
}

/**
 * Open a new file using pio. The default fill mode will be used (FILL
 * for netCDF and netCDF-4 formats, NOFILL for pnetcdf.) Input
 * parameters are read on comp task 0 and ignored elsewhere.
 *
 * @param iosysid : A defined pio system descriptor (input)
 * @param cmode : The netcdf mode for the create operation.
 * @param filename : The filename to open
 * @param ncidp : A pio file descriptor (output)
 * @return 0 for success, error code otherwise.
 * @ingroup PIO_create
 * @author Ed Hartnett
 */
int PIOc_create(int iosysid, const char *filename, int cmode, int *ncidp)
{
    int iotype;            /* The PIO IO type. */

    /* Set the default iotype. */
#ifdef _NETCDF
    iotype = PIO_IOTYPE_NETCDF;
#else /* Assume that _PNETCDF is defined. */
    iotype = PIO_IOTYPE_PNETCDF;
#endif

    /* Figure out the iotype. */
    if (cmode & NC_NETCDF4)
    {
#ifdef _NETCDF4
        if (cmode & NC_MPIIO || cmode & NC_MPIPOSIX)
            iotype = PIO_IOTYPE_NETCDF4P;
        else
            iotype = PIO_IOTYPE_NETCDF4C;
#endif
    }
    else
    {
#ifdef _PNETCDF
        if (cmode & NC_PNETCDF || cmode & NC_MPIIO)
            iotype = PIO_IOTYPE_PNETCDF;
#endif
    }

    return PIOc_createfile_int(iosysid, ncidp, &iotype, filename, cmode);
}

/* Internal helper function to perform sync operations
 * ncid : the ncid of the file to sync
 * Returns PIO_NOERR for success, error code otherwise
 */
static int sync_file(int ncid)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    bool flushtodisk = false;

    LOG((1, "sync_file ncid = %d", ncid));

    /* Get the file info from the ncid. */
    if ((ierr = pio_get_file(ncid, &file)))
    {
        return pio_err(NULL, NULL, ierr, __FILE__, __LINE__,
                        "Syncing file (ncid=%d) failed. Invalid file id. Unable to find internal structure associated with the file id", ncid);
    }

#ifdef _ADIOS
    if (file->iotype == PIO_IOTYPE_ADIOS)
        return PIO_NOERR;
#endif

    ios = file->iosystem;

    /* Flush data buffers on computational tasks. */
    if (!ios->async || !ios->ioproc)
    {
        if (file->mode & PIO_WRITE)
        {
            wmulti_buffer *wmb, *twmb;

            LOG((3, "sync_file checking buffers"));
            wmb = &file->buffer;
            while (wmb)
            {
                /* If there are any data arrays waiting in the
                 * multibuffer, flush it to IO tasks. */
                if (wmb->num_arrays > 0)
                    flush_buffer(ncid, wmb, flushtodisk);
                twmb = wmb;
                wmb = wmb->next;
                if (twmb == &file->buffer)
                {
                    twmb->ioid = -1;
                    twmb->next = NULL;
                }
                else
                {
                    free(twmb);
                }
            }
        }
    }

    /* If async is in use, send message to IO master tasks. */
    if (ios->async)
    {
        int msg = PIO_MSG_SYNC;

        PIO_SEND_ASYNC_MSG(ios, msg, &ierr, ncid);
        if (ierr != PIO_NOERR)
        {
            return pio_err(ios, NULL, ierr, __FILE__, __LINE__,
                            "Syncing file %s (ncid=%d) failed. Unable to send asynchronous message, PIO_MSG_SYNC, on iosystem (iosysid=%d)", pio_get_fname_from_file(file), ncid, ios->iosysid);
        }
    }

    /* Call the sync function on IO tasks.

       We choose not to call ncmpi_sync() for PIO_IOTYPE_PNETCDF,
       as it has been confirmed to have a very high cost on some
       systems. ncmpi_sync() itself does nothing but simply calls
       MPI_File_sync(), which usually incurs a huge performance
       penalty by calling POSIX sync internally. It is designed
       to ensure the data is safely stored on the disk hardware,
       before the function returns. People use it for extremely
       cautious behavior only.
     */
    if (file->mode & PIO_WRITE)
    {
        if (ios->ioproc)
        {
            switch (file->iotype)
            {
#ifdef _NETCDF4
            case PIO_IOTYPE_NETCDF4P:
                ierr = nc_sync(file->fh);
                break;
            case PIO_IOTYPE_NETCDF4C:
#endif
#ifdef _NETCDF
            case PIO_IOTYPE_NETCDF:
                if (ios->io_rank == 0)
                    ierr = nc_sync(file->fh);
                break;
#endif
#ifdef _PNETCDF
            case PIO_IOTYPE_PNETCDF:
                flush_output_buffer(file, true, 0);
                break;
#endif
            default:
                return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__,
                                "Syncing file %s (ncid=%d) failed. Invalid/Unsupported iotype (%s:%d) provided", pio_get_fname_from_file(file), ncid, pio_iotype_to_string(file->iotype), file->iotype);
            }
        }
        LOG((2, "sync_file ierr = %d", ierr));
    }

    ierr = check_netcdf(ios, NULL, ierr, __FILE__, __LINE__);
    if (ierr != PIO_NOERR)
    {
        LOG((1, "nc*_sync failed, ierr = %d", ierr));
        return ierr;
    }

    return PIO_NOERR;
}

/* Close the file ("hard close")
 * @param ios: Pointer to the iosystem_desc
 * @param file: Pointer to the file_desc for the file
 * @returns PIO_NOERR for success, a pio error code otherwise
 */
int PIO_hard_closefile(iosystem_desc_t *ios, file_desc_t *file)
{
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS;
#ifdef _ADIOS
    char outfilename[PIO_MAX_NAME + 1];
    size_t len = 0;
#endif

    assert(ios && file);

    /* ADIOS: assume all procs are also IO tasks */
#ifdef _ADIOS
    if (file->iotype == PIO_IOTYPE_ADIOS)
    {
        if (file->adios_fh != -1)
        {
            LOG((2, "ADIOS close file %s", file->filename));
            adios_define_attribute_byvalue(file->adios_group, "/__pio__/fillmode", "",
                                           adios_integer, 1, &file->fillmode);
            ierr = adios_close(file->adios_fh);
            file->adios_fh = -1;
        }

        if (file->adios_group != -1)
        {
            adios_free_group(file->adios_group);
            file->adios_group = -1;
        }

        for (int i = 0; i < file->num_dim_vars; i++)
        {
            free(file->dim_names[i]);
            file->dim_names[i] = NULL;
        }

        file->num_dim_vars = 0;

        for (int i = 0; i < file->num_vars; i++)
        {
            free(file->adios_vars[i].name);
            file->adios_vars[i].name = NULL;
            free(file->adios_vars[i].gdimids);
            file->adios_vars[i].gdimids = NULL;
        }

        file->num_vars = 0;

        /* Track attributes */
        for (int i = 0; i < file->num_attrs; i++)
        {
            free(file->adios_attrs[i].att_name);
            file->adios_attrs[i].att_name = NULL;
        }

        file->num_attrs = 0;

#ifdef _ADIOS_BP2NC_TEST /* Comment out for large scale run */
#ifdef _PNETCDF
        char conv_iotype[] = "pnetcdf";
#else
        char conv_iotype[] = "netcdf";
#endif

        /* Convert XXXX.nc.bp to XXXX.nc */
        len = strlen(file->filename);
        assert(len > 6 && len <= PIO_MAX_NAME);
        strncpy(outfilename, file->filename, len - 3);
        outfilename[len - 3] = '\0';
        LOG((1, "CONVERTING: %s", file->filename));
        C_API_ConvertBPToNC(file->filename, outfilename, conv_iotype, ios->union_comm);
        LOG((1, "DONE CONVERTING: %s", file->filename));
#endif

        free(file->filename);
        ierr = 0;
    }
#endif

#if PIO_ENABLE_SOFT_SYNC
    /* Hard sync - the previous syncs were "soft sync"s */
    ierr = PIO_hard_sync(ios, file);
    if(ierr != PIO_NOERR)
    {
        return pio_err(ios, file, ierr, __FILE__, __LINE__,
                        "Closing file (%s, ncid=%d) failed. A hard sync on the file failed (iotype=%s)", pio_get_fname_from_file(file), file->pio_ncid, pio_iotype_to_string(file->iotype));
    }
#endif

    /* If this is an IO task, then call the netCDF function. */
    if (ios->ioproc)
    {
        switch (file->iotype)
        {
#ifdef _NETCDF4
        case PIO_IOTYPE_NETCDF4P:
            ierr = nc_close(file->fh);
            break;
        case PIO_IOTYPE_NETCDF4C:
#endif /* _NETCDF4 */
        case PIO_IOTYPE_NETCDF:
            if (ios->io_rank == 0)
                ierr = nc_close(file->fh);
            break;
#ifdef _PNETCDF
        case PIO_IOTYPE_PNETCDF:
            if ((file->mode & PIO_WRITE)){
                ierr = ncmpi_buffer_detach(file->fh);
            }
            ierr = ncmpi_close(file->fh);
            break;
#endif /* _PNETCDF */
#ifdef _ADIOS
        case PIO_IOTYPE_ADIOS: /* Needed to avoid default case and error. */
            ierr = 0;
            break;
#endif
        default:
            return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__,
                            "Closing file (%s, ncid=%d) failed. Unsupported iotype (%d) specified", pio_get_fname_from_file(file), file->pio_ncid, file->iotype);
        }
    }

    ierr = check_netcdf(NULL, file, ierr, __FILE__, __LINE__);
    if(ierr != PIO_NOERR){
        LOG((1, "nc*_close failed, ierr = %d", ierr));
        return pio_err(NULL, file, ierr, __FILE__, __LINE__,
                        "Closing file (%s, ncid=%d) failed. Underlying I/O library (iotype=%s) call failed", pio_get_fname_from_file(file), file->pio_ncid, pio_iotype_to_string(file->iotype));
    }

    /* Delete file from our list of open files. */
    ierr = pio_delete_file_from_list(file->pio_ncid);
    if(ierr != PIO_NOERR)
    {
        return pio_err(ios, file, ierr, __FILE__, __LINE__,
                        "Closing file (%s, ncid=%d) failed. Deleting file from the list of open files failed", pio_get_fname_from_file(file), file->pio_ncid);
    }

    return PIO_NOERR;
}

/* "Soft close" the file
 * The function assumes that only writes are pending on this file
 * @param ios: Pointer to the iosystem_desc
 * @param file: Pointer to the file_desc for the file
 * @returns PIO_NOERR for success, a pio error code otherwise
 */
int PIO_soft_closefile(iosystem_desc_t *ios, file_desc_t *file)
{
    int ierr = PIO_NOERR;

    assert(ios && file);
    if(file->nasync_pend_ops > 0){
#if PIO_USE_ASYNC_WR_THREAD
        /* Queue up the pending async ops on the file to the async thread pool */
        ierr = pio_tpool_async_pend_op_add(ios, PIO_ASYNC_FILE_WRITE_OPS,
                (void *)file);
        if(ierr != PIO_NOERR)
        {
            return pio_err(ios, file, ierr, __FILE__, __LINE__,
                            "Closing file (%s, ncid=%d) failed. Soft close of the file failed (iosysid=%d). Internal error while adding pending asynchronous write operations for this file into the thread pool", pio_get_fname_from_file(file), file->pio_ncid, ios->iosysid);
        }
#else
        /* Queue up the pending async ops on the file to the iosystem */
        ierr = pio_iosys_async_pend_op_add(ios, PIO_ASYNC_FILE_WRITE_OPS,
                (void *)file);
        if(ierr != PIO_NOERR)
        {
            return pio_err(ios, file, ierr, __FILE__, __LINE__,
                            "Closing file (%s, ncid=%d) failed. Soft close of the file failed (iosysid=%d). Internal error while adding pending asynchronous write operations for this file into the iosystem queue", pio_get_fname_from_file(file), file->pio_ncid, ios->iosysid);
        }
#endif
    }
    else{
        ierr = PIO_hard_closefile(ios, file);
        if(ierr != PIO_NOERR){
            return pio_err(ios, file, ierr, __FILE__, __LINE__,
                            "Closing file (%s, ncid=%d) failed. Internal error while closing the file (iosysid=%d).", pio_get_fname_from_file(file), file->pio_ncid, ios->iosysid);
        }
    }

    return PIO_NOERR;
}

/**
 * Close a file previously opened with PIO.
 *
 * @param ncid: the file pointer
 * @returns PIO_NOERR for success, error code otherwise.
 * @author Jim Edwards, Ed Hartnett
 */
int PIOc_closefile(int ncid)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    file_desc_t *file;     /* Pointer to file information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS, mpierr2;  /* Return code from MPI function codes. */
#ifdef _ADIOS
    char outfilename[PIO_MAX_NAME + 1];
    size_t len = 0;
#endif

#ifdef TIMING
    GPTLstart("PIO:PIOc_closefile");
#endif
    LOG((1, "PIOc_closefile ncid = %d", ncid));

    /* Find the info about this file. */
    if ((ierr = pio_get_file(ncid, &file)))
    {
        return pio_err(NULL, NULL, ierr, __FILE__, __LINE__,
                        "Closing file (ncid=%d) failed. Invalid file id. Unable to find internal structure associated with the file id", ncid);
    }
    ios = file->iosystem;

#ifdef TIMING
    if (file->mode & PIO_WRITE)
        GPTLstart("PIO:PIOc_closefile_write_mode");
#endif

    /* Sync changes before closing on all tasks if async is not in
     * use, but only on non-IO tasks if async is in use. */
    if (!ios->async || !ios->ioproc)
        if (file->mode & PIO_WRITE)
            sync_file(ncid);

    /* If async is in use and this is a comp tasks, then the compmaster
     * sends a msg to the pio_msg_handler running on the IO master and
     * waiting for a message. Then broadcast the ncid over the intercomm
     * to the IO tasks. */
    if (ios->async)
    {
        int msg = PIO_MSG_CLOSE_FILE;

        PIO_SEND_ASYNC_MSG(ios, msg, &ierr, ncid);
        if(ierr != PIO_NOERR)
        {
            return pio_err(ios, file, ierr, __FILE__, __LINE__,
                        "Closing file %s (ncid=%d) failed. Unable to send asynchronous message, PIO_MSG_CLOSE_FILE, on iosystem (iosysid=%d)", pio_get_fname_from_file(file), ncid, ios->iosysid);
        }
    }

#if PIO_ENABLE_SOFT_CLOSE
    /* A "soft close" does not sync all data to the disk */
    ierr = PIO_soft_closefile(ios, file);
    if(ierr != PIO_NOERR)
    {
        return pio_err(ios, file, ierr, __FILE__, __LINE__,
                        "Closing file %s (ncid=%d) failed. Internal error while performing soft close (asynchronous close) on the file", pio_get_fname_from_file(file), ncid);
    }

#else /* PIO_ENABLE_SOFT_CLOSE */
    ierr = PIO_hard_closefile(ios, file);
    if(ierr != PIO_NOERR)
    {
        return pio_err(ios, file, ierr, __FILE__, __LINE__,
                        "Closing file %s (ncid=%d) failed. Internal error while  closing the file", pio_get_fname_from_file(file), ncid);
    }
#endif /* PIO_ENABLE_SOFT_CLOSE */

#ifdef TIMING
    if (file->mode & PIO_WRITE)
        GPTLstop("PIO:PIOc_closefile_write_mode");
#endif

#ifdef TIMING
    GPTLstop("PIO:PIOc_closefile");
#endif
    return ierr;
}

/**
 * Delete a file.
 *
 * @param iosysid a pio system handle.
 * @param filename a filename.
 * @returns PIO_NOERR for success, error code otherwise.
 * @author Jim Edwards, Ed Hartnett
 */
int PIOc_deletefile(int iosysid, const char *filename)
{
    iosystem_desc_t *ios;  /* Pointer to io system information. */
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    int mpierr = MPI_SUCCESS;  /* Return code from MPI function codes. */
     int msg = PIO_MSG_DELETE_FILE;
    size_t len;

#ifdef TIMING
    GPTLstart("PIO:PIOc_deletefile");
#endif
    LOG((1, "PIOc_deletefile iosysid = %d filename = %s", iosysid, filename));

    /* Get the IO system info from the id. */
    if (!(ios = pio_get_iosystem_from_id(iosysid)))
    {
        return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__,
                        "Deleting file (%s) failed. Invalid I/O system id (iosysid=%d) specified.", (filename) ? filename : "NULL", iosysid);
    }

    /* If async is in use, send message to IO master task. */
    if (ios->async)
    {
        len = strlen(filename) + 1;

        PIO_SEND_ASYNC_MSG(ios, msg, &ierr, len, filename);
        if(ierr != PIO_NOERR)
        {
            return pio_err(ios, NULL, ierr, __FILE__, __LINE__,
                        "Deleting file (%s) failed. Sending async message, PIO_MSG_DELETE_FILE, failed", (filename) ? filename : "NULL");
        }
    }

    /* If this is an IO task, then call the netCDF function. The
     * barriers are needed to assure that no task is trying to operate
     * on the file while it is being deleted. IOTYPE is not known, but
     * nc_delete() will delete any type of file. */
    if (ios->ioproc)
    {
        mpierr = MPI_Barrier(ios->io_comm);

        if (!mpierr && ios->io_rank == 0)
#ifdef _NETCDF
             ierr = nc_delete(filename);
#else /* Assume that _PNETCDF is defined. */
             ierr = ncmpi_delete(filename, MPI_INFO_NULL);
#endif

        if (!mpierr)
            mpierr = MPI_Barrier(ios->io_comm);
    }
    LOG((2, "PIOc_deletefile ierr = %d", ierr));

    ierr = check_netcdf(ios, NULL, ierr, __FILE__, __LINE__);
    if(ierr != PIO_NOERR){
        return pio_err(ios, NULL, ierr, __FILE__, __LINE__,
                    "Deleting file (%s) failed. Internal I/O library call failed.", (filename) ? filename : "NULL");
    }

#ifdef TIMING
    GPTLstop("PIO:PIOc_deletefile");
#endif
    return ierr;
}

/** Sync data on this file to the disk
 * @param ios Pointer to the iosystem_desc
 * @param file Pointer to the file_desc
 * @returns PIO_NOERR on success, a pio error code otherwise
 */
int PIO_hard_sync(iosystem_desc_t *ios, file_desc_t *file)
{
    int ierr = PIO_NOERR;
    int mpierr = MPI_SUCCESS, mpierr2;

    assert(ios && file);

    /* Call the sync function on IO tasks. */
    if (file->mode & PIO_WRITE)
    {
        if (ios->ioproc)
        {
            switch(file->iotype)
            {
#ifdef _NETCDF4
            case PIO_IOTYPE_NETCDF4P:
                ierr = nc_sync(file->fh);
                break;
            case PIO_IOTYPE_NETCDF4C:
#endif
            case PIO_IOTYPE_NETCDF:
                if (ios->io_rank == 0)
                    ierr = nc_sync(file->fh);
                break;
#ifdef _PNETCDF
            case PIO_IOTYPE_PNETCDF:
                /* force flush data */
                flush_output_buffer(file, true, 0);
                ierr = ncmpi_sync(file->fh);
                break;
#endif
            default:
                return pio_err(ios, file, PIO_EBADIOTYPE, __FILE__, __LINE__,
                                "Syncing file (%s, ncid=%d) failed. Unsupported iotype (%d) specified", pio_get_fname_from_file(file), file->pio_ncid, file->iotype);
            }
        }
        LOG((2, "PIOc_sync ierr = %d", ierr));
    }

    /* Broadcast and check the return code. */
    if ((mpierr = MPI_Bcast(&ierr, 1, MPI_INT, ios->ioroot, ios->my_comm)))
        return check_mpi(ios, NULL, mpierr, __FILE__, __LINE__);

    return ierr;
}

/**
 * PIO interface to nc_sync This routine is called collectively by all
 * tasks in the communicator ios.union_comm.
 *
 * Refer to the <A
 * HREF="http://www.unidata.ucar.edu/software/netcdf/docs/modules.html"
 * target="_blank"> netcdf </A> documentation.
 *
 * @param ncid the ncid of the file to sync.
 * @returns PIO_NOERR for success, error code otherwise.
 * @author Jim Edwards, Ed Hartnett
 */
int PIOc_sync(int ncid)
{
    int ierr = PIO_NOERR;  /* Return code from function calls. */
    bool flushtodisk = true;

#ifdef TIMING
    GPTLstart("PIO:PIOc_sync");
#endif

    LOG((1, "PIOc_sync ncid = %d", ncid));

    ierr = sync_file(ncid);

#ifdef TIMING
    GPTLstop("PIO:PIOc_sync");
#endif
    return ierr;
}
