// For simplicity, just feed AHF the high-res DM particles right now
#include "io_hdf5.h"
#include "hdf5.h"
#include "io.h"
#include "io_gadget_def.h"

extern io_mhdf5_t io_mhdf5_open(io_logging_t log, char *fname,
                              io_file_swap_t swapped, io_file_mode_t mode,
                              uint32_t reader) {

    if (mode == IO_FILE_WRITE) {
        io_logging_fatal(log, "HDF5 writing is not yet supported");
        return NULL;
    }

    io_mhdf5_t f = malloc(sizeof(io_mhdf5_t));
    f->ftype = IO_FILE_MHDF5;

    // Read multiple blocks -- note that I'm going to do this in parallel
    
    /* Split the filename in path and stem */
    f->path = NULL;
    f->stem = NULL;
    if (  io_util_split_pathfname(fname, &(f->path), &(f->stem)) == NULL) {
        io_logging_fatal(log, "io_mhdf5_open(): Could not split %s in path and filename.", fname);
    if(f->path) free(f->path);
    if(f->stem) free(f->stem);
        free(f);
        return NULL;
    }
    io_logging_msg(log, INT32_C(1), "Will look in %s for %s", f->path, f->stem);
    /* Get the filenames */
    f->numfiles = io_util_findfiles(f->path, f->stem, "%i", "", &fnames);
    
    if (f->numfiles <= 0) {
      io_logging_fatal(log, "io_mhdf5_open(): Could not open anything starting with %s in %s.", f->stem, f->path);
      if(f->path)  free(f->stem);
      if(f->stem)  free(f->path);
      free(f);
      return NULL;
    }

    /* Glue the files into the MHDF5 structure */
    f->files = (io_hdf5_t *)malloc( sizeof(io_hdf5_t)*(f->numfiles));
    if (f->files == NULL) {
        io_logging_memfatal(log,  "io_mhdf5 structure (2)");
        for (i=0; i<f->numfiles; i++)
            free(fnames[i]);
        free(fnames);
        free(f->stem);
        free(f->path);
        free(f);
        return NULL;
    }

    for (i=0; i<f->numfiles; i++) {
        io_logging_msg(log,6,"\ntrying to open file %s (%d of %d) for reading ... ", fnames[i],i+1,f->numfiles);
        (f->files)[i] = io_hdf5_open(log, fnames[i], swapped, mode, reader);
        if ((f->files)[i] == NULL) {
            int32_t j;
            for (j=i; i<f->numfiles; j++)
                free(fnames[j]);
            free(fnames);
            while (i>0) {
                i--;
                io_gadget_close(log, &((f->files)[i]));
            }
          free(f->stem);
          free(f->path);
          free(f);
          return NULL;
        }

    }
    /* the whole char *fnames[] is not needed anymore */
    for (i=0; i<f->numfiles; i++) {
    free(fnames[i]);
    }
    free(fnames);
    return f;
}

extern void io_mhdf5_init(io_logging_t log, io_hdf5_t f) {}

extern bool io_mhdf5_get(io_logging_t log, io_hdf5_t f, io_file_get_t what,
                        void *res) {
    switch (what) {
    case IO_FILE_GET_NOPART:
    case IO_FILE_GET_NOPART_IN_FILE:
        *(long *)res = (long)f->header.np[1];
        break;
    case IO_FILE_GET_NOVPART:
        break;
    case IO_FILE_GET_NOSPECIES:
        *(int *)res = 1;
        break;
    case IO_FILE_GET_BOXSIZE:
        *(double *)res = f->header.boxsize * f->posscale;
        break;
    case IO_FILE_GET_PMASS:
        *(double *)res = f->header.massarr[1] * f->weightscale;
        break;
    case IO_FILE_GET_ZINITIAL:
    case IO_FILE_GET_Z:
        *(double *)res = f->header.redshift;
        break;
    case IO_FILE_GET_AINITIAL:
    case IO_FILE_GET_A:
        *(double *)res = f->header.expansion;
        break;
    case IO_FILE_GET_OMEGA0:
        *(double *)res = f->header.omega0;
        break;
    case IO_FILE_GET_OMEGAL:
        *(double *)res = f->header.omegalambda;
        break;
    case IO_FILE_GET_H:
        *(double *)res = f->header.hubbleparameter;
        break;
    case IO_FILE_GET_DOUBLE:
        *(int *)res = 1;
        break;
    case IO_FILE_GET_MMASS:
        *(int *)res = 0;
        break;
    case IO_FILE_GET_NOTSTEP:
        *(int32_t *)res = 0;
        break;
    case IO_FILE_GET_TSTEP:
        *(double *)res = 0.0;
        break;
    case IO_FILE_GET_HEADERSTR:
        *(char **)res = "No header string.";
        break;
    case IO_FILE_GET_MINWEIGHT:
        *(double *)res = f->minweight / f->mmass;
        break;
    case IO_FILE_GET_MAXWEIGHT:
        *(double *)res = f->maxweight / f->mmass;
        break;
    }

    return true;
}

extern bool io_mhdf5_set(io_logging_t log, io_hdf5_t f, io_file_get_t what,
                        void *res) {
    if (f == NULL)
        return false;

    switch (what) {
    case IO_FILE_GET_BOXSIZE:
        f->header.boxsize = *((double *)res);
        break;
    case IO_FILE_GET_PMASS:
        f->header.massarr[1] = *((double *)res);
        break;
    case IO_FILE_GET_Z:
        f->header.redshift = *((double *)res);
        break;
    case IO_FILE_GET_A:
        f->header.expansion = *((double *)res);
        break;
    case IO_FILE_GET_OMEGA0:
        f->header.omega0 = *((double *)res);
        break;
    case IO_FILE_GET_OMEGAL:
        f->header.omegalambda = *((double *)res);
        break;
    case IO_FILE_GET_H:
        f->header.hubbleparameter = *((double *)res);
        break;
    default:
        io_logging_fatal(log, "Requesting something unkown in %s.", __func__);
        return false;
    }

    return true;
}

extern uint64_t io_mhdf5_readpart(io_logging_t log, io_hdf5_t f, uint64_t pskip,
                                 uint64_t pread, io_file_strg_struct_t strg) {

    for (int i = pskip; i < (pskip + pread); ++i) {
        *(float *)strg.posx.val = f->coordinates[3 * i];
        strg.posx.val += strg.posx.stride;
        if ((i == 0) || (f->coordinates[3 * i] < f->minpos[0]))
            f->minpos[0] = f->coordinates[3 * i];
        if ((i == 0) || (f->coordinates[3 * i] > f->maxpos[0]))
            f->maxpos[0] = f->coordinates[3 * i];

        *(float *)strg.posy.val = f->coordinates[3 * i + 1];
        strg.posx.val += strg.posx.stride;
        if ((i == 0) || (f->coordinates[3 * i + 1] < f->minpos[1]))
            f->minpos[1] = f->coordinates[3 * i + 1];
        if ((i == 0) || (f->coordinates[3 * i + 1] > f->maxpos[1]))
            f->maxpos[1] = f->coordinates[3 * i + 1];

        *(float *)strg.posz.val = f->coordinates[3 * i + 2];
        strg.posx.val += strg.posx.stride;
        if ((i == 0) || (f->coordinates[3 * i + 2] < f->minpos[2]))
            f->minpos[2] = f->coordinates[3 * i + 2];
        if ((i == 0) || (f->coordinates[3 * i + 2] > f->maxpos[2]))
            f->maxpos[2] = f->coordinates[3 * i + 2];

        *(float *)strg.momx.val = f->velocities[3 * i];
        strg.momx.val += strg.posx.stride;

        *(float *)strg.momy.val = f->velocities[3 * i + 1];
        strg.momy.val += strg.posx.stride;

        *(float *)strg.momz.val = f->velocities[3 * i + 2];
        strg.momz.val += strg.posz.stride;

        if (strg.weight.val != NULL) {
            *(float *)strg.weight.val = f->masses[i];
            strg.weight.val += strg.weight.stride;
            if ((i == 0) || (f->masses[i] < f->minweight))
                f->minweight = f->masses[i];
            if ((i == 0) || (f->masses[i] > f->maxweight))
                f->maxweight = f->masses[i];
        }

        *(int *)strg.id.val = f->particleids[i];
        strg.id.val += strg.id.stride;

        *(float *)strg.u.val = -1.0;
        strg.u.val += strg.u.stride;
    }

    return pread;
}

extern void io_mhdf5_close(io_logging_t log, io_hdf5_t *f) {
    free((*f)->coordinates);
    free((*f)->velocities);
    free((*f)->masses);
    free((*f)->particleids);
    free(*f);
    f = NULL;
}

extern void io_mhdf5_log(io_logging_t log, io_hdf5_t f) {}
