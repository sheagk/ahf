// For simplicity, just feed AHF the high-res DM particles right now
#include "io_hdf5.h"
#include "hdf5.h"
#include "io.h"
#include "io_gadget_def.h"

// This assumes that the data is contiguious, if not we get UB
void *HDF5_READ(hid_t handle, const char *path) {
    hid_t dset = H5Dopen(handle, path, H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    hsize_t size = H5Sget_simple_extent_npoints(space);
    hid_t dset_type = H5Dget_type(dset);
    size_t type_size = H5Tget_size(dset_type);

    printf("%llu %zu\n", size, type_size);

    void *buf = malloc(size * type_size);
    H5Dread(dset, dset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    H5Tclose(dset_type);
    H5Sclose(space);
    H5Dclose(dset);

    return buf;
}

extern io_hdf5_t io_hdf5_open(io_logging_t log, char *fname,
                              io_file_swap_t swapped, io_file_mode_t mode,
                              uint32_t reader) {

    if (mode == IO_FILE_WRITE) {
        io_logging_fatal(log, "HDF5 writing is not yet supported");
        return NULL;
    }

    io_hdf5_t f = malloc(sizeof(io_hdf5_t));
    f->ftype = IO_FILE_HDF5;
    hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    f->coordinates = (float *)HDF5_READ(file, "PartType1/Coordinates");
    f->velocities = (float *)HDF5_READ(file, "PartType1/Velocities");
    f->masses = (float *)HDF5_READ(file, "PartType1/Masses");
    f->particleids = (unsigned int *)HDF5_READ(file, "PartType1/ParticleIDs");

    // TODO: Record the max and min values for params to do scaling

    hid_t dset = H5Gopen(file, "Header", H5P_DEFAULT);
    hid_t attr;

    attr = H5Aopen(dset, "NumPart_ThisFile", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &f->header.np);
    H5Aclose(attr);

    attr = H5Aopen(dset, "MassTable", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.massarr);
    H5Aclose(attr);

    attr = H5Aopen(dset, "HubbleParam", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.expansion);
    H5Aclose(attr);

    attr = H5Aopen(dset, "Time", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.expansion);
    H5Aclose(attr);

    attr = H5Aopen(dset, "Redshift", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.redshift);
    H5Aclose(attr);

    attr = H5Aopen(dset, "BoxSize", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.boxsize);
    H5Aclose(attr);

    attr = H5Aopen(dset, "Omega0", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.omega0);
    H5Aclose(attr);

    attr = H5Aopen(dset, "OmegaLambda", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &f->header.omegalambda);
    H5Aclose(attr);

    H5Fclose(file);
    return f;
}

extern void io_hdf5_init(io_logging_t log, io_hdf5_t f) {}

extern bool io_hdf5_get(io_logging_t log, io_hdf5_t f, io_file_get_t what,
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

extern bool io_hdf5_set(io_logging_t log, io_hdf5_t f, io_file_get_t what,
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

extern uint64_t io_hdf5_readpart(io_logging_t log, io_hdf5_t f, uint64_t pskip,
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

extern void io_hdf5_close(io_logging_t log, io_hdf5_t *f) {
    free((*f)->coordinates);
    free((*f)->velocities);
    free((*f)->masses);
    free((*f)->particleids);
    free(*f);
    f = NULL;
}

extern void io_hdf5_log(io_logging_t log, io_hdf5_t f) {}
