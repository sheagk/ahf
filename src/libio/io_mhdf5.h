#ifndef IO_MHDF5_H
#define IO_MHDF5_H

// For simplicity, just feed AHF the high-res DM particles right now
#include "io.h"
#include "io_mgadget_def.h"

typedef struct _io_mhdf5 {
    io_file_type_t ftype;
    io_gadget_header_struct_t header;
    float *coordinates;
    float *velocities;
    float *masses;
    int *particleids;
    double minpos[3];
    double maxpos[3];
    double posscale;
    double weightscale;
    double minweight;
    double maxweight;
    double mmass;
} io_mhdf5;

typedef io_mhdf5 *io_mhdf5_t;

extern io_mhdf5_t io_mhdf5_open(io_logging_t log, char *fname,
                              io_file_swap_t swapped, io_file_mode_t mode,
                              uint32_t reader);

extern void io_mhdf5_init(io_logging_t log, io_mhdf5_t f);

extern bool io_mhdf5_get(io_logging_t log, io_mhdf5_t f, io_file_get_t what,
                        void *res);

extern bool io_mhdf5_set(io_logging_t log, io_mhdf5_t f, io_file_get_t what,
                        void *res);

extern uint64_t io_mhdf5_readpart(io_logging_t log, io_mhdf5_t f, uint64_t pskip,
                                 uint64_t pread, io_file_strg_struct_t strg);

extern void io_mhdf5_close(io_logging_t log, io_mhdf5_t *f);

extern void io_mhdf5_log(io_logging_t log, io_mhdf5_t f);

#endif
