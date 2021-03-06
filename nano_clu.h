#ifndef _NANO_CLU_H_
#define _NANO_CLU_H_

#include "debwt.h"
#include "bntseq.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#ifdef __DEBUG__
#define CHUNK_READ_N 1
#else
#define CHUNK_READ_N 100000
#endif
#define CHUNK_SIZE  10000000

#define NANO_MEM_LEN 14
#define NANO_UNI_OCC_THD 2

typedef struct {
    int n, m;
    int *vote_id, *hit, *vote_score;
} vote_t; // cluster

typedef struct {
    int uni_n, uni_m;
    uni_sa_t uni_id;  // array XXX
    ref_off_t uni_off, uni_loc_len; // array XXX
    int read_off, read_loc_len;
} uni_loc_t; // for LOB

typedef struct {
    int lob_flag; // -1:NULL, 0/1: 1 LOB
    int cur_i;
    uni_loc_t lob[2];
} lob_t;

typedef struct {
    uni_sa_t uid;
    ref_off_t uni_off, uni_loc_len;
    int read_off, read_loc_len;
} loc_t; // MEM/LOB

typedef struct {
    int n, m; // number of seeds
    loc_t *loc; int *clu_w;
} loc_clu_t;

typedef struct {
    int n, m;   // number of seeds
    loc_t *loc; // sort by seeds' read_off
} seed_loc_t;

typedef struct {
    int n_thread;  // number of threads

    int mem_len;  //
    int debwt_uni_occ_thd;     
    char pre[1024];
} nano_clu_para;

typedef struct {
    vote_t *v;
} nano_seq_t;

typedef struct {
    int tid;

    debwt_t *db;       // index
    uint8_t *pac;
    bntseq_t *bns;

    int n_seqs;         // read seqs
    kseq_t *w_seqs;
    vote_t *v;          // vote for each cluster

    FILE *clu, *unclu;

    nano_clu_para *cp; // clu parameters
    // aux data during alignment
    // alignment result
} nano_aux_t;


int nano_clu(int argc, char *argv[]);
void realloc_seed_loc(seed_loc_t *s);

#endif
