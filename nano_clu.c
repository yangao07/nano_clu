#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>
#include "nano_clu.h"
#include "debwt_clu.h"
#include "debwt.h"
#include "utils.h"
#include "kseq.h"

extern const char PROG[20];

int nano_clu_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s clu [option] <gene_ref.fa> <read.fa/fq>\n\n", PROG);
    err_printf("Input Options:\n\n");
    err_printf("         -t --thread      [INT]    number of thread. [%d]\n", 1);
    err_printf("         -l --mem-len     [INT]    minimum length of mem seed to count. [%d]\n", NANO_MEM_LEN);
    err_printf("         -o --uni-occ     [INT]    maximum occurrence of seed's hit to end the debwt backtracking. [%d]\n", NANO_UNI_OCC_THD);
    err_printf("         -O --out-prefix  [STR]    prefix of output files.\n");
    err_printf("\n");
	return 1;
}

nano_clu_para *nano_init_cp(void)
{
    // XXX init para
    nano_clu_para *cp = (nano_clu_para*)calloc(1, sizeof(nano_clu_para));
    cp->n_thread = 1;
    cp->mem_len = NANO_MEM_LEN;
    cp->debwt_uni_occ_thd = NANO_UNI_OCC_THD;
    return cp;
}

vote_t *init_vote(int v_n)
{
    vote_t *v = (vote_t*)_err_malloc(v_n * sizeof(vote_t));
    int i; for (i = 0; i < v_n; ++i) {
        v[i].n = 0, v[i].m = 10;
        v[i].vote_id = (int*)_err_malloc(10 * sizeof(int));
        v[i].vote_score = (int*)_err_malloc(10 * sizeof(int));
    }
    return v;
}

void free_vote(vote_t *v, int v_n)
{
    int i;
    for (i = 0; i < v_n; ++i) {
        free(v[i].vote_id); free(v[i].vote_score);
    }
    free(v);
}

void aux_free(nano_aux_t *aux)
{
    int i;
    debwt_index_free(aux->db); free(aux->db);
    free(aux->pac); bns_destroy(aux->bns);

    for (i = 0; i < CHUNK_READ_N; ++i) {
        if (aux->w_seqs+i != NULL) {
            free((aux->w_seqs+i)->name.s);
            free((aux->w_seqs+i)->comment.s);
            free((aux->w_seqs+i)->seq.s);
            free((aux->w_seqs+i)->qual.s);
        }
    }
    ks_destroy(aux->w_seqs->f);
    free(aux->w_seqs);
    free_vote(aux->v, CHUNK_READ_N); 
    fclose(aux->clu); fclose(aux->unclu);
    free(aux->cp);
    free(aux);
}

int nano_read_seq(kseq_t *read_seq, int chunk_read_n)
{
    kseq_t *s = read_seq; int n = 0;
    while (kseq_read(s+n) >= 0) {
        if (chunk_read_n-1 == n++) break;
    }
    return n;
}

int nano_output_clu(nano_aux_t *aux)
{
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; vote_t *v = aux->v;
    int i, j, max, max_id, sec, sec_id;
    for (i = 0; i < n_seqs; ++i) {
        kseq_t *seqs = w_seqs + i;
        if (v[i].n > 0) {
            max = 0; max_id = v[i].vote_id[0];
            sec = 0; sec_id = -1;
            for (j = 0; j < v[i].n; ++j) {
                if (v[i].vote_score[j] > max) {
                    max = v[i].vote_score[j];
                    max_id = v[i].vote_id[j];
                } else if (v[i].vote_score[j] > sec) {
                    sec = v[i].vote_score[j];
                    sec_id = v[i].vote_id[j];
                }
            }
            // output to .clu.fa
#ifdef __DEBUG__
            fprintf(stdout, "%s\nmax: %d\tmax_socre: %d\tsec: %d\tsec_score: %d\n", seqs->name.s, max_id, max, sec_id, sec);
#endif
            fprintf(aux->clu, ">%s %s\n", seqs->name.s, aux->bns->anns[max_id].name);
            fprintf(aux->clu, "%s\n", seqs->seq.s);
        } else { // output to .unclu.fa
            fprintf(aux->unclu, ">%s\n", seqs->name.s);
            fprintf(aux->unclu, "%s\n", seqs->seq.s);
        }
    }
    return 0;
}

int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

int nano_cal_clu(debwt_t *db, bntseq_t *bns, uint8_t *pac, kseq_t *seqs, nano_clu_para *cp, vote_t *v)
{
    uint32_t i;
    uint8_t *bseq = (uint8_t*)_err_malloc(seqs->seq.l * sizeof(uint8_t));
    for (i = 0; i < seqs->seq.l; ++i) bseq[i] = nst_nt4_table[(int)(seqs->seq.s[i])];


    // seeding and locating
    v->n = 0;
#ifdef __DEBUG__
    stdout_printf("%s\n%s\n", seqs->name.s, seqs->seq.s);
#endif
    debwt_gen_loc_clu(bseq, seqs->seq.l, db, bns, pac, cp, v);
    // output read to gene cluster file

    /* exact match test
    debwt_count_t ok, ol, l;
    uint32_t uid, off, m;
    l = debwt_exact_match(db, seqs->seq.l, bseq, &ok, &ol);
    for (i = 0; i < l; ++i) {
        uid = debwt_sa(db, ok+i, &off);
        stdout_printf("UID: #%d\n", uid);
        for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m)
            stdout_printf("%c%d\t%d\n", "+-"[_debwt_get_strand(db->uni_pos_strand, m)], (int)db->uni_pos[m]+1, off);
    }
    */
    free(bseq);
    return 0;
}

int nano_main_clu(nano_aux_t *aux)
{
    debwt_t *db = aux->db; bntseq_t *bns = aux->bns; uint8_t *pac = aux->pac;
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; 
    nano_clu_para *cp = aux->cp; vote_t *w_v = aux->v;
    int i_seq=0;
    while (i_seq < n_seqs) {
        if (i_seq == n_seqs) break;
        kseq_t *seqs = w_seqs+i_seq;
        vote_t *v = w_v+i_seq;
        nano_cal_clu(db, bns, pac, seqs, cp, v);
        i_seq++;
    }
    return 0;
}

static void *nano_thread_main_clu(void *a)
{
    nano_aux_t *aux = (nano_aux_t*)a;
    debwt_t *db = aux->db; bntseq_t *bns = aux->bns; uint8_t *pac = aux->pac;
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; 
    nano_clu_para *cp = aux->cp; vote_t *w_v = aux->v;
    int i_seq;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i_seq = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i_seq >= n_seqs) break;
        kseq_t *seqs = w_seqs+i_seq;
        vote_t *v = w_v+i_seq;
        nano_cal_clu(db, bns, pac, seqs, cp, v);
    }
    return 0;
}

int nano_clu_core(const char *ref_fn, const char *read_fn, nano_clu_para *nano_cp)
{
    /* load index */
    err_printf("[nano_clu_core] Restoring ref-indices ... ");
    debwt_t *db_idx = debwt_restore_index(ref_fn);
    bntseq_t *bns = bns_restore(ref_fn);
    uint8_t *pac = (uint8_t*)_err_calloc(bns->l_pac/4+1, 1);
    fread(pac, 1, bns->l_pac/4+1, bns->fp_pac); err_printf("done!\n");

    // print debwt
    // int i, rid, pos;
    // for (i = 0; i < db_idx->n_unipath; ++i) {
    //     uint32_t m;
    //     for (m = db_idx->uni_pos_c[i]; m < db_idx->uni_pos_c[i+1]; ++m) {
    //         pos = db_idx->uni_pos[m];
    //         rid = bns_pos2rid(bns, pos);
    //         pos -= (bns->anns[rid].offset + rid);
    //         printf("uni: %d len: %d\t", pos, db_idx->uni_len[i]);
    //         printf("%s\t%c\n", bns->anns[rid].name, "+-"[_debwt_get_strand(db_idx->uni_pos_strand, m)]);
    //     }
    // }
    
    /* output cluster file */

    /* read read.fa/fq */
    int n_seqs, i;
    gzFile readfp; kseq_t *read_seqs;
    readfp = err_xzopen_core(__func__, read_fn, "r");
    kstream_t *fs = ks_init(readfp);
    read_seqs = (kseq_t*)_err_calloc(CHUNK_READ_N, sizeof(kseq_t));
    for (i = 0; i < CHUNK_READ_N; ++i) read_seqs[i].f = fs;

    // alloc and init for auxiliary data
    if (nano_cp->n_thread < 1) nano_cp->n_thread = 1;
    nano_aux_t *aux = (nano_aux_t*)_err_calloc(nano_cp->n_thread, sizeof(nano_aux_t));
    for (i = 0; i < nano_cp->n_thread; ++i) {
        aux[i].tid = i; aux[i].cp = nano_cp;
        aux[i].db = db_idx; aux[i].bns = bns; aux[i].pac = pac;
    }
    vote_t *v = init_vote(CHUNK_READ_N);
    char clu_name[1024], unclu_name[1024];
    strcpy(clu_name, nano_cp->pre); strcat(clu_name, ".clu.fa");
    strcpy(unclu_name, nano_cp->pre); strcat(unclu_name, ".unclu.fa");
    FILE *clu = fopen(clu_name, "w"), *unclu = fopen(unclu_name, "w");
 
    if (nano_cp->n_thread <= 1) {
        while ((n_seqs = nano_read_seq(read_seqs, CHUNK_READ_N)) != 0) { 
            aux->n_seqs = n_seqs;
            aux->w_seqs = read_seqs;
            aux->v = v;
            aux->clu = clu;
            aux->unclu = unclu;
            nano_main_clu(aux);
            nano_output_clu(aux);
        }
    } else { // multi-threads
        pthread_rwlock_init(&RWLOCK, NULL);
        while ((n_seqs = nano_read_seq(read_seqs, CHUNK_READ_N)) != 0) { 
            THREAD_READ_I = 0;
            pthread_t *tid; pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            tid = (pthread_t*)_err_calloc(nano_cp->n_thread, sizeof(pthread_t));
            int j;
            for (j = 0; j < nano_cp->n_thread; ++j) {
                aux[j].n_seqs = n_seqs;
                aux[j].w_seqs = read_seqs;
                aux[j].v = v;
                aux[j].clu = clu;
                aux[j].unclu = unclu;
                pthread_create(&tid[j], &attr, nano_thread_main_clu, aux+j);
            }
            for (j = 0; j < nano_cp->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
            nano_output_clu(aux);
        }
        pthread_rwlock_destroy(&RWLOCK);
    }
    aux_free(aux);
    err_gzclose(readfp);
    return 0;
}

int nano_clu(int argc, char *argv[])
{
    int c;
    nano_clu_para *nano_cp = nano_init_cp();

    while ((c = getopt(argc, argv, "t:l:o:O:")) >= 0) {
        switch (c)
        {
            case 't': nano_cp->n_thread = atoi(optarg); break;
            case 'l': nano_cp->mem_len = atoi(optarg); break;
            case 'o': nano_cp->debwt_uni_occ_thd = atoi(optarg); break;
            case 'O': strcpy(nano_cp->pre, optarg); break;
            default: return nano_clu_usage();
        }
    }
    if (argc - optind != 2) return nano_clu_usage();

    if (strlen(nano_cp->pre) == 0) strcpy(nano_cp->pre, argv[optind+1]);

    nano_clu_core(argv[optind], argv[optind+1], nano_cp);
    return 0;
}
