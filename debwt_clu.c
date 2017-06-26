#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debwt_index.h"
#include "nano_clu.h"
#include "debwt_clu.h"
#include "debwt.h"
#include "kmer_hash.h"
#include "bntseq.h"
#include "utils.h"

#define MEM_LEN 14
#define LOB_LEN 11
#define LOB_DIS 3


seed_loc_t *init_seed_loc()
{
    seed_loc_t *loc = (seed_loc_t*)_err_malloc(sizeof(seed_loc_t));
    loc->n = 0, loc->m = 100;
    loc->loc = (loc_t*)_err_malloc(loc->m * sizeof(loc_t));
    int i; 
    for (i = 0; i < loc->m; ++i) {
        loc->loc[i].uid = 0;
        loc->loc[i].uni_off = 0, loc->loc[i].uni_loc_len = 0;
        loc->loc[i].read_off = 0, loc->loc[i].read_loc_len = 0;
    }
    return loc;
}
void realloc_seed_loc(seed_loc_t *s)
{
    s->m <<= 1;
    s->loc = (loc_t*)_err_realloc(s->loc, s->m * sizeof(loc_t));
}

void free_seed_loc(seed_loc_t *loc) { free(loc->loc); free(loc); }


void uni_pos_print(uni_sa_t uid, debwt_t *db, const bntseq_t *bns)
{
    uint32_t m; uint32_t pos;
    int rid;
    for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m) {
        pos = db->uni_pos[m];
        rid = bns_pos2rid(bns, pos);
        if (rid == -1) err_fatal(__func__, "pos: %d\n", pos+1);
        stdout_printf("%s\t%c\n", bns->anns[rid].name, "+-"[_debwt_get_strand(db->uni_pos_strand, m)]);
    }
}

#define _add_vote(v, l, p, pc1, pc2, n, m) { \
    int _flag, _k_i;    \
    _bin_insert_idx(v, p, n, m, int, _flag, _k_i)   \
    if (_flag == 0) {                \
        if (n == m) {               \
            int _m1 = m, _m2 = m; \
            _realloc(p, _m1, int)    \
            _realloc(pc1, _m2, int)    \
            _realloc(pc2, m, int)    \
        }                           \
        if (_k_i <= n-1) {             \
            memmove(p+_k_i+1, p+_k_i, (n-_k_i)*sizeof(int));  \
            memmove(pc1+_k_i+1, pc1+_k_i, (n-_k_i)*sizeof(int));  \
            memmove(pc2+_k_i+1, pc2+_k_i, (n-_k_i)*sizeof(int));  \
        }   \
        (p)[_k_i] = v;               \
        (pc1)[_k_i] = l; \
        (pc2)[_k_i] = 1; \
        (n)++;                      \
    } else {                        \
        ((pc1)[_k_i])+=l;   \
        ((pc2)[_k_i])+=1;   \
    }       \
}

void nano_add_vote(vote_t *v, int len, uni_sa_t uid, debwt_t *db, const bntseq_t *bns)
{
    uint32_t m, uni_n = db->uni_pos_c[uid+1] - db->uni_pos_c[uid]; uint32_t pos;
    int rid;
    for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m) {
        pos = db->uni_pos[m];
        rid = bns_pos2rid(bns, pos);
        if (rid == -1) err_fatal(__func__, "pos: %d\n", pos+1);
        if (_debwt_get_strand(db->uni_pos_strand, m) == 1) rid = -rid;
        _add_vote(rid, len/uni_n, v->vote_id, v->vote_score, v->hit, v->n, v->m)
#ifdef __DEBUG__
        stdout_printf("vote: %s\t%c\t%d\n", bns->anns[rid].name, "+-"[_debwt_get_strand(db->uni_pos_strand, m)], db->uni_pos[m]);
#endif
    }
}

int bi_extend(uint8_t *seq1, uint8_t *seq2, int off1, int off2, int len, int *l1, int *l2)
{
    int i, len1=0, len2=0;

    for (i = off1-1; i >=0; --i) {
        if (seq1[i] != seq2[i]) break;
        ++len1;
    }
    for (i = off2+1; i < len; ++i) {
        if (seq1[i] != seq2[i]) break;
        ++len2;
    }
    *l1 = len1, *l2 = len2;
    return len1+len2+off2-off1+1;
}

/*
 * uni_mem: directly match with pac XXX
 */
/*
 * read         [----- off1 ----- off2 -----]
 * unipath [ ... ----- off1 ----- off2 ---------- ... ]
 */
int uni_mem(uint8_t *read_seq, int read_len, uint32_t read_off1, uint32_t read_off2, ref_off_t uni_off1, ref_off_t uni_off2, debwt_t *db, uni_sa_t uid, bntseq_t *bns, uint8_t *pac, int *l1, int *l2)
{
    uni_sa_t uni_pos_i = db->uni_pos_c[uid];
    ref_off_t uni_pos = db->uni_pos[uni_pos_i], who_uni_len = db->uni_len[uid];
    ref_off_t uni_start, uni_len, uni_len1, uni_len2, off1, off2; // off1/2: for both read and unipath
    ref_off_t pac_coor;
    uint8_t *_read_seq = read_seq;

    if (uni_off1 >= read_off1) {
        off1 = read_off1;
        uni_len1 = read_off1; uni_start = uni_off1-read_off1;
        who_uni_len -= (uni_off1-read_off1), uni_off2 -= (uni_off1-read_off1);
    } else {
        off1 = uni_off1; 
        uni_len1 = uni_off1; uni_start = 0;
        _read_seq = read_seq+(read_off1-uni_off1);
        read_len -= (read_off1-uni_off1), read_off2 -= (read_off1-uni_off1);
    }
    off2 = read_off2; // read_off2 == uni_off2
    if ((who_uni_len-uni_off2) >= (read_len-read_off2)) {
        uni_len2 = read_len-read_off2-1;
    } else {
        uni_len2 = who_uni_len-uni_off2-1;
    }
    uni_len = uni_len1+off2-off1+1+uni_len2;

    uint8_t is_rev = _debwt_get_strand(db->uni_pos_strand, uni_pos_i);
    if (is_rev) pac_coor = uni_pos-uni_start-uni_len+1;
    else pac_coor = uni_pos+uni_start;
    uint8_t *uni_seq = _bns_get_seq(bns->l_pac, pac, pac_coor, uni_len, is_rev);
    int mem_l = bi_extend(_read_seq, uni_seq, off1, off2, uni_len, l1, l2);
    free(uni_seq);
    return mem_l;
}

int push_loc(seed_loc_t *loc, uni_loc_t uloc)
{
    if (loc->n == loc->m) realloc_seed_loc(loc);

    loc->loc[loc->n].uid = uloc.uni_id;
    loc->loc[loc->n].uni_off = uloc.uni_off, loc->loc[loc->n].uni_loc_len = uloc.uni_loc_len;
    loc->loc[loc->n].read_off = uloc.read_off, loc->loc[loc->n].read_loc_len = uloc.read_loc_len;
    ++loc->n;
    return uloc.read_off;
}

int push_lob(seed_loc_t *loc, lob_t lob) // flag == 0
{
    if (loc->n == loc->m) realloc_seed_loc(loc);

    int lob_i = lob.cur_i;

    loc->loc[loc->n].uid = lob.lob[lob_i].uni_id;
    loc->loc[loc->n].uni_off = lob.lob[lob_i].uni_off, loc->loc[loc->n].read_off = lob.lob[lob_i].read_off;
    // calcu len1,len2 and next_seed_i
    loc->loc[loc->n].uni_loc_len = lob.lob[lob_i].uni_loc_len, loc->loc[loc->n].read_loc_len = lob.lob[lob_i].read_loc_len;
    ++loc->n;
    return lob.lob[lob_i].read_off;
}

int lob_check(uni_loc_t l1, uni_loc_t l2, debwt_t *db)
{
    //int i, j;
    //for (i = 0; i < l1.uni_n; ++i) { // XXX
        //for (j = 0; j < l2.uni_n; ++j) { //XXX
            if (l1.uni_id == l2.uni_id) { // same unipath
                if (abs(l1.read_off-(l2.read_off + l2.read_loc_len)) <= LOB_DIS) { 
                    int uni_off2 = l2.uni_off, uni_off1 = l1.uni_off, uni_loc_len2 = l2.uni_loc_len;
                    if (abs((l1.read_off-(l2.read_off+l2.read_loc_len)) - (uni_off1-(uni_off2+uni_loc_len2))) <= LOB_DIS) return 1; // XXX LOB_DIS ?? LOB_DIFF
                }
            } else { // different unipath
                //uint32_t m, n;
                //for (m = db->uni_pos_c[l1.uni_id[i]]; m < db->uni_pos_c[l1.uni_id[i]+1]; ++m) {
                //    for (n = db->uni_pos_c[l2.uni_id[j]]; n < db->uni_pos_c[l2.uni_id[j]+1]; ++n) {
                //        if (abs((int)(db->uni_pos[m] + l1.off + l1.len) - (int)(db->uni_pos[n] + l2.off)) <= LOB_DIS) return 1;
                //    }
                //}
            }
        //}
    //}
    return 0;
}

// return value:
#define ONE_LOB 0
#define TWO_LOB 1
#define UNMERG_LOB 2
#define MERGED_LOB 3
// [empty] flag = -1
// [lob] flag = 0
// [LOB] flag = 1
int push_1lob(lob_t *lob, uni_loc_t uni_loc, debwt_t *db)
{
    int cur_i = lob->cur_i;
    if (lob->lob_flag == -1) { // NULL
        lob->lob[0] = uni_loc;
        lob->lob_flag = 0;
        lob->cur_i = 0;
        return ONE_LOB;
    } else if (lob->lob_flag == 0 || lob->lob_flag == 1) {
        // check_lob(lob[0/1] and new loc)
        if (lob_check(lob->lob[cur_i], uni_loc, db)) {
            // update LOB
            lob->lob[cur_i].uni_off = uni_loc.uni_off;
            lob->lob[cur_i].read_off = uni_loc.read_off;
            lob->lob[cur_i].uni_loc_len = lob->lob[cur_i].uni_off-uni_loc.uni_off+lob->lob[cur_i].uni_loc_len;
            lob->lob[cur_i].read_loc_len = lob->lob[cur_i].read_off-uni_loc.read_off+lob->lob[cur_i].read_loc_len;
            lob->lob_flag = 1;
            return MERGED_LOB;
        } else {
            lob->lob[1-cur_i] = uni_loc;
            lob->cur_i = 1 - cur_i;
            if (lob->lob_flag == 0) return TWO_LOB;
            else { lob->lob_flag = 0; return UNMERG_LOB; }
        }
    } else {
        err_printf("[push_lob] Error: unknown lob flag: %d.\n", lob->lob_flag);
        exit(-1);
    }
}

void set_uni_loc(uni_loc_t *uni_loc, int read_off, int read_loc_len, uni_sa_t uid, ref_off_t uni_off, ref_off_t uni_loc_len)
{
    uni_loc->read_off = read_off;
    uni_loc->read_loc_len = read_loc_len;
    uni_loc->uni_id = uid;
    uni_loc->uni_off = uni_off;
    uni_loc->uni_loc_len = uni_loc_len;
    uni_loc->uni_n = 1;
    uni_loc->uni_m = 1;
}

/*int debwt_gen_loc_clu(uint8_t *bseq, int seq_len, debwt_t *db, bntseq_t *bns, uint8_t *pac, nano_clu_para *cp, vote_t *v)
{
    int cur_i, old_i, old_lob_i;
    debwt_count_t i, uni_occ_thd = cp->debwt_uni_occ_thd, k, l, il;
    seed_loc_t *loc_clu = init_seed_loc();
    lob_t *lob = (lob_t*)_err_malloc(sizeof(lob_t)); lob->lob_flag = -1;
    uni_loc_t uni_loc;

    for (cur_i = seq_len - _BWT_HASH_K; cur_i >= 0; --cur_i) {
        printf("cur_i: %d\n", cur_i);
        old_i = cur_i;
        // debwt hash
        k = db->bwt_hash[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        il = db->bwt_hash_il[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        if (il == 0) continue;
        l = k + il - 1;
        // bwt backtrack
        while (il > uni_occ_thd && cur_i >= 1) {
            il = debwt_exact_match_alt(db, 1, bseq+cur_i-1, &k, &l);
            if (il == 0 || cur_i == 1) break;
            --cur_i;
        }
        if (il == 0) continue;

        // base extend (back/forward) // XXX NOT span unipaths
        int max_len = 0, m_len; ref_off_t uni_off, max_uni_off=0, max_loc_len1=0; uni_sa_t uid, max_uid=-1;
        int l1, l2, max_read_off=0, max_loc_len2=0;
        for (i = 0; i < il; ++i) { // 0 < il <= M
            uid = debwt_sa(db, k+i, &uni_off);
            // XXX
            m_len = uni_mem(bseq, seq_len, cur_i, old_i+_BWT_HASH_K-1, uni_off, uni_off+(old_i-cur_i)+_BWT_HASH_K-1, db, uid, bns, pac, &l1, &l2);
            // for one unipath, only one MEM location
            // for one seed, keep one longest MEM-unipath and (THD-1) secondary MEM-unipath(>l/2)
            // need a heap XXX
            if (m_len > max_len) {
                max_uid = uid;
                max_uni_off = uni_off-l1, max_loc_len1 = m_len;
                max_read_off = cur_i-l1, max_loc_len2 = m_len;
                max_len = m_len;
            }
        }
        // next loop
        cur_i = old_i;
        if (max_len > 0) set_uni_loc(&uni_loc, max_read_off, max_loc_len2, max_uid, max_uni_off, max_loc_len1);
        if (max_len >= MEM_LEN) { // MEM seed
            if (lob->lob_flag == 1) {
                lob->cur_i = 1 - lob->cur_i;
                push_lob(loc_clu, *lob);
                loc_t l = loc_clu->loc[loc_clu->n-1];
                stdout_printf("LOB id: %d, uni_off: %d, read_off: %d, len1: %d, len2: %d\n", l.uid, l.uni_off, l.read_off, l.uni_loc_len, l.read_loc_len);
        nano_add_vote(v, l.uni_loc_len, l.uid, db, bns);
            }
            lob->lob_flag = -1;
            cur_i = push_loc(loc_clu, uni_loc) - _BWT_HASH_K; // push mem loc
            cur_i = uni_loc.read_off - _BWT_HASH_K;
            stdout_printf("MEM: id: %d, uni_off: %d, read_off: %d, len: %d\n", uni_loc.uni_id, uni_loc.uni_off, uni_loc.read_off, uni_loc.uni_loc_len);
            nano_add_vote(v, uni_loc.uni_loc_len, uni_loc.uni_id, db, bns);
        } else if (max_len >= LOB_LEN) {
            // XXX bug exists
            int res = push_1lob(lob, uni_loc, db);
            if (res == ONE_LOB) { // lob_flag == -1
                cur_i = uni_loc.read_off - _BWT_HASH_K;
                old_lob_i = old_i;
            } else if (res == TWO_LOB) { // lob_flag == 0
                cur_i = old_lob_i;
            } else if (res == UNMERG_LOB) { // push lob->lob[1-cur_i]
                cur_i = push_lob(loc_clu, *lob) - _BWT_HASH_K;
                old_lob_i = cur_i;
                loc_t l = loc_clu->loc[loc_clu->n-1];
                stdout_printf("LOB id: %d, uni_off: %d, read_off: %d, len1: %d, len2: %d\n", l.uid, l.uni_off, l.read_off, l.uni_loc_len, l.read_loc_len);
                nano_add_vote(v, l.uni_loc_len, l.uid, db, bns);
            } else { // MERGED_LOB
                cur_i = uni_loc.read_off - _BWT_HASH_K;
                old_lob_i = cur_i;
            }
        }
    }
    if (lob->lob_flag == 1) {
        lob->cur_i = 1 - lob->cur_i;
        push_lob(loc_clu, *lob);
        loc_t l = loc_clu->loc[loc_clu->n-1];
        stdout_printf("LOB id: %d, uni_off: %d, read_off: %d, len1: %d, len2: %d\n", l.uid, l.uni_off, l.read_off, l.uni_loc_len, l.read_loc_len);
        nano_add_vote(v, l.uni_loc_len, l.uid, db, bns);
    }
    
    free(lob); free_seed_loc(loc_clu);
    return 0;
}*/

int debwt_gen_loc_clu(uint8_t *bseq, int seq_len, debwt_t *db, bntseq_t *bns, uint8_t *pac, nano_clu_para *cp, vote_t *v)
{
    int cur_i, old_i;
    debwt_count_t i, uni_occ_thd = cp->debwt_uni_occ_thd, k, l, il;
    uni_loc_t uni_loc;

    for (cur_i = seq_len - _BWT_HASH_K; cur_i >= 0; --cur_i) {
        //printf("cur_i: %d\n", cur_i);
        old_i = cur_i;
        // debwt hash
        k = db->bwt_hash[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        il = db->bwt_hash_il[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        if (il == 0) continue;
        l = k + il - 1;
        // bwt backtrack
        while (il > uni_occ_thd && cur_i >= 1) {
            il = debwt_exact_match_alt(db, 1, bseq+cur_i-1, &k, &l);
            if (il == 0 || cur_i == 1) break;
            --cur_i;
        }
        if (il == 0) continue;

        // base extend (back/forward) // XXX NOT span unipaths
        int max_len = 0, m_len; ref_off_t uni_off, max_uni_off=0, max_loc_len1=0; uni_sa_t uid, max_uid=-1;
        int l1, l2, max_read_off=0, max_loc_len2=0;
        for (i = 0; i < il; ++i) { // 0 < il <= M
            uid = debwt_sa(db, k+i, &uni_off);
            // XXX
            m_len = uni_mem(bseq, seq_len, cur_i, old_i+_BWT_HASH_K-1, uni_off, uni_off+(old_i-cur_i)+_BWT_HASH_K-1, db, uid, bns, pac, &l1, &l2);
            // for one unipath, only one MEM location
            // for one seed, keep one longest MEM-unipath and (THD-1) secondary MEM-unipath(>l/2)
            // need a heap XXX
            if (m_len > max_len) {
                max_uid = uid;
                max_uni_off = uni_off-l1, max_loc_len1 = m_len;
                max_read_off = cur_i-l1, max_loc_len2 = m_len;
                max_len = m_len;
            }
        }
        // next loop
        cur_i = old_i;
        if (max_len > 0) set_uni_loc(&uni_loc, max_read_off, max_loc_len2, max_uid, max_uni_off, max_loc_len1);
        if (max_len >= cp->mem_len) { // MEM seed
            cur_i = uni_loc.read_off - _BWT_HASH_K; // push mem loc
#ifdef __DEBUG__
            stdout_printf("MEM: id: %d, uni_off: %d, read_off: %d, len: %d\n", uni_loc.uni_id, uni_loc.uni_off, uni_loc.read_off, uni_loc.uni_loc_len);
#endif
            nano_add_vote(v, uni_loc.uni_loc_len, uni_loc.uni_id, db, bns);
        }
    }
    
    return 0;
}
