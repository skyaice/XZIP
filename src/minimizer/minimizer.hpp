/*
 * minimizer.hpp
 *
 *  Created on: 2022年4月2日
 *      Author: fenghe
 */

#ifndef MINIMIZER_HPP_
#define MINIMIZER_HPP_

#include <vector>
#include <cstring>
extern "C"{
	#include "../clib/utils.h"
	#include "../clib/ksort.h"
}

struct mm128_t
{
	uint64_t x, y;
	void print(FILE* o){
		fprintf(stderr, "hash value: %ld, kmerSpan %ld, rid %ld, lastPos %ld strand %ld\n", (x >> 8), (x & 0xff), (y >> 32), ((y & 0xffffffff) >> 1), (y & 1) );
	}
} ;

class Minimizer_generater{

public:
	void init(int window_size_, int mini_kmer_len_, int is_hpc_){
		window_size = window_size_;
		mini_kmer_len = mini_kmer_len_;
		is_hpc = is_hpc_;
	}

	void print_minimizer_list(std::vector<mm128_t> &p){
		for(auto & m : p)
			m.print(stderr);
	}

private:

	int window_size;
	int mini_kmer_len;
	int is_hpc;

	static inline uint64_t hash64(uint64_t key, uint64_t mask)
	{
		key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
		key = key ^ key >> 24;
		key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
		key = key ^ key >> 14;
		key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
		key = key ^ key >> 28;
		key = (key + (key << 31)) & mask;
		return key;
	}

	typedef struct { // a simplified version of kdq
		int front, count;
		int a[32];
	} tiny_queue_t;

	static inline void tq_push(tiny_queue_t *q, int x)
	{
		q->a[((q->count++) + q->front) & 0x1f] = x;
	}

	static inline int tq_shift(tiny_queue_t *q)
	{
		int x;
		if (q->count == 0) return -1;
		x = q->a[q->front++];
		q->front &= 0x1f;
		--q->count;
		return x;
	}
public:
	/**
	 * COPY from Minimape (Heng Li)
	 * Find symmetric (w,k)-minimizers on a DNA sequence
	 *
	 * @param seq    DNA sequence
	 * @param w      find a minimizer for every $w consecutive k-mers
	 * @param k      k-mer size
	 * @param rid    reference ID; will be copied to the output $p array
	 * @param is_hpc homopolymer-compressed or not
	 * @param p      minimizers
	 *               p->a[i].x = kMer<<8 | kmerSpan
	 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
	 *               where lastPos is the position of the last base of the i-th minimizer,
	 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
	 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
	 */
	void mm_sketch(const char *seq, int seq_len, std::vector<mm128_t> &p);
	void generate(int argc, char *argv[]);

};
#endif /* MINIMIZER_HPP_ */

//
///*
// * minimizer_idx.hpp
// *
// */
//
//#ifndef MINIMIZER_IDX_HPP_
//#define MINIMIZER_IDX_HPP_
//
//#include "minimizer.hpp"
//#include <stdint.h>
//#include <stdlib.h>
//#include <string>
//extern "C" {
//#include "../clib/kthread.h"
//#include "../clib/ksort.h"
//}
//
//#define MM_F_NO_DIAG       0x001 // no exact diagonal hit
//#define MM_F_NO_DUAL       0x002 // skip pairs where query name is lexicographically larger than target name
//#define MM_F_CIGAR         0x004
//#define MM_F_OUT_SAM       0x008
//#define MM_F_NO_QUAL       0x010
//#define MM_F_OUT_CG        0x020
//#define MM_F_OUT_CS        0x040
//#define MM_F_SPLICE        0x080 // splice mode
//#define MM_F_SPLICE_FOR    0x100 // match GT-AG
//#define MM_F_SPLICE_REV    0x200 // match CT-AC, the reverse complement of GT-AG
//#define MM_F_NO_LJOIN      0x400
//#define MM_F_OUT_CS_LONG   0x800
//#define MM_F_SR            0x1000
//#define MM_F_FRAG_MODE     0x2000
//#define MM_F_NO_PRINT_2ND  0x4000
//#define MM_F_2_IO_THREADS  0x8000
//#define MM_F_LONG_CIGAR    0x10000
//#define MM_F_INDEPEND_SEG  0x20000
//#define MM_F_SPLICE_FLANK  0x40000
//#define MM_F_SOFTCLIP      0x80000
//#define MM_F_FOR_ONLY      0x100000
//#define MM_F_REV_ONLY      0x200000
//#define MM_F_HEAP_SORT     0x400000
//#define MM_F_ALL_CHAINS    0x800000
//#define MM_F_OUT_MD        0x1000000
//#define MM_F_COPY_COMMENT  0x2000000
//#define MM_F_EQX           0x4000000 // use =/X instead of M
//#define MM_F_PAF_NO_HIT    0x8000000 // output unmapped reads to PAF
//#define MM_F_NO_END_FLT    0x10000000
//#define MM_F_HARD_MLEVEL   0x20000000
//#define MM_F_SAM_HIT_ONLY  0x40000000
//#define MM_F_RMQ           (0x80000000LL)
//#define MM_F_QSTRAND       (0x100000000LL)
//#define MM_F_NO_INV        (0x200000000LL)
//#define MM_F_NO_HASH_NAME  (0x400000000LL)
//
//#define MM_I_HPC          0x1
//#define MM_I_NO_SEQ       0x2
//#define MM_I_NO_NAME      0x4
//
//#define MM_IDX_MAGIC   "MMI\2"
//
//#define MM_MAX_SEG       255
//
//#define MM_CIGAR_MATCH      0
//#define MM_CIGAR_INS        1
//#define MM_CIGAR_DEL        2
//#define MM_CIGAR_N_SKIP     3
//#define MM_CIGAR_SOFTCLIP   4
//#define MM_CIGAR_HARDCLIP   5
//#define MM_CIGAR_PADDING    6
//#define MM_CIGAR_EQ_MATCH   7
//#define MM_CIGAR_X_MISMATCH 8
//
//#define MM_CIGAR_STR  "MIDNSHP=XB"
//
//// emulate 128-bit integers and arrays
//typedef struct mm_idx_bucket_s {
//	std::vector<mm128_t> a;
//	int32_t n;   // size of the _p_ array
//	uint64_t *p; // position array for minimizers appearing >1 times
//
//	std::unordered_map<uint64_t, uint64_t> h; // hash table indexing _p_ and minimizers appearing once ???
//} mm_idx_bucket_t;
//
//typedef struct {
//	int32_t st, en, max; // max is not used for now
//	int32_t score :30, strand :2;
//} mm_idx_intv1_t;
//
//typedef struct mm_idx_intv_s {
//	int32_t n, m;
//	mm_idx_intv1_t *a;
//} mm_idx_intv_t;
//
//typedef struct {
//	uint32_t rid;
//	std::string name;
//	std::string seq;
//} mm_bseq1_t;
//
//// minimap2 index
//struct mm_idx_seq_t {
//
//	mm_idx_seq_t(int flag, mm_bseq1_t &s, uint64_t sum_len) {
//		name.clear();
//		if (!(flag & MM_I_NO_NAME))
//			name = s.name;
//		len = s.seq.size();
//		offset = sum_len;
//		is_alt = 0;
//	}
//
//	std::string name; // name of the db sequence
//	uint64_t offset; // offset in mm_idx_t::S
//	uint32_t len;    // length
//	uint32_t is_alt;
//};
//
//typedef struct {
//	int32_t b, w, k, flag;
//
//	int32_t index;
//	int32_t n_alt;
//
//	std::vector<mm_idx_seq_t> seq;
//	std::vector<uint32_t> S;   // 4-bit packed sequence
//	struct mm_idx_bucket_s *B; // index (hidden)
//	struct mm_idx_intv_s *I;   // intervals (hidden)
//	//khash_t(str) *h;
//	//hash table
//} mm_idx_t;
//
//// minimap2 alignment
//typedef struct {
//	uint32_t capacity;                  // the capacity of cigar[]
//	int32_t dp_score, dp_max, dp_max2; // DP score; score of the max-scoring segment; score of the best alternate mappings
//	uint32_t n_ambi :30, trans_strand :2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
//	uint32_t n_cigar;                   // number of cigar operations in cigar[]
//	uint32_t cigar[];
//} mm_extra_t;
//
//typedef struct {
//	int32_t id;             // ID for internal uses (see also parent below)
//	int32_t cnt;            // number of minimizers; if on the reverse strand
//	int32_t rid; // reference index; if this is an alignment from inversion rescue
//	int32_t score;          // DP alignment score
//	int32_t qs, qe, rs, re; // query start and end; reference start and end
//	int32_t parent, subsc; // parent==id if primary; best alternate mapping score
//	int32_t as;             // offset in the a[] array (for internal uses only)
//	int32_t mlen, blen; // seeded exact match length; seeded alignment block length
//	int32_t n_sub;          // number of suboptimal mappings
//	int32_t score0;    // initial chaining score (before chain merging/spliting)
//	uint32_t mapq :8, split :2, rev :1, inv :1, sam_pri :1, proper_frag :1,
//			pe_thru :1, seg_split :1, seg_id :8, split_inv :1, is_alt :1,
//			strand_retained :1, dummy :5;
//	uint32_t hash;
//	float div;
//	mm_extra_t *p;
//} mm_reg1_t;
//
//// indexing and mapping options
//typedef struct {
//	short k, w, flag, bucket_bits;
//	int64_t mini_batch_size;
//	uint64_t batch_size;
//} mm_idxopt_t;
//
//typedef struct {
//	int64_t flag;    // see MM_F_* macros
//	int seed;
//	int sdust_thres; // score threshold for SDUST; 0 to disable
//
//	int max_qlen;    // max query length
//
//	int bw, bw_long; // bandwidth
//	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
//	int max_frag_len;
//	int max_chain_skip, max_chain_iter;
//	int min_cnt;         // min number of minimizers on each chain
//	int min_chain_score; // min chaining score
//	float chain_gap_scale;
//	float chain_skip_scale;
//	int rmq_size_cap, rmq_inner_dist;
//	int rmq_rescue_size;
//	float rmq_rescue_ratio;
//
//	float mask_level;
//	int mask_len;
//	float pri_ratio;
//	int best_n;      // top best_n chains are subjected to DP alignment
//
//	float alt_drop;
//
//	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
//	int sc_ambi; // score when one or both bases are "N"
//	int noncan;      // cost of non-canonical splicing sites
//	int junc_bonus;
//	int zdrop, zdrop_inv; // break alignment if alignment score drops too fast along the diagonal
//	int end_bonus;
//	int min_dp_max; // drop an alignment if the score of the max scoring segment is below this threshold
//	int min_ksw_len;
//	int anchor_ext_len, anchor_ext_shift;
//	float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio
//
//	int rank_min_len;
//	float rank_frac;
//
//	int pe_ori, pe_bonus;
//
//	float mid_occ_frac;  // only used by mm_mapopt_update(); see below
//	float q_occ_frac;
//	int32_t min_mid_occ, max_mid_occ;
//	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
//	int32_t max_occ, max_max_occ, occ_dist;
//	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
//	int64_t max_sw_mat;
//	int64_t cap_kalloc;
//
//	const char *split_prefix;
//} mm_mapopt_t;
//
//// index reader
//typedef struct {
//	int is_idx, n_parts;
//	int64_t idx_size;
//	mm_idxopt_t opt;
//	FILE *fp_out;
//	union {
//		struct mm_bseq_file_s *seq;
//		FILE *idx;
//	} fp;
//} mm_idx_reader_t;
//
///**
// * Initialize an index reader
// *
// * @param fn         index or fasta/fastq file name (this function tests the file type)
// * @param opt        indexing parameters
// * @param fn_out     if not NULL, write built index to this file
// *
// * @return an index reader on success; NULL if fail to open _fn_
// */
//mm_idx_reader_t* mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt,
//		const char *fn_out);
//
///**
// * Read/build an index
// *
// * If the input file is an index file, this function reads one part of the
// * index and returns. If the input file is a sequence file (fasta or fastq),
// * this function constructs the index for about mm_idxopt_t::batch_size bases.
// * Importantly, for a huge collection of sequences, this function may only
// * return an index for part of sequences. It needs to be repeatedly called
// * to traverse the entire index/sequence file.
// *
// * @param r          index reader
// * @param n_threads  number of threads for constructing index
// *
// * @return an index on success; NULL if reaching the end of the input file
// */
//mm_idx_t* mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);
//
///**
// * Destroy/deallocate an index reader
// *
// * @param r          index reader
// */
//void mm_idx_reader_close(mm_idx_reader_t *r);
//
//int mm_idx_reader_eof(const mm_idx_reader_t *r);
//
///**
// * Check whether the file contains a minimap2 index
// *
// * @param fn         file name
// *
// * @return the file size if fn is an index file; 0 if fn is not.
// */
//int64_t mm_idx_is_idx(const char *fn);
//
///**
// * Load a part of an index
// *
// * Given a uni-part index, this function loads the entire index into memory.
// * Given a multi-part index, it loads one part only and places the file pointer
// *  * at the end of that part.
// *
// * @param fp         pointer to FILE object
// *
// * @return minimap2 index read from fp
// */
////
////mm_idx_t *mm_idx_load(FILE *fp)
////{
////	char magic[4];
////	uint32_t x[5], i;
////	uint64_t sum_len = 0;
////	mm_idx_t *mi;
////
////	if (fread(magic, 1, 4, fp) != 4) return 0;
////	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
////	if (fread(x, 4, 5, fp) != 5) return 0;
////	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
////	mi->n_seq = x[3];
////	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
////	for (i = 0; i < mi->n_seq; ++i) {
////		uint8_t l;
////		mm_idx_seq_t *s = &mi->seq[i];
////		fread(&l, 1, 1, fp);
////		if (l) {
////			s->name = (char*)kmalloc(mi->km, l + 1);
////			fread(s->name, 1, l, fp);
////			s->name[l] = 0;
////		}
////		fread(&s->len, 4, 1, fp);
////		s->offset = sum_len;
////		s->is_alt = 0;
////		sum_len += s->len;
////	}
////	for (i = 0; i < 1<<mi->b; ++i) {
////		mm_idx_bucket_t *b = &mi->B[i];
////		uint32_t j, size;
////		khint_t k;
////		idxhash_t *h;
////		fread(&b->n, 4, 1, fp);
////		b->p = (uint64_t*)malloc(b->n * 8);
////		fread(b->p, 8, b->n, fp);
////		fread(&size, 4, 1, fp);
////		if (size == 0) continue;
////		b->h = h = kh_init(idx);
////		kh_resize(idx, h, size);
////		for (j = 0; j < size; ++j) {
////			uint64_t x[2];
////			int absent;
////			fread(x, 8, 2, fp);
////			k = kh_put(idx, h, x[0], &absent);
////			assert(absent);
////			kh_val(h, k) = x[1];
////		}
////	}
////	if (!(mi->flag & MM_I_NO_SEQ)) {
////		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
////		fread(mi->S, 4, (sum_len + 7) / 8, fp);
////	}
////	return mi;
////}
//
///*************
// * index I/O *
// *************/
///**
// * Append an index (or one part of a full index) to file
// *
// * @param fp         pointer to FILE object
// * @param mi         minimap2 index
// */
////
////void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
////{
////	uint64_t sum_len = 0;
////	uint32_t x[5], i;
////
////	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
////	fwrite(MM_IDX_MAGIC, 1, 4, fp);
////	fwrite(x, 4, 5, fp);
////	for (i = 0; i < mi->n_seq; ++i) {
////		if (mi->seq[i].name) {
////			uint8_t l = strlen(mi->seq[i].name);
////			fwrite(&l, 1, 1, fp);
////			fwrite(mi->seq[i].name, 1, l, fp);
////		} else {
////			uint8_t l = 0;
////			fwrite(&l, 1, 1, fp);
////		}
////		fwrite(&mi->seq[i].len, 4, 1, fp);
////		sum_len += mi->seq[i].len;
////	}
////	for (i = 0; i < 1<<mi->b; ++i) {
////		mm_idx_bucket_t *b = &mi->B[i];
////		khint_t k;
////		idxhash_t *h = (idxhash_t*)b->h;
////		uint32_t size = h? h->size : 0;
////		fwrite(&b->n, 4, 1, fp);
////		fwrite(b->p, 8, b->n, fp);
////		fwrite(&size, 4, 1, fp);
////		if (size == 0) continue;
////		for (k = 0; k < kh_end(h); ++k) {
////			uint64_t x[2];
////			if (!kh_exist(h, k)) continue;
////			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
////			fwrite(x, 8, 2, fp);
////		}
////	}
////	if (!(mi->flag & MM_I_NO_SEQ))
////		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
////	fflush(fp);
////}
//
///**
// * Create an index from strings in memory
// *
// * @param w            minimizer window size
// * @param k            minimizer k-mer size
// * @param is_hpc       use HPC k-mer if true
// * @param bucket_bits  number of bits for the first level of the hash table
// * @param n            number of sequences
// * @param seq          sequences in A/C/G/T
// * @param name         sequence names; could be NULL
// *
// * @return minimap2 index
// */
//mm_idx_t* mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n,
//		const char **seq, const char **name);
//
///**
// * Print index statistics to stderr
// *
// * @param mi         minimap2 index
// */
//void mm_idx_stat(const mm_idx_t *idx);
//
///**
// * Destroy/deallocate an index
// *
// * @param r          minimap2 index
// */
//void mm_idx_destroy(mm_idx_t *mi) {
//	uint32_t i;
//	if (mi == 0)
//		return;
//	if (mi->B) {
//		for (i = 0; i < 1U << mi->b; ++i) {
//			free(mi->B[i].p);
//			mi->B[i].a.clear();
//			mi->B[i].h.clear();
//		}
//	}
//	if (mi->I) {
//		for (i = 0; i < mi->seq.size(); ++i)
//			free(mi->I[i].a);
//		free(mi->I);
//	}
//	free(mi->B);
//	mi->S.clear();
//	free(mi);
//}
//// query sequence name and sequence in the minimap2 index
//int mm_idx_index_name(mm_idx_t *mi);
//int mm_idx_name2id(const mm_idx_t *mi, const char *name);
//int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en,
//		uint8_t *seq);
//
//int mm_idx_alt_read(mm_idx_t *mi, const char *fn);
//int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc);
//int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en,
//		uint8_t *s);
//
//// deprecated APIs for backward compatibility
//void mm_mapopt_init(mm_mapopt_t *opt);
//mm_idx_t* mm_idx_build(const char *fn, int w, int k, int flag, int n_threads);
//
///******************
// * Generate index *
// ******************/
//
//#include <string.h>
//#include <zlib.h>
//
//typedef struct {
//	int mini_batch_size;
//	uint64_t batch_size, sum_len;
//	gzFile fp;
//	mm_idx_t *mi;
//} pipeline_t;
//
//typedef struct {
//	std::vector<mm_bseq1_t> seq;
//	std::vector<mm128_t> a;
//} step_t;
//
//static void mm_idx_add(mm_idx_t *mi, std::vector<mm128_t> &a) {
//	int i, mask = (1 << mi->b) - 1;
//	for (i = 0; i < a.size(); ++i) {
//		std::vector<mm128_t> &p = mi->B[a[i].x >> 8 & mask].a;
//		p.emplace_back(a[i]);
//	}
//}
//
////blank functions
////loading alt strings into memory
//void fa_read_demo(gzFile fp, int64_t chunk_size, std::vector<mm_bseq1_t> &seq_v) {
//	seq_v.clear();
//	int64_t total_load = 0;
//	char *temp = new char[MAX_LINE_LENGTH];//10M
//	char *analysis_line = new char[MAX_LINE_LENGTH];//10M
//	std::vector<std::string> item_value;
//	bool is_the_final_block = false;
//	while(1){
//		if(NULL == gzgets(fp, analysis_line, MAX_LINE_LENGTH)){
//			is_the_final_block = true; break;
//		}
//		else{
//			//remove the final \n
//			analysis_line[strlen(analysis_line) - 1] = 0;//todo:: need check whether to remove
//			if(analysis_line[0] == '>')//skip contig NAME line in .fa format
//				analysis_line[0] = 0;
//			item_value.clear();
//			split_string(item_value, temp, analysis_line, "N");
//			for(uint64_t i = 0; i < item_value.size(); i++){
//				seq_v.emplace_back();
//				seq_v.back().seq = item_value[i];
//				total_load += item_value[i].size();
//			}
//			if(total_load >= chunk_size){ break; }
//		}
//	}
//
//	free(temp);
//	free(analysis_line);
//}
//
//#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))
//#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)
//unsigned char seq_nt4_table[256] = { 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4 };
//
////data SHARED by all threads
//struct MINI_IDX_PRE {
//	int mini_batch_size;
//	uint64_t batch_size, sum_len;
//	gzFile fp;
//	mm_idx_t *mi;
//	std::vector<step_t> s_l;
//	Minimizer_generater mg; //todo INIT::
//};
//
//static void* copy_alt_strings(MINI_IDX_PRE &idx, int tid) {
//	step_t &s = idx.s_l[tid];
//	s.a.clear();
//	s.seq.clear();
//	if (idx.sum_len > idx.batch_size)
//		return 0;
//	fa_read_demo(idx.fp, idx.mini_batch_size, s.seq); // read a mini-batch
//	if (!s.seq.empty()) {
//		uint64_t m = idx.mi->seq.size() + s.seq.size();
//		idx.mi->seq.reserve(m);
//
//		//copy seq in s into seq in p
//		//malloc for seq, 4 bit for per base
//		// make room for idx.mi->S
//		if (!(idx.mi->flag & MM_I_NO_SEQ)) {
//			uint64_t sum_len, max_len;
//			for (int i = 0, sum_len = 0; i < s.seq.size(); ++i)
//				sum_len += s.seq[i].seq.size();
//			max_len = (idx.sum_len + sum_len + 7) / 8;
//			idx.mi->S.reserve(max_len);
//			idx.mi->S.resize(max_len);
//		}
//		// populate idx.mi->seq
//		for (int i = 0; i < s.seq.size(); ++i) {
//			//copy basic info to index
//			idx.mi->seq.emplace_back(idx.mi->flag, s.seq[i], idx.sum_len);
//			mm_idx_seq_t &c_idx_seq = idx.mi->seq.back();
//			// copy the sequence
//			if (!(idx.mi->flag & MM_I_NO_SEQ)) {
//				for (int j = 0; j < c_idx_seq.len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
//					uint64_t o = idx.sum_len + j;
//					int c = seq_nt4_table[(uint8_t) s.seq[i].seq[j]];
//					mm_seq4_set(idx.mi->S, o, c);
//				}
//			}
//			// update idx.sum_len and idx.mi->n_seq
//			idx.sum_len += c_idx_seq.len;
//			s.seq[i].rid = idx.mi->seq.size();
//		}
//	}
//	return 0;
//}
//
//static void* run_mm_sketch(Minimizer_generater &mg, step_t &s) {
//	//void mm_sketch(std::string & seq, uint32_t rid, std::vector<mm128_t> &p);
//	for (int i = 0; i < s.seq.size(); ++i) {
//		mm_bseq1_t &t = s.seq[i];
//		if (!t.seq.empty()) {
//			mg.mm_sketch(t.seq.c_str(),t.seq.size() , s.a);
//		}else
//			fprintf(stderr,
//					"[WARNING] the length database sequence '%s' is 0\n",
//					t.name.c_str());
//	}
//	return 0;
//}
//
////#define MAX_read_size 100000000//100M
////void *classify_pipeline(void *shared, int step, int tid, void *_data) {
////	CLASSIFY_SHARE_DATA * s = (CLASSIFY_SHARE_DATA*) shared;
////	//step0: read read data from files; step1: process; step2: output result
////	if 		(step == 0)	{ 	if((s->data[tid].readNum = load_reads(s->_fp1, s->data[tid].seqs1, s->data[tid].seqs2, N_NEEDED, s->o, s->buff))) 	return (void *)1; }
////	else if (step == 1)	{	kt_for(s->o->thread_n, worker_for, s->data + tid, s->data[tid].readNum);  							return (void *)1; }
////	else if (step == 2)	{	output_results(s->data[tid].readNum, s->idx, s->output_file, s->data[tid].b1, s->data[tid].b2,  s->output_file_ori, s->data[tid].ori_b1, s->data[tid].ori_b2);    	return (void *)1; }
////	return 0;
////}
//
//static void* worker_pipeline(void *shared, int step, int tid, void *_data) {
//	MINI_IDX_PRE *p = (MINI_IDX_PRE*) shared;
//	if (step == 0) {		copy_alt_strings(*p, tid);				return (void*) 1;	} // step 0: read and copy sequences
//	else if (step == 1) {	run_mm_sketch(p->mg, p->s_l[tid]);		return (void*) 1;	} // step 1: compute sketch
//	else if (step == 2) {	mm_idx_add(p->mi, p->s_l[tid].a);		return (void*) 1;	} // dispatch sketch to buckets
//	return 0;
//}
//
//mm_idx_t* mm_idx_init(int w, int k, int b, int flag) {
//	mm_idx_t *mi;
//	if (k * 2 < b)
//		b = k * 2;
//	if (w < 1)
//		w = 1;
//	mi = (mm_idx_t*) calloc(1, sizeof(mm_idx_t));
//	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
//	mi->B = (mm_idx_bucket_t*) calloc(1 << b, sizeof(mm_idx_bucket_t));
//	return mi;
//}
//
///*********************************
// * Sort and generate hash tables *
// *********************************/
//
//static void worker_post(void *g, long i, int tid) {
////	int n, n_keys;
////	size_t j, start_a, start_p;
////
////	mm_idx_t *mi = (mm_idx_t*) g;
////	mm_idx_bucket_t *b = &mi->B[i];
////	std::unordered_map<uint64_t, uint64_t> &h = b->h; //todo:: register ？？？
////	if (b->a.empty())
////		return;
////
////	// sort by minimizer
////	radix_sort_128x(&(b->a[0]), &(b->a[0]) + b->a.size());
////
////	// count and preallocate
////	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.size(); ++j) {
////		if (j == b->a.size() || b->a[j].x >> 8 != b->a[j - 1].x >> 8) {
////			++n_keys;
////			if (n > 1)
////				b->n += n;
////			n = 1;
////		} else
////			++n;
////	}
////	h.reserve(n_keys);
////	b->p = (uint64_t*) calloc(b->n, 8);
////
////	// create the hash table
////	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.size(); ++j) {
////		if (j == b->a.size() || b->a[j].x >> 8 != b->a[j - 1].x >> 8) {
////			mm128_t *p = &(b->a[j - 1]);
////			uint64_t hash_key = p->x >> 8 >> mi->b << 1;
////			//check absence ??
////			auto itr = h.find(hash_key);
////			xassert(itr == h.end() && j == start_a + n, "");
////			if (n == 1) {
////				h[hash_key | 1] = p->y;
////			} else {
////				int k;
////				for (k = 0; k < n; ++k)
////					b->p[start_p + k] = b->a[start_a + k].y;
////				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
////				h[hash_key | 1] = (uint64_t) start_p << 32 | n;
////				start_p += n;
////			}
////			start_a = j, n = 1;
////		} else
////			++n;
////	}
////	b->h = h;
////	assert(b->n == (int32_t )start_p);
////	b->a.clear();
//}
//
//static void mm_idx_post(mm_idx_t *mi, int n_threads) {
//	kt_for(n_threads, worker_post, mi, 1 << mi->b);
//}
//
//mm_idx_t* mm_idx_gen(gzFile fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size) {
//	pipeline_t pl;
//	if (fp == NULL) return 0;
//	memset(&pl, 0, sizeof(pipeline_t));
//	pl.mini_batch_size =
//			(uint64_t) mini_batch_size < batch_size ?
//					mini_batch_size : batch_size;
//	pl.batch_size = batch_size;
//	pl.fp = fp;
//	pl.mi = mm_idx_init(w, k, b, flag);
//
//	kt_pipeline(3, worker_pipeline, &pl, 3);
//	fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__);
//
//	mm_idx_post(pl.mi, n_threads);
//		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__);
//
//	return pl.mi;
//}
//
//mm_idx_t* mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
//		{
//	gzFile fp;
//	mm_idx_t *mi;
//	fp =  gzopen(fn, "rb");
//	if (fp == 0)
//		return 0;
//	mi = mm_idx_gen(fp, w, k, 14, flag, 1 << 18, n_threads, UINT64_MAX);
//	gzclose(fp);
//	return mi;
//}
//
//#endif /* MINIMIZER_IDX_HPP_ */
//
