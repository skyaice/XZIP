/*
 * bwt.hpp
 *
 *  Created on: 2022年5月26日
 *      Author: fenghe
 */

#ifndef BWT_HPP_
#define BWT_HPP_

#include <vector>
#include <string.h>
#include <stdlib.h>
extern "C"{
	#include "../clib/utils.h"
	//#include "../clib/bam_file.h"
	#include "../clib/desc.h"
}
#include "../htslib/htslib/sam.h"
typedef uint64_t bwtint_t;
typedef uint8_t ubyte_t;

// requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

// The following two lines are ONLY correct when OCC_INTERVAL==0x80
#define bwt_bwt(k) (bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
//inline bwtint_t bwt_bwt(bwtint_t k){ return (bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)]);}

//inline uint32_t * bwt_occ_intv(bwtint_t k){ return (bwt + ((k)>>7<<4)); }
#define bwt_occ_intv(k) (bwt + ((k)>>7<<4))
/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(k) (bwt_bwt(k)>>((~(k)&0xf)<<1)&3)
//inline bwtint_t bwt_B0(bwtint_t k){	return (bwt_bwt(k)>>((~(k)&0xf)<<1)&3);	}
#define bwt_set_intv(c, ik) ((ik).x[0] = L2[(int)(c)]+1, (ik).x[2] = L2[(int)(c)+1]-L2[(int)(c)], (ik).x[1] = L2[3-(c)]+1, (ik).info = 0)

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	int32_t is_alt;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

struct bntseq_t{
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;

	inline int64_t sa_depos(uint64_t pos, bool &is_rev)
	{
		is_rev = (pos >= (l_pac));
		return (is_rev)? (l_pac*2) - 1 - pos : pos;
	}

	int bns_pos2rid_buff(int64_t pos_f, uint32_t & offset, uint32_t &old_chrID)
	{
		if(pos_f >= anns[old_chrID].offset && pos_f < anns[old_chrID + 1].offset){
			offset = pos_f - anns[old_chrID].offset;
			return old_chrID;
		}
		int chr_ID = bns_pos2rid(pos_f, offset);
		if(chr_ID > 0) old_chrID = chr_ID;
		return chr_ID;
	}

	int bns_pos2rid(int64_t pos_f, uint32_t & offset)
	{
		int left, mid, right; offset = 0;
		if (pos_f >= l_pac) return -1;
		left = 0; mid = 0; right = n_seqs;
		while (left < right) { // binary search
			mid = (left + right) >> 1;
			if (pos_f >= anns[mid].offset) {
				if (mid == n_seqs - 1) break;
				if (pos_f < anns[mid+1].offset) break; // bracketed
				left = mid + 1;
			} else right = mid;
		}
		offset = pos_f - anns[mid].offset;
		return mid;
	}

};

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename);
bntseq_t *bns_restore(const char *prefix);
void bns_destroy(bntseq_t *bns);
void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line, bam_hdr_t *header, char **sn);
//char *bwa_set_rg(const char *sample_name, const char *readset, const char *pl,const char *cn);
char *bwa_set_rg(const char *s);
char *bwa_insert_header(const char *s, char *hdr);

class BWT_IDX_loader{

public:
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;

	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

	void bwt_gen_cnt_table()
	{
		int i, j;
		for (i = 0; i != 256; ++i) {
			uint32_t x = 0;
			for (j = 0; j != 4; ++j)
				x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
			cnt_table[i] = x;
		}
	}

	static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
	{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
		const int bufsize = 0x1000000; // 16M block
		bwtint_t offset = 0;
		while (size) {
			int x = bufsize < size? bufsize : size;
			if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
			size -= x; offset += x;
		}
		return offset;
	}

public:
	uint8_t *get_pac(){return pac;}
	uint64_t get_seq_len(){return seq_len;}
	inline bool is_seq_len(bwtint_t i){return seq_len == i;};

	bntseq_t *bns;

	void restore_idx(const char *prefix, bool load_BWT, bool load_SA, bool load_BNS, bool load_PAC){
		char str[1024];
		if(load_BWT){
			strcpy(str, prefix); strcat(str, ".bwt");  bwt_restore_bwt(str);
		}
		if(load_SA){
			strcpy(str, prefix); strcat(str, ".sa");   bwt_restore_sa(str);
		}
		if(load_BNS){
			bns = bns_restore(prefix);
		}
		if(load_BNS && load_PAC){
			pac = (uint8_t *)xcalloc(bns->l_pac/4+1, 1);
			err_fread_noeof(pac, 1, bns->l_pac/4+1, bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(bns->fp_pac);
		}
	}

public:
	void bwt_restore_bwt(const char *fn){
		FILE *fp = xopen(fn, "rb");
		err_fseek(fp, 0, SEEK_END);
		bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
		bwt = (uint32_t*)xcalloc(bwt_size, 4);
		err_fseek(fp, 0, SEEK_SET);
		err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
		L2[0] = 0;
		err_fread_noeof(L2+1, sizeof(bwtint_t), 4, fp);
		fread_fix(fp, bwt_size<<2, bwt);
		seq_len = L2[4];
		err_fclose(fp);
		bwt_gen_cnt_table();
	}

	void bwt_restore_sa(const char *fn)
	{
		char skipped[256];
		FILE *fp;
		bwtint_t primary;

		fp = xopen(fn, "rb");
		err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
		xassert(primary == primary, "SA-BWT inconsistency: primary is not the same.");
		err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
		err_fread_noeof(&sa_intv, sizeof(bwtint_t), 1, fp);
		err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
		xassert(primary == seq_len, "SA-BWT inconsistency: seq_len is not the same.");

		n_sa = (seq_len + sa_intv) / sa_intv;
		sa = (bwtint_t*)xcalloc(n_sa, sizeof(bwtint_t));
		sa[0] = -1;

		fread_fix(fp, sizeof(bwtint_t) * (n_sa - 1), sa + 1);
		err_fclose(fp);
	}

	void bwt_destroy(){
		if(sa != NULL) { free(sa);  sa = NULL; }
		if(bwt != NULL){ free(bwt); bwt = NULL; }
		if(pac != NULL){ free(pac); pac = NULL; }
	}

	static inline int __occ_aux(uint64_t y, int c)
	{
		// reduce nucleotide counting to bits counting
		y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
		// count the number of 1s in y
		y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
		return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
	}
	bntann1_t *get_ann(int i){
		bntann1_t *p = bns->anns + i;
		return p;
	}

	bntann1_t *get_ann_ref(int i){
		bntann1_t *p = bns->anns + bns->n_seqs - 25 + i;
		return p;
	}

	bwtint_t bwt_occ(bwtint_t k, ubyte_t c){
		bwtint_t n;
		uint32_t *p, *end;

		if (k == seq_len) return L2[c+1] - L2[c];
		if (k == (bwtint_t)(-1)) return 0;
		k -= (k >= primary); // because $ is not in bwt

		// retrieve Occ at k/OCC_INTERVAL
		n = ((bwtint_t*)(p = bwt_occ_intv(k)))[c];
		p += sizeof(bwtint_t); // jump to the start of the first BWT cell

		// calculate Occ up to the last k/32
		end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
		for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

		// calculate Occ
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		return n;
	}

	// an analogy to bwt_occ() but more efficient, requiring k <= l
	void bwt_2occ(bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
	{
		bwtint_t _k, _l;
		_k = (k >= primary)? k-1 : k;
		_l = (l >= primary)? l-1 : l;
		if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
			*ok = bwt_occ(k, c);
			*ol = bwt_occ(l, c);
		} else {
			bwtint_t m, n, i, j;
			uint32_t *p;
			if (k >= primary) --k;
			if (l >= primary) --l;
			n = ((bwtint_t*)(p = bwt_occ_intv(k)))[c];
			p += sizeof(bwtint_t);
			// calculate *ok
			j = k >> 5 << 5;
			for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
				n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
			m = n;
			n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
			if (c == 0) n -= ~k&31; // corrected for the masked bits
			*ok = n;
			// calculate *ol
			j = l >> 5 << 5;
			for (; i < j; i += 32, p += 2)
				m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
			m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
			if (c == 0) m -= ~l&31; // corrected for the masked bits
			*ol = m;
		}
	}

	#define __occ_aux4(b)								\
		(cnt_table[(b)&0xff] + cnt_table[(b)>>8&0xff]		\
		 + cnt_table[(b)>>16&0xff] + cnt_table[(b)>>24])

	void bwt_occ4(bwtint_t k, bwtint_t cnt[4])
	{
		bwtint_t x;
		uint32_t *p, tmp, *end;
		if (k == (bwtint_t)(-1)) {
			memset(cnt, 0, 4 * sizeof(bwtint_t));
			return;
		}
		k -= (k >= primary); // because $ is not in bwt
		p = bwt_occ_intv(k);
		memcpy(cnt, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
		for (x = 0; p < end; ++p) x += __occ_aux4(*p);
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(tmp) - (~k&15);
		cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
	}

	// an analogy to bwt_occ4() but more efficient, requiring k <= l
	void bwt_2occ4(bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
	{
		bwtint_t _k, _l;
		_k = k - (k >= primary);
		_l = l - (l >= primary);
		if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
			bwt_occ4(k, cntk);
			bwt_occ4(l, cntl);
		} else {
			bwtint_t x, y;
			uint32_t *p, tmp, *endk, *endl;
			k -= (k >= primary); // because $ is not in bwt
			l -= (l >= primary);
			p = bwt_occ_intv(k);
			memcpy(cntk, p, 4 * sizeof(bwtint_t));
			p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
			// prepare cntk[]
			endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
			endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
			for (x = 0; p < endk; ++p) x += __occ_aux4(*p);
			y = x;
			tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
			x += __occ_aux4(tmp) - (~k&15);
			// calculate cntl[] and finalize cntk[]
			for (; p < endl; ++p) y += __occ_aux4(*p);
			tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
			y += __occ_aux4(tmp) - (~l&15);
			memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
			cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
			cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
		}
	}

	inline bwtint_t bwt_invPsi(bwtint_t k) // compute inverse CSA
	{
		bwtint_t x = k - (k > primary);
		x = bwt_B0(x);
		x = L2[x] + bwt_occ(k, x);
		return k == primary? 0 : x;
	}

	inline bwtint_t debug_bwt_invPsi(bwtint_t k) // compute inverse CSA
	{
		bwtint_t x = k - (k > primary);
		x = bwt_B0(x);
		fprintf(stderr, "%c", "ACGT"[x]);
		x = L2[x] + bwt_occ(k, x);
		return k == primary? 0 : x;
	}

//	bwtint_t bwt_sa_region(bwtint_t k, bwtint_t l)
//	{
//		bwtint_t sa_int = 0, mask = sa_intv - 1;
//		while (k & mask) {
//			++sa_int;
//
//			bwtint_t x = k - (k > primary);
//
//			x = ((bwt[((x)>>7<<4) + sizeof(bwtint_t) + (((x)&0x7f)>>4)])>>((~(x)&0xf)<<1)&3);
//
//			//x = bwt_B0(x);
//
//			x = L2[x] + bwt_occ(k, x);
//			k= (k == primary)? 0 : x;
//
//		}
//		/* without setting bwt->sa[0] = -1, the following line should be
//		   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
//		return sa_int + sa[k/sa_intv];
//	}

// private:
// 	std::vector<uint64_t> tree_bg_stack;
// 	std::vector<uint64_t> tree_ed_stack;
// 	std::vector<uint16_t> tree_bwt_depth;
// 	std::vector<char> tree_char_stack;
// 	int tree_max_length;
// public:
// 	std::vector<char> tree_string;
// 	void bwt_tree_init(int max_length_){
// 		tree_bg_stack.clear();
// 		tree_ed_stack.clear();
// 		tree_bwt_depth.clear();

// 		tree_char_stack.clear();
// 		tree_char_stack.emplace_back('$');

// 		tree_string.clear();
// 		tree_string.emplace_back(0);
// 		tree_bg_stack.emplace_back(0);
// 		tree_ed_stack.emplace_back(seq_len);

// 		//tree_bg_stack.emplace_back(1);
// 		//tree_ed_stack.emplace_back(495980);


// 		tree_bwt_depth.emplace_back(0);
// 		tree_max_length = max_length_;
// 	}

// 	void bwt_tree_next(uint64_t & bwt_bg, uint64_t & bwt_ed){
// 		tree_string.pop_back();//remove 0:
// 		if(tree_bwt_depth.back() < tree_max_length){
// 			bwtint_t cntk[4]; bwtint_t cntl[4];
// 			while(!tree_bwt_depth.empty()){
// 				uint16_t father_bwt_depth = tree_bwt_depth.back();
// 				tree_bwt_depth.pop_back();
// 				uint64_t father_bwt_bg = tree_bg_stack.back();
// 				tree_bg_stack.pop_back();
// 				uint64_t father_bwt_ed = tree_ed_stack.back();
// 				tree_ed_stack.pop_back();

// 				if(tree_string.size() > father_bwt_depth){
// 					tree_string.resize(father_bwt_depth);
// 				}

// 				tree_string.emplace_back(tree_char_stack.back());
// 				tree_char_stack.pop_back();

// 				bwt_2occ4(father_bwt_bg - 1, father_bwt_ed, cntk, cntl);
// 				for(int i = 3; i >= 0; i--){
// 					if(cntk[i] < cntl[i]){//with rst
// 						tree_bwt_depth.emplace_back(father_bwt_depth + 1);
// 						tree_bg_stack.emplace_back(L2[i] + cntk[i] + 1);
// 						tree_ed_stack.emplace_back(L2[i] + cntl[i]);
// 						tree_char_stack.emplace_back("ACGT"[i]);
// 					}
// 				}
// 				if(father_bwt_depth + 1 == tree_max_length){
// 					break;
// 				}
// 			}
// 		}

// 		if(tree_bwt_depth.empty()){
// 			bwt_bg = bwt_ed = MAX_uint64_t;
// 		}
// 		else{
// 			bwt_bg = tree_bg_stack.back();
// 			tree_bg_stack.pop_back();
// 			bwt_ed = tree_ed_stack.back();
// 			tree_ed_stack.pop_back();
// 			tree_bwt_depth.pop_back();
// 			if(tree_string.size() == tree_max_length + 1){
// 				tree_string.pop_back();
// 			}
// 			//tree_string.pop_back();
// 			tree_string.emplace_back(tree_char_stack.back());
// 			tree_char_stack.pop_back();
// 		}
// 		tree_string.emplace_back(0);
// 	}

#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
	uint8_t *bns_get_seq(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len)
	{
		uint8_t *seq = 0;
		if (end < beg) end ^= beg, beg ^= end, end ^= beg; // if end is smaller, swap
		if (end > l_pac<<1) end = l_pac<<1;
		if (beg < 0) beg = 0;
		if (beg >= l_pac || end <= l_pac) {
			int64_t k, l = 0;
			*len = end - beg;
			seq = (uint8_t *)xmalloc(end - beg);
			if (beg >= l_pac) { // reverse strand
				int64_t beg_f = (l_pac<<1) - 1 - end;
				int64_t end_f = (l_pac<<1) - 1 - beg;
				for (k = end_f; k > beg_f; --k)
					seq[l++] = 3 - _get_pac(pac, k);
			} else { // forward strand
				for (k = beg; k < end; ++k)
					seq[l++] = _get_pac(pac, k);
			}
		} else *len = 0; // if bridging the forward-reverse boundary, return nothing
		return seq;
	}

	void bwt_same_check(){
		//rst:
		uint8_t seq1[300]; memset(seq1,0,300);
		uint8_t seq2[300]; memset(seq2,0,300);
		uint8_t* str = seq1;
		uint8_t* old_str = seq2;
		int str_same_count = 0;
		for(uint64_t i = 0; i < seq_len; i++){
			std::swap(str, old_str);
			bool is_rev;
			uint64_t golbal_offset = bwt_sa(i);
			//uint64_t golbal_offset = random() % seq_len;
			golbal_offset = sa_depos(golbal_offset, is_rev);
			//show string:
			int64_t l = 0;
			if(is_rev){
				int64_t beg = golbal_offset -150;
				int64_t end = golbal_offset;
				for (int64_t k = end; k > beg; --k)
					str[l++] = 3 - _get_pac(pac, k);
			}else{
				int64_t beg = golbal_offset;
				int64_t end = golbal_offset + 150;
				for (int64_t k = beg; k < end; ++k)
					str[l++] = _get_pac(pac, k);
			}
			//compare
			bool str_same = true;
			for(int i = 150 - 1; i >= 0; i--){
				if(str[i] != old_str[i]){
					str_same = false;
					break;
				}
			}
			if(str_same == false){//when the string is not same
				if(false){
					for(int i = 0;i < 150; i++){
						fprintf(stderr, "%c", "ACGT"[old_str[i]]);
					}
					fprintf(stderr, "golbal_offset %ld str_same %d, str_same_count %d\n", golbal_offset, str_same, str_same_count);
				}
				//reset count
				str_same_count = 1;
			}
			else //when the string is same
			{
				str_same_count ++;



			}



//
//			for(int i = 0;i < 150; i++){
//				fprintf(stderr, "%c", "ACGT"[str[i]]);
//			}

//			//
			if(i % 1000000 == 0){
				fprintf(stderr, "i %ld golbal_offset %ld is_rev %d\n", i, golbal_offset, is_rev);
			}
		}


	}

	inline int64_t sa_depos(uint64_t pos, bool &is_rev)
	{
		is_rev = (pos >= (seq_len>>1));
		return (is_rev)? (seq_len) - 1 - pos : pos;
	}


	bwtint_t bwt_sa(bwtint_t k)
	{
		bwtint_t sa_int = 0, mask = sa_intv - 1;
		while (k & mask) {
			++sa_int;
			k = bwt_invPsi(k);
		}
		/* without setting bwt->sa[0] = -1, the following line should be
		   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
		return sa_int + sa[k/sa_intv];
	}

	int bwt_match_exact_demo(int len, const ubyte_t *str, bwtint_t &sa_begin, bwtint_t &sa_end)
	{
		bwtint_t k, l, ok, ol;
		int i;
		k = 0; l = seq_len;
		for (i = len - 1; i >= 0; --i) {
			fprintf(stderr, "i %d k %ld l %ld \n", i, k, l);
			ubyte_t c = str[i];
			if (c > 3) return 0; // no match
			bwt_2occ(k - 1, l, c, &ok, &ol);
			k = L2[c] + ok + 1;
			l = L2[c] + ol;
			if (k > l) break; // no match
		}
		if (k > l) return 0; // no match
		sa_begin = k;
		sa_end = l;
		return l - k + 1;
	}

	int bwt_match_exact_forward_demo(int len, const ubyte_t *str, bwtint_t &sa_begin, bwtint_t &sa_end, int &last_base)
	{
		bwtint_t k, l, ok, ol;
		int i;
		k = 0; l = seq_len;
		for (i = 0; i < len; i++) {
			ubyte_t c = str[i];
			if (c > 3){last_base = i;return 0;}  // no match
			c = 3-c;
			bwt_2occ(k - 1, l, c, &ok, &ol);
			k = L2[c] + ok + 1;
			l = L2[c] + ol;
			if (k > l) break; // no match
		}
		last_base = i;
		if (k > l) return 0; // no match
		sa_begin = k;
		sa_end = l;


		return l - k + 1;
	}

	int bwt_compact_seq(int len, const ubyte_t *str, const char *qual, bwtint_t &sa_begin, bwtint_t &sa_end, int &early_break, uint8_t min_base_qual, int read_end_id, bool early_break_flag, int &path_N)
	{
		//double cpu_time = cputime(); double rtime = realtime();	
		bwtint_t k, l, ok, ol;
		bwtint_t pre_k, pre_l;
		int i;
		k = 0; l = seq_len;
		//for (i = len - 1; i >= 0; --i) {
		for (i = 0; i< len; ++i){
			bwtint_t old_k = k, old_l = l;
			ubyte_t c = str[i];
			if(c<=3){
				c = 3-c;
				bwt_2occ(k - 1, l, c, &ok, &ol);
				k = L2[c] + ok + 1;
				l = L2[c] + ol;
				c = 3-c;
				
				
			}
			if(c > 3){
				if(qual[i] - '#' <= min_base_qual){
					bwtint_t cntk[4]; bwtint_t cntl[4];
					bwt_2occ4(old_k - 1, old_l, cntk, cntl);
					path_N = 0;
					for(int j = 0; j < 4; j++){
						if(cntk[j] < cntl[j]){
							path_N++;
							k = L2[j] + cntk[j] + 1;
							l = L2[j] + cntl[j];
						}
					}
					if(path_N != 1){ k = 1; l = 0; break; }
				}
				else{ 
					break; 
				}
			}
			if(k <= l){
				pre_k = k;
				pre_l = l;
				early_break = i;
			}else{
				break;
				}
		}
		sa_begin = pre_k;
		sa_end = pre_l;
		//fprintf(stderr, "inside %d bwt time  : cpu time %.3f sec, real time %.3f sec\n" ,k<l, cputime() - cpu_time, realtime()-rtime);
		return l - k + 1;
	}

	int calculate_lowq(int len, const char *str, const char *qual,  uint8_t min_base_qual, int read_end_id, uint64_t count[4][150], uint64_t read_qual_count[150]){
		int i = 0;
		int lowq_count = 0;
		char qual_pre = 'T';
		int run_length_count = 1;
		for (i = 0; i< len; ++i){
			if(qual[i]!=qual_pre){
				if(qual_pre !='T'){
					read_qual_count[run_length_count] ++;
					run_length_count = 1;
				}
			}else{
				run_length_count ++;
			}
			qual_pre = qual[i];
			if(i==len-1){
				// if(run_length_count == 150){
				// 	fprintf(stderr, "a");
				// }
				if(qual[i]!=qual_pre) 
					read_qual_count[1]++;
				else 
					read_qual_count[run_length_count]++;
			}
			switch (qual[i])
			{
			case '#':

				count[0][i] += 1;
				lowq_count += 1;
				break;
			case ',':
				count[1][i] += 1;
				lowq_count += 1;
				break;
			case ':':
				count[2][i] += 1;
				lowq_count +=1;
				break;
			case 'F':
				count[3][i] += 1;
				
				break;
			
			default:
				exit(1);
				break;
			}
		}
		// if(lowq_count>=50) 
		// 	lowq_count = 49;
		// if(lowq_count > 50)
		// 	fprintf(stderr, "");
		//read_qual_count[lowq_count] += 1;

	}

	int bwt_match_exact_skip_low_mapq_base(int len, const ubyte_t *str, const char *qual, bwtint_t &sa_begin, bwtint_t &sa_end, int &early_break, uint8_t min_base_qual, int read_end_id, bool early_break_flag, int &path_N)
	{
		//double cpu_time = cputime(); double rtime = realtime();	
		bwtint_t k, l, ok, ol;
		int i;
		k = 0; l = seq_len;
		//for (i = len - 1; i >= 0; --i) {
		for (i = 0; i< len; ++i){
			bwtint_t old_k = k, old_l = l;
			ubyte_t c = str[i];
			if(c<=3){
				c = 3-c;
				bwt_2occ(k - 1, l, c, &ok, &ol);
				k = L2[c] + ok + 1;
				l = L2[c] + ol;
				c = 3-c;
			}
			if(early_break_flag && k == l ){
				if(i<80){
					early_break = i+1;
					break;
				}
				
			}
			if(k > l || c > 3){
				if(qual[i] - '#' <= min_base_qual){
					bwtint_t cntk[4]; bwtint_t cntl[4];
					bwt_2occ4(old_k - 1, old_l, cntk, cntl);
					path_N = 0;
					for(int j = 0; j < 4; j++){
						if(cntk[j] < cntl[j]){
							path_N++;
							k = L2[j] + cntk[j] + 1;
							l = L2[j] + cntl[j];
						}
					}
					if(path_N != 1){ k = 1; l = 0; break; }
				}
				else{ 
					k = 1; l = 0; 
					// code for stats
					
					break; 
				}
			}
		}
		sa_begin = k;
		sa_end = l;
		//fprintf(stderr, "inside %d bwt time  : cpu time %.3f sec, real time %.3f sec\n" ,k<l, cputime() - cpu_time, realtime()-rtime);
		return l - k + 1;
	}


	int bwt_match_exact_forward_skip_low_mapq_base(int len, const ubyte_t *str, const char *qual, bwtint_t &sa_begin, bwtint_t &sa_end, uint8_t min_base_qual)
	{
		bwtint_t k, l, ok, ol;
		int i;
 		k = 0; l = seq_len;
		for (i = 0; i < len; i++) {
			bwtint_t old_k = k, old_l = l;
			ubyte_t c = str[i];
			if(c <= 3){
				c = 3-c;
				bwt_2occ(k - 1, l, c, &ok, &ol);
				k = L2[c] + ok + 1;
				l = L2[c] + ol;
				c = 3-c;
			}
			if(k > l || c > 3){
				if(qual[i] - '#' <= min_base_qual){
					bwtint_t cntk[4]; bwtint_t cntl[4];
					bwt_2occ4(old_k - 1, old_l, cntk, cntl);
					int path_N = 0;
					for(int j = 0; j < 4; j++){
						if(cntk[j] < cntl[j]){
							path_N++;
							k = L2[j] + cntk[j] + 1;
							l = L2[j] + cntl[j];
						}
					}
					if(path_N != 1){ k = 1; l = 0; break; }
				}
				else{ k = 1; l = 0;	break; }
			}
		}
		sa_begin = k;
		sa_end = l;
		return l - k + 1;
	}

	//deBUG code:
	int bwt_match_exact_forward_skip_low_mapq_base_show_error(int len, const ubyte_t *str, const char *qual, bwtint_t &sa_begin, bwtint_t &sa_end, uint8_t min_base_qual)
		{
			bwtint_t k, l, ok, ol;
			int i;
			k = 0; l = seq_len;
			for (i = 0; i < len; i++) {
				bwtint_t old_k = k, old_l = l;
				ubyte_t c = str[i];
				if(c <= 3){
					c = 3-c;
					bwt_2occ(k - 1, l, c, &ok, &ol);
					k = L2[c] + ok + 1;
					l = L2[c] + ol;
					c = 3-c;
				}
				if(k > l || c > 3){
					if(qual[i] - '#' <= min_base_qual){
						bwtint_t cntk[4]; bwtint_t cntl[4];
						bwt_2occ4(old_k - 1, old_l, cntk, cntl);
						int path_N = 0;
						for(int j = 0; j < 4; j++){
							if(cntk[j] < cntl[j]){
								path_N++;
								k = L2[j] + cntk[j] + 1;
								l = L2[j] + cntl[j];
							}
						}
						if(path_N != 1){ k = 1; l = 0; break; }
						fprintf(stderr, "BASE change： str[i] %c qual[i] %c \n", "ACGT"[str[i]], qual[i]);
					}
					else{ k = 1; l = 0;	break; }
				}
			}
			sa_begin = k;
			sa_end = l;
			return l - k + 1;
		}


	//with rst: l-k+1 > 0
	//without rst: l-k+1 <= 0
	//without rst: l+1 <= k

	//k = L2[c] + ok + 1;
	//l = L2[c] + ol;

	//without rst:  ol - (ok) <= 0
	//with rst:  ol - (ok) > 0


	int bwt_match_exact_forward_skip_low_mapq_base_debug(int len, const ubyte_t *str, const char *qual, bwtint_t &sa_begin, bwtint_t &sa_end, uint8_t min_base_qual, int &last_base, uint8_t &stop_reason)
	{
		stop_reason = 1;
		bwtint_t k, l, ok, ol;
		int i;
		k = 0; l = seq_len;
		for (i = 0; i < len; i++) {
			bwtint_t old_k = k, old_l = l;
			ubyte_t c = str[i];
			if(c <= 3){
				c = 3-c;
				bwt_2occ(k - 1, l, c, &ok, &ol);
				k = L2[c] + ok + 1;
				l = L2[c] + ol;
				c = 3-c;
			}
			if(k > l || c > 3){
				if(qual[i] - '#' <= min_base_qual){
					bwtint_t cntk[4]; bwtint_t cntl[4];
					bwt_2occ4(old_k - 1, old_l, cntk, cntl);
					int path_N = 0;
					for(int j = 0; j < 4; j++){
						if(cntk[j] < cntl[j]){//with rst:
							path_N++;
							k = L2[j] + cntk[j] + 1;
							l = L2[j] + cntl[j];
						}
					}
					if(path_N != 1){ k = 1; l = 0; stop_reason = path_N; break; }
				}
				else{ k = 1; l = 0;	break; }
			}
		}
		last_base = i;
		if (k > l) {if(stop_reason == 1)stop_reason = 0; return 0;} // no match
		sa_begin = k;
		sa_end = l;
		return l - k + 1;
	}

	void seq_reverse(int len, ubyte_t *seq, int is_comp)
	{
	  int i;
	  if (is_comp) {
	    for (i = 0; i < len>>1; ++i) {
	      char tmp = seq[len-1-i];
	      if (tmp < 4) tmp = 3 - tmp;
	      seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	      seq[i] = tmp;
	    }
	    if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	  } else {
	    for (i = 0; i < len>>1; ++i) {
	      char tmp = seq[len-1-i];
	      seq[len-1-i] = seq[i]; seq[i] = tmp;
	    }
	  }
	}

	void seq_bin(int len, char * seq, ubyte_t *seq_bin){
		for (int i = 0; i < len; ++i) {
			switch(seq[i]){
			case 'A': seq_bin[i] = 0; break;
			case 'C': seq_bin[i] = 1; break;
			case 'G': seq_bin[i] = 2; break;
			case 'T':  seq_bin[i] = 3; break;
			default :
				 seq_bin[i] = 4; break;
			}
		}
	}


private:

};



#endif /* BWT_HPP_ */
