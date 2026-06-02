/*
 * idx_loader.hpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_PROCESS0_IDX_LOADER_HPP_
#define SRC_BWT_ONLINE_PROCESS0_IDX_LOADER_HPP_

#include "aln_para.hpp"
#include "BWT_idx/var_map.hpp"
#include "BWT_idx/bwt.hpp"
#include "BWT_idx/bwt_cpt_sa.hpp"
#include "index_building/haplotype_online.hpp"

namespace Aln_online{

#define SEED_LEN 25

struct Gap_Aln_SEED_t{
	uint8_t *seed_str_p;
	int seed_len;
	int seed_bg_in_read;
	Gap_Aln_SEED_t(int read_bg, int seed_length, uint8_t* read_string){
		seed_bg_in_read = read_bg;
		seed_len = seed_length;
		seed_str_p = read_string;
	}

	void init(int read_bg, int seed_length, std::vector<uint8_t> & read_string){
		seed_bg_in_read = read_bg;
		seed_len = seed_length;
		seed_str_p = &read_string[read_bg];
	}
};

//47 bit in total, 8 byte
struct Gap_aln_hit_t{
	uint32_t wb_id;//22 bit
	uint32_t suggent_read_to_WB_offset:16;//
	uint32_t seed_id:7;
	uint32_t seed_len:8;
	uint32_t is_rev:1;
	Gap_aln_hit_t(){}

	Gap_aln_hit_t(Gap_Aln_SEED_t &seed, uint8_t ori, uint32_t wb_id, uint32_t wb_offset, uint8_t seed_id, uint WB_step_len){
		//this->seed_p = &seed;
		this->wb_id = wb_id;
		this->is_rev = ori;
		this->seed_len = seed.seed_len;
		if(this->is_rev){
			this->suggent_read_to_WB_offset = wb_offset + seed.seed_bg_in_read;
		}else{
			if((int)wb_offset < seed.seed_bg_in_read){
				this->suggent_read_to_WB_offset = wb_offset + WB_step_len - seed.seed_bg_in_read;
				if(this->wb_id > 0) this->wb_id -= 1;
			}else
				this->suggent_read_to_WB_offset = wb_offset - seed.seed_bg_in_read;
		}
		this->seed_id = seed_id;
	}

	static inline int cmp_by_wb_id(const Gap_aln_hit_t &a, const Gap_aln_hit_t &b){
		if(a.wb_id != b.wb_id)	return a.wb_id < b.wb_id;
		if(a.is_rev != b.is_rev)	return a.is_rev < b.is_rev;
		return a.seed_id < b.seed_id;
	}

};

class IDX_loader{

public:
	BWT_IDX_loader *bwt_idx;
	std::vector<Window_t> var_map;
	BWT_SA_compactor *sa_idx;
	ALN_ONLINE::MM_idx_loader *wb_idx;
public:

	void restore_idx(OL_PAR * p){
		char str[1024];
		//load bwt index
		//bwt_idx = (BWT_IDX_loader*)xmalloc(sizeof(BWT_IDX_loader));
		//bwt_idx->restore_idx(strcat(strcpy(str, p->idx_dir), "bwt_idx/pan_ref.fa"), true, false, true, true);
		//load var map file
		//Window_t::load_variant_map(strcat(strcpy(str, p->idx_dir), "var_map.txt"), 0, 100000000, var_map);
		//load SA index
		//sa_idx = (BWT_SA_compactor*)xmalloc(sizeof(BWT_SA_compactor));
		//sa_idx->load(strcat(strcpy(str, p->idx_dir), "bwt_idx/pan_ref.fa"));//load SA idx
		//load WB index
		wb_idx = (ALN_ONLINE::MM_idx_loader*)xcalloc(1, sizeof(ALN_ONLINE::MM_idx_loader));
		wb_idx->load_all_index(strcat(strcpy(str, p->idx_dir), "wb_index/"), NULL, true);
		//for sa benchmark test
		// load simdata read pos info
	}

	void seq_bin(int len, char * seq, ubyte_t *seq_bin, int is_rev){

		
			for (int i = 0; i < len; ++i) {
				//seq_bin[i] = nst_nt4_table[seq[i]];
				switch(seq[i]){
				case 'A': seq_bin[i] = 0; break;
				case 'C': seq_bin[i] = 1; break;
				case 'G': seq_bin[i] = 2; break;
				case 'T':  seq_bin[i] = 3; break;
				default :
					 seq_bin[i] = 4; break;
				}
			}
			if(is_rev){
				for (int i = 0; i < len>>1; ++i) {
	      			char tmp = seq_bin[len-1-i];
	      			seq_bin[len-1-i] = seq_bin[i]; seq_bin[i] = tmp;
	   			 }
			}
	}

};

}

#endif /* SRC_BWT_ONLINE_PROCESS0_IDX_LOADER_HPP_ */
