/*
 * bwt_cpt_sa.hpp
 *
 *  Created on: 2022年6月15日
 *      Author: fenghe
 */

#ifndef BWT_CPT_SA_HPP_
#define BWT_CPT_SA_HPP_

extern "C"{
	#include "../clib/utils.h"
}
#include "bwt.hpp"
#include "../index_building/haplotype_online.hpp"
#include "../CPPLIB/tools.hpp"
#include <vector>
#include <algorithm>

#define SA_IDX_OFFSET 8
#define SA_IDX_MASK 0xff

//for WB and POS
#define POS_OFFSET 22
#define POS_MASK 0x3fffff
#define WB_POS_MAX ((POS_MASK)-0x3ff)
#define REF_POS_MAX (POS_MASK)

struct BASIC_MAP_rst{
	int32_t wb_ID;
	//POS: for exact alignment, this value is the position of read in WB
	//POS: for gap alignment, this value is -1
	int32_t aln_position_in_WB;
	int8_t is_rev;
	int8_t is_gap_aln;
	//SCORE: for exact alignment, this value is the read length,
	//for gap alignment, this value is the length of all seed/hit
	int16_t single_read_score;

	void set_exact_result(int32_t wb_ID, int32_t POS, int8_t ori, int read_length){
		this->wb_ID = wb_ID;
		this->aln_position_in_WB = POS;
		this->is_rev = ori;
		if(POS < 0 || POS > 0x3fffff){ fprintf(stderr, "FATAL ERROR 1: wb_ID %d, POS %d\n", wb_ID, aln_position_in_WB); xassert(0,"");}
		this->single_read_score = read_length;//todo:: re-set, the score of exact alignment is not always the read length
		this->is_gap_aln = false;
	}

	void set_gap_result(int32_t wb_ID, int8_t ori, int cluster_score, int aln_position_in_WB){
		this->wb_ID = wb_ID;
		this->is_rev = ori;
		this->single_read_score = cluster_score;
		this->is_gap_aln = true;
		this->aln_position_in_WB = aln_position_in_WB;
	}

	void set_unmapped_result(int32_t wb_ID, int8_t ori){
		this->wb_ID = wb_ID;
		this->is_rev = ori;
		this->single_read_score = 0;
		this->is_gap_aln = true;
		this->aln_position_in_WB = -1;
	}

	bool is_unmapped_result(){
		if(this->aln_position_in_WB==-1){
			return true;
		}
		return false;
	}

	BASIC_MAP_rst(){
		wb_ID = -1;
		aln_position_in_WB = -1;
		is_rev = -1;
		single_read_score = -1;
		is_gap_aln = 0;
	}

	static inline int cmp_by_wb_id(const BASIC_MAP_rst &a, const BASIC_MAP_rst &b){
		return a.wb_ID < b.wb_ID;
	}
	//sort:
	static inline int cmp_by_single_read_score(const BASIC_MAP_rst &a, const BASIC_MAP_rst &b){
		return a.single_read_score > b.single_read_score;
	}

	void print(){ fprintf(stderr, "chr_ID %d, POS %d, ori %d single_read_score %d, is_gap_aln %d; "
			, wb_ID, aln_position_in_WB, is_rev, single_read_score, is_gap_aln); }
//	void set(std::vector<BASIC_MAP_rst> & l, int idx){
//		if(idx < 0){
//			wb_ID = -1;
//			exact_aln_position_in_WB = -1;
//			is_rev = -1;
//		}else{
//			wb_ID = l[idx].wb_ID;
//			exact_aln_position_in_WB = l[idx].exact_aln_position_in_WB;
//			is_rev = l[idx].is_rev;
//		}
//	}
	int inline pair_score(BASIC_MAP_rst &r2, int SV_wb_bg){
		int wb_diff = wb_ID - r2.wb_ID;
		//SV check
		int score = 16 + r2.single_read_score + single_read_score;
		if(wb_ID >= SV_wb_bg){
			score -= 8;
		}
		if(r2.wb_ID >= SV_wb_bg){
			score -= 8;
		}
		//POSITION check
		if(wb_diff <= 6 && wb_diff >= -6){
			score += 80;
		}
		//direction check
		if(is_rev != r2.is_rev){
			if((!is_rev && wb_diff <= 0) || (is_rev && wb_diff >= 0)){
				score += 40;
			}
		}
		return score;
	}

	int inline pair_score(int SV_wb_bg){
		if(wb_ID >= SV_wb_bg){
			return single_read_score;
		}
		return 8 + single_read_score;
	}
};

struct BWT_MAP_rst_CPT_SA{
	int64_t bwt_idx;
	uint32_t wb_ID;
	int32_t POS;
	int8_t is_rev;
	int8_t is_REF;

	BWT_MAP_rst_CPT_SA(int64_t bwt_idx_, uint32_t wb_ID_, uint32_t POS_, int8_t ori_, int8_t is_REF_){
		bwt_idx = bwt_idx_;
		wb_ID = wb_ID_;
		POS = POS_;
		is_rev = ori_;
		is_REF = is_REF_;
		if(POS_ < 0 || POS_ > 0x3fffff){ fprintf(stderr, "FATAL ERROR 1: wb_ID %d, POS %d\n", wb_ID, POS); xassert(0,"");}
	}

	static inline int cmp_by_position(const BWT_MAP_rst_CPT_SA &a, const BWT_MAP_rst_CPT_SA &b){
		if(a.wb_ID == b.wb_ID)	return a.POS > b.POS;
		else							return a.wb_ID < b.wb_ID;
	}

	static inline int cmp_by_bwt_idx(const BWT_MAP_rst_CPT_SA &a, const BWT_MAP_rst_CPT_SA &b){
		return a.bwt_idx < b.bwt_idx;
	}

	void print(){ fprintf(stderr, "wb_ID %d, POS %d, ori %d is_REF %d\n", wb_ID, POS, is_rev, is_REF); }
};

struct BWT_same_Checker{
	struct bwt_same_SA{
		int64_t bwt;
		int64_t golbal_offset;
		bool is_rev;
		bwt_same_SA(int64_t bwt_, int64_t golbal_offset_,bool is_rev_){ bwt = bwt_; golbal_offset = golbal_offset_; is_rev = is_rev_;}
		void show(){	fprintf(stderr, "bwt %ld, golbal_offset  %ld, is_rev %d\n", bwt, golbal_offset, is_rev);}
	};

	void show_c_same_bwt(){
		if(true){
			for(int i = 0;i < search_length; i++){ fprintf(stderr, "%c", "ACGT"[old_str[i]]); }
			fprintf(stderr, "golbal_offset %ld, str_same_count %ld\n", golbal_offset, rst.size());
//			for(auto & r: rst){ r.show(); }
		}
	}

	uint8_t *seq1;
	uint8_t *seq2;
	int64_t bwt_begin;
	int64_t bwt_end;
	uint8_t *str;
	uint8_t *old_str;
	BWT_IDX_loader *idx;
	std::vector<bwt_same_SA> rst;

	uint8_t *pac = NULL;
	//old rst
	bool is_rev;
	int64_t golbal_offset;
	int search_length;

	void init(BWT_IDX_loader *idx_, int search_length_) {
		search_length = search_length_;
		seq1 = (uint8_t*) xcalloc(search_length, 1);
		seq2 = (uint8_t*) xcalloc(search_length, 1);
		str = seq1;
		old_str = seq2;
		bwt_begin = 0;
		bwt_end = idx_->get_seq_len();
		pac = idx_->get_pac();
		idx = idx_;
		rst.clear();

		is_rev = 0;
		golbal_offset = 0;
		next();

	}

	bool next(){
		if(bwt_begin == bwt_end) return false;
		rst.clear();
		if(!idx->is_seq_len(golbal_offset)){
			rst.emplace_back(bwt_begin-1,golbal_offset,is_rev);
		}

		for(;bwt_begin < bwt_end; bwt_begin++){
			std::swap(str, old_str);
			golbal_offset = idx->bwt_sa(bwt_begin);
			//uint64_t golbal_offset = random() % seq_len;
			golbal_offset = idx->sa_depos(golbal_offset, is_rev);
			//show string:
 			int64_t l = 0;
			if(is_rev){
				int64_t beg = golbal_offset - search_length;
				int64_t end = golbal_offset;
				for (int64_t k = end; k > beg; --k)
					str[l++] = 3 - _get_pac(pac, k);
			}else{
				int64_t beg = golbal_offset;
				int64_t end = golbal_offset + search_length;
				for (int64_t k = beg; k < end; ++k)
					str[l++] = _get_pac(pac, k);
			}
			//compare
			bool str_same = true;
			for(int i = search_length - 1; i >= 0; i--){
				if(str[i] != old_str[i]){ str_same = false; break; }
			}
			if(str_same == false){
				bwt_begin++;
				break;//when the string is not same
			}
			else
				rst.emplace_back(bwt_begin,golbal_offset,is_rev); //when the string is same
		}
		return true;
	}
};

class BWT_SA_compactor{
	uint64_t* SA_data;
	uint64_t SA_data_size;
	uint64_t* SA_bucket_idx;
	uint64_t SA_bucket_idx_size;

public:
	void load(char *bwa_idx_dir){
		std::vector<int64_t> cpt_sa_info_v;
		{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa_info"); load_int64_v_from_file(str, cpt_sa_info_v); fprintf(stderr, "compact %s" , str); }
		SA_data_size = cpt_sa_info_v[2];
		SA_bucket_idx_size = cpt_sa_info_v[3];
		fprintf(stderr, "SA data size %ld \n SA index size %ld \n" , SA_data_size, SA_bucket_idx_size);
		SA_data = (uint64_t*) xcalloc(SA_data_size, sizeof(uint64_t));
		SA_bucket_idx = (uint64_t*) xcalloc(SA_bucket_idx_size, sizeof(uint64_t));

		//load basic info
		{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa"); FILE * f = xopen(str, "rb"); xread(SA_data,sizeof(uint64_t), SA_data_size,f); fprintf(stderr, "sa load success \n");}
		{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa_idx"); FILE * f = xopen(str, "rb"); xread(SA_bucket_idx,sizeof(uint64_t), SA_bucket_idx_size,f);  fprintf(stderr, "sa idx load success \n");}

		if(false){
			for(int k = 800; k < 900; k++){
				uint64_t l = SA_bucket_idx[k];
				uint64_t r = SA_bucket_idx[k+1]-1;

				fprintf(stderr, "BG, %ld %ld; !!!\t", l, r);
				for(uint64_t i = l; i <= r + 1; i++){
					fprintf(stderr, "%ld \t", (SA_data[i] & SA_IDX_MASK));
				}
				fprintf(stderr, "\n");
			}
		}
	}

	struct SA_MAP_RST{
		int64_t bwt_idx;
		BASIC_MAP_rst r;
		SA_MAP_RST(int64_t bwt_idx_, int32_t wb_ID_, int32_t POS_, int8_t ori_){
			r.wb_ID = wb_ID_;
			r.aln_position_in_WB = POS_;
			r.is_rev = ori_;
			if(POS_ < 0 || POS_ > 0x3fffff){ fprintf(stderr, "FATAL ERROR 1: wb_ID %d, POS %d\n", r.wb_ID, r.aln_position_in_WB); xassert(0,"");}
			bwt_idx = bwt_idx_;
		}
	};


	//INPUT:int bwa_r_ID, and local_offset
	//output: wb_ID; new local_offset and is_REF
	static int get_wb_ID(ALN_ONLINE::MM_idx_loader *hap_idx, int bwa_r_ID, uint &local_offset, bool &is_REF, bool is_rev){
		int wb_ID = 0;
		if(bwa_r_ID >= (int)hap_idx->wb_info_size){
			bwa_r_ID -= hap_idx->wb_info_size;
			if(false){ fprintf(stderr, "REF: bwa_r_ID %d , local_offset %d\n", bwa_r_ID, local_offset); }
			if(bwa_r_ID >= hap_idx->total_chr_number){ return -1; }//skiping mapping result after chrM
			//get the begin of reads in reference
			uint32_t read_bg;
			if(is_rev){
				if(local_offset > 150)
					read_bg = local_offset - ( 150 -1);
				else{
					read_bg = 0;
				}
			}else{
				read_bg = local_offset;
			}
			wb_ID = hap_idx->convert_ref_coordinate_to_wb_coordinate(bwa_r_ID, read_bg, WB_POS_MAX);
			if(is_rev){//restore local_offset using read_bg
				local_offset = read_bg + (150 -1);
			}else{
				local_offset = read_bg;
			}
			is_REF = true;
		}else{
			wb_ID = bwa_r_ID;
			is_REF = false;
		}
		return wb_ID;
	}

	void remove_WB_result_duplication(std::vector<BWT_MAP_rst_CPT_SA> *rst_v, std::vector<BWT_MAP_rst_CPT_SA> &all_rst){
		rst_v->clear();
		std::sort(all_rst.begin(), all_rst.end(),BWT_MAP_rst_CPT_SA::cmp_by_position);
		if(false){ /*DEBUG: show all ori results:*/ for(auto & A: all_rst){ fprintf(stderr, "After sort\t"); A.print();	}}

		if(all_rst.empty()){/*DO NOTHING*/}
		else if(all_rst.size() == 1){
			rst_v->emplace_back(all_rst[0]);
		}else{
			rst_v->emplace_back(all_rst[0]);//always keep the first one
			uint old_wb = all_rst[0].wb_ID;
			for(uint rID = 1; rID < all_rst.size(); rID++){
				if(old_wb == all_rst[rID].wb_ID){
					if(all_rst[rID].is_REF) //always keep the REF
						rst_v->emplace_back(all_rst[rID]);
				}
				else{
					rst_v->emplace_back(all_rst[rID]); //always keep the first one
					old_wb = all_rst[rID].wb_ID;
				}
			}
		}
		if(false){ /*DEBUG: show all result*/ for(auto & A: *rst_v){ 	fprintf(stderr, "After select\t"); A.print(); } }
	}

	int generate_SA_index_from_BWA_SA(int argc, char *argv[]) {
			const uint64_t BWA_SA_BUFF_SIZE = 1000000;
			uint32_t target_read_length = 150;
			uint32_t window_block_size = 300;

			char *bwa_idx_dir = argv[optind];
			char *wb_dir = argv[optind + 1];
			//idx used for alignment
			BWT_IDX_loader bwt_idx;
			//load BWT
			bwt_idx.restore_idx(bwa_idx_dir, true, true, true, true);

			FILE * cpt_sa = NULL, * cpt_sa_idx = NULL;
			{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa");  cpt_sa = xopen(str, "wb"); }
			{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa_idx"); cpt_sa_idx = xopen(str, "wb"); }

			//{//load wb index
			ALN_ONLINE::MM_idx_loader hap_idx;
			hap_idx.load_all_index(wb_dir, NULL, true);
			//loading hap_string size
			std::vector<std::vector<uint32_t>> hap_end_offset_table;
			{
				//get all length of haplotype
				ALN_ONLINE::hap_string_loader_single_thread hl_r1;
				std::string ref_str;
				ref_str.resize(window_block_size);
				memset(&(ref_str[0]), 'A', window_block_size);
				//count the wb data
				hap_end_offset_table.resize(hap_idx.wb_info_size);
				for (uint64_t window_ID = 0; window_ID < hap_idx.wb_info_size;
						window_ID++) {
					std::vector<std::string> &hap_seq_v = hl_r1.get_string_list(
							window_ID, ref_str, &hap_idx);
					if (!hap_seq_v.empty()) {
						hap_end_offset_table[window_ID].resize(hap_seq_v.size());
						hap_end_offset_table[window_ID][0] = hap_seq_v[0].size();
						for (uint32_t i = 1; i < hap_seq_v.size(); i++)
							hap_end_offset_table[window_ID][i] = hap_seq_v[i].size()
									+ hap_end_offset_table[window_ID][i - 1];
					}
				}
			}

			BWT_same_Checker ck;
			ck.init(&bwt_idx, target_read_length);
			uint bwa_r_ID_old_buff = 0;

			SA_bucket_idx_size = (((bwt_idx.bns->l_pac + 255) >> 8) + 1) << 1;
			SA_bucket_idx = (uint64_t *)xcalloc(SA_bucket_idx_size,sizeof(uint64_t));
			std::vector<uint64_t> SA_store_buff;
			SA_store_buff.resize(BWA_SA_BUFF_SIZE);
			uint32_t cur_buff_size = 0;
			uint64_t block_idx = 0;
			uint32_t write_idx = 0;
			std::vector<BWT_MAP_rst_CPT_SA> all_rst;
			std::vector<BWT_MAP_rst_CPT_SA> rm_dup_rst;

			int ccc = 0;

			while (ck.next()) {
				all_rst.clear();
				rm_dup_rst.clear();
				uint local_offset;
				int bwa_r_ID;
				bool is_REF = false;

				for(auto &r : ck.rst){
					uint64_t golbal_offset = r.golbal_offset;
					bool is_rev = r.is_rev;
					int bwa_r_ID = bwt_idx.bns->bns_pos2rid_buff(golbal_offset, local_offset, bwa_r_ID_old_buff);
					//check if the result is reference
					int wb_ID = get_wb_ID(&hap_idx, bwa_r_ID, local_offset, is_REF, is_rev);
//					if(is_REF){
//						fprintf(stderr, " ");
//					}
					if(wb_ID == -1)	continue;
//					if(r.bwt==328332547){
//						fprintf(stderr, "a");}

					//haplotype check, skip wrong haplotypes
					if(!is_right_haplotype(hap_end_offset_table, bwa_r_ID, local_offset, r.is_rev, is_REF)){
						if(false){ fprintf(stderr, "result is skip\t"); fprintf(stderr, "wb_ID %d, local_offset %d, ori %d is_REF %d\n", wb_ID, local_offset, is_rev, is_REF); }
						continue;
					}
//					//new code:
//					if(is_rev && !is_REF && local_offset < 149){
//						fprintf(stderr, "FATAL ERROR un_compact_sa_rst1, wb_ID %d @ POS %d\n ", wb_ID, local_offset);
//						if(!is_right_haplotype(hap_end_offset_table, bwa_r_ID, local_offset, r.is_rev, is_REF)){
//							if(false){ fprintf(stderr, "result is skip\t"); fprintf(stderr, "wb_ID %d, local_offset %d, ori %d is_REF %d\n", wb_ID, local_offset, is_rev, is_REF); }
//							continue;
//						}
//					}
					all_rst.emplace_back(r.bwt, wb_ID, local_offset, is_rev, is_REF);
				}
				remove_WB_result_duplication(&rm_dup_rst, all_rst);
				//sort
				std::sort(rm_dup_rst.begin(), rm_dup_rst.end(),BWT_MAP_rst_CPT_SA::cmp_by_bwt_idx);

				if(rm_dup_rst.size() == 2 && rm_dup_rst[0].POS < 10000){
					//fprintf(stderr, " ");
					ccc++;
				}

				for (uint32_t i = 0; i < rm_dup_rst.size(); i++) {
					uint64_t bg_bwt = rm_dup_rst[i].bwt_idx;
					uint64_t idx_bucket = bg_bwt >> SA_IDX_OFFSET;			  //256 per block
					uint64_t bwt_left_value = bg_bwt & SA_IDX_MASK;

					uint64_t store_value = rm_dup_rst[i].wb_ID;
					store_value <<= 22;
					store_value += rm_dup_rst[i].POS;
					store_value <<= 1;
					store_value += rm_dup_rst[i].is_rev;
					store_value <<= SA_IDX_OFFSET;
					store_value += bwt_left_value;
					//dump store value
					SA_bucket_idx[idx_bucket]++;
					SA_store_buff[cur_buff_size++] = store_value;

//					if(true){
//						fprintf(stderr, "BG, %ld %ld %ld; !!!\t", idx_bucket, store_value & SA_IDX_MASK, cur_buff_size - 1 );
//						fprintf(stderr, "\n");
//					}

					if (cur_buff_size == BWA_SA_BUFF_SIZE) {
						//
						fwrite(&(SA_store_buff[0]), sizeof(uint64_t), cur_buff_size, cpt_sa);
						fprintf(stderr, "\n write idx %ld, %d @ rst_v[i].bwt_idx %ld\n", write_idx++, cur_buff_size, rm_dup_rst[0].bwt_idx);
						fprintf(stderr, "ccc %d\n", ccc);
						cur_buff_size = 0;
					}
				}
			}

			//store the final part
			fwrite(&(SA_store_buff[0]), sizeof(uint64_t), cur_buff_size, cpt_sa);
			fprintf(stderr, "\n write idx %ld \n", write_idx++);
			//final process
			SA_data_size = 0;
			for (uint i = 0; i < SA_bucket_idx_size - 1; i++) {
				SA_data_size += SA_bucket_idx[i];
				SA_bucket_idx[i] = SA_data_size - SA_bucket_idx[i];
			}
			SA_bucket_idx[SA_bucket_idx_size - 1] = SA_data_size;
			fwrite(SA_bucket_idx, sizeof(uint64_t), SA_bucket_idx_size, cpt_sa_idx);

			//store data in info:
			std::vector<int64_t> cpt_sa_info_v;
			cpt_sa_info_v.emplace_back(window_block_size);
			cpt_sa_info_v.emplace_back(target_read_length);
			cpt_sa_info_v.emplace_back(SA_data_size);
			cpt_sa_info_v.emplace_back(SA_bucket_idx_size);

			{	char str[1024];	strcpy(str, bwa_idx_dir); strcat(str, ".compact_sa_info"); dump_int64_v_to_file(str, cpt_sa_info_v); }

			fclose(cpt_sa_idx);
			fclose(cpt_sa);

			return 0;
		}

	void sa_search(uint64_t bwt_bg, uint64_t bwt_ed, std::vector<uint64_t> &SA_pos_result)
	{
		int64_t max_load = 2000;
		SA_pos_result.clear();
		if(bwt_bg > bwt_ed)
			return;
		int64_t sa_idx_start;
		int64_t m; uint64_t tmp; int64_t r_bg;

		uint64_t idx_bucket_begin = bwt_bg >> SA_IDX_OFFSET;
		uint64_t idx_bucket_end = bwt_ed >> SA_IDX_OFFSET;
		{
			int64_t l = SA_bucket_idx[idx_bucket_begin];
			int64_t r = r_bg = SA_bucket_idx[idx_bucket_begin+1]-1;
			sa_idx_start = r + 1;
			uint64_t key = bwt_bg & SA_IDX_MASK;

			{
				while (l <= r){
					m = (r + l)>> 1;
					tmp = SA_data[m] & SA_IDX_MASK ;
					if(tmp > key){
						sa_idx_start = m;
						r = m - 1;
					}
					else if(tmp < key){
						l = m + 1;
					}else{
						sa_idx_start = m;
						break;
					}
				}
			}
		}
		//
		if(((idx_bucket_begin) == (idx_bucket_end))){//with in same bucket
			uint64_t end_key = bwt_ed & SA_IDX_MASK;
			for(int64_t idx=sa_idx_start;idx <= r_bg; idx ++){
				if(max_load < idx - sa_idx_start){
					SA_pos_result.clear(); break;
				}
				if((SA_data[idx] & SA_IDX_MASK) <= end_key){
					SA_pos_result.emplace_back(SA_data[idx] >> SA_IDX_OFFSET);
				}else{
					break;
				}
			}
		}
		else
		{//NOT same bucket
			int64_t l = SA_bucket_idx[idx_bucket_end];
			int64_t r = SA_bucket_idx[idx_bucket_end+1]-1;

			int64_t sa_idx_end = l-1;
			uint64_t key = bwt_ed & SA_IDX_MASK;

			{
				while (l <= r){
					m = (r + l)>> 1;
					tmp = SA_data[m] & SA_IDX_MASK ;
					if(tmp > key){
						r = m - 1;
					}
					else if(tmp < key){
						sa_idx_end = m;
						l = m + 1;
					}else{
						sa_idx_end = m;
						break;
					}
				}
			}

			if(max_load < sa_idx_end - sa_idx_start){
				SA_pos_result.clear(); return;
			}
			for(int64_t idx=sa_idx_start; idx <=sa_idx_end; idx ++){
				SA_pos_result.emplace_back(SA_data[idx] >> SA_IDX_OFFSET);
			}
		}
	}

	bool is_right_haplotype(std::vector<std::vector<uint32_t>> &hap_end_offset_table,
			int wb_ID, uint local_offset, bool is_rev, bool is_REF){
		if(is_REF) // keep all result in reference
			return true;
		std::vector<uint32_t> &cur_wb_hap_end_offset_list = hap_end_offset_table[wb_ID];
		uint32_t read_bg, read_ed;
		if(is_rev){
			if(local_offset < ( 150 -1)){
				return false;
			}
			read_bg = local_offset - ( 150 -1);
			read_ed = local_offset;
		}else{
			read_bg = local_offset;
			read_ed = local_offset + 150 -1;
		}
		//haplotype check:
		//uint read_ed = local_offset + 150 - 1;//[ ]//616, 600,900
		for(uint i = 0; i < cur_wb_hap_end_offset_list.size(); i++){
			uint hap_bg = (i == 0)?0:(cur_wb_hap_end_offset_list[i-1]);//0, 300, 600, 900~~~~
			uint hap_ed = cur_wb_hap_end_offset_list[i];//300,600,900~~~~
			if(false)
				fprintf(stderr, " hap_bg %d , hap_ed %d, read_bg %d,\t", hap_bg, hap_ed, read_bg);
			if(read_ed > hap_ed)				{ /*do nothing*/}
			else if(read_bg < hap_bg)
				return false;
			else if(read_bg == hap_bg)
				return false;//remove the read at the just beginning of WB
			else
				return true;
		}
		return false;
	}
};

#endif /* BWT_CPT_SA_HPP_ */


