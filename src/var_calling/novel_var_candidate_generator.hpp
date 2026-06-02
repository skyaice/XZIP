/*
 * novel_var_caller.hpp
 *
 *  Created on: 2023年2月12日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_NOVEL_VAR_CALLER_HPP_
#define SRC_BWT_ONLINE_NOVEL_VAR_CALLER_HPP_

#include "../aln_RST_DEF.hpp"
#include "../index_building/haplotype_online.hpp"
#include "../CPPLIB/Assembler/assembler.hpp"
#include "../BWT_idx/var_map.hpp"
#include "GT_TYPE.hpp"
#include <map>
#include <string>
#include <algorithm>

extern "C"
{
#include "../clib/kswlib/kalloc.h"
#include "../clib/kswlib/ksw2.h"
}

namespace Aln_online{

#define Position_depth_item_Novel_SNP 0
#define Position_depth_item_Novel_INS 1
#define Position_depth_item_Novel_DEL 2

#define Position_depth_item_Novel_CHANGE_DETAIL 0
#define Position_depth_item_Novel_ASSEMBLY 1

struct novel_var_supp_item{
	uint32_t wb_ID;
	uint8_t var_type;
	std::vector<uint8_t> alt_bin;//blank in the beginning

    std::vector<char> ref_char;
    std::vector<char> alt_char;

	int var_len;
	int ref_supp;
	int var_supp;
	int signal_TYPE;//describe how the signal is generated, from assembly or "change detail" ?
	//FOR snp
	novel_var_supp_item(){}
	void set(uint32_t wb_ID, uint8_t var_type, int var_length, int ref_supp, int var_supp, int signal_TYPE){
		this->wb_ID = wb_ID;
		this->var_type = var_type;
		this->var_len = var_length;
		this->ref_supp = ref_supp;
		this->var_supp = var_supp;
		this->signal_TYPE = signal_TYPE;
	}
};
//it is to store the mapping results
//the result in the same reference position is stored in one "Position_depth_item_Novel"
class Position_depth_item_Novel{

	uint32_t ref_supp_number_sum;
	uint32_t var_pos_in_wb;
	uint8_t ref_char_bin;

	public:
	std::vector<novel_var_supp_item> variant_supp;

	bool is_signal_change_detail(){
		if(variant_supp.empty())
			return false;
		return (variant_supp[0].signal_TYPE == Position_depth_item_Novel_CHANGE_DETAIL);
	}

	void init(uint32_t var_pos_in_wb, uint8_t ref_char_bin){
		this->var_pos_in_wb = var_pos_in_wb;
		ref_supp_number_sum = 0;
		this->ref_char_bin = ref_char_bin;
	}

	void add_ref_supp(int ref_supp_number){
		ref_supp_number_sum += ref_supp_number;
	}

	void add_candidate_by_assembly(uint32_t wb_ID, uint32_t change_type, uint8_t *alt, int var_length, int var_supp){
		//check duplication:
		for(novel_var_supp_item & check_v:variant_supp){
			if(check_v.var_type == change_type && check_v.var_len == var_length){
				bool is_same = true;
				if(change_type == Position_depth_item_Novel_SNP){
					is_same = (check_v.alt_bin[0] == alt[0]);
				}else if(change_type == Position_depth_item_Novel_INS){
					for(int i = 0; i < var_length; i++)
						if(check_v.alt_bin[i+1] != alt[i])
							is_same = false;
				}
				if(is_same){
					check_v.var_supp = MAX(check_v.var_supp, var_supp);
					return;
				}
			}
		}

		variant_supp.emplace_back();
		variant_supp.back().set(wb_ID, change_type, var_length, 0, var_supp, Position_depth_item_Novel_ASSEMBLY);
		variant_supp.back().alt_bin.clear();
		if(change_type == Position_depth_item_Novel_SNP){
			variant_supp.back().alt_bin.emplace_back(alt[0]);
		}else if(change_type == Position_depth_item_Novel_INS){
			variant_supp.back().alt_bin.emplace_back(ref_char_bin);
			for(int i = 0; i < var_length; i++){
				variant_supp.back().alt_bin.emplace_back(alt[i]);
			}
		}
	}

	void add_candidate_by_change_detail(uint32_t wb_ID, uint32_t change_detail, int var_supp){
		uint32_t change_type = (change_detail >> 30);
		//the var length value is ONLY meaningful for INDELs;
		//for SNPs, it is always set 0
		uint32_t var_length = ((change_detail >> 26) & 0xf);//t
		//parse the variants
		//for CPX stored result
		//      000:0 INS: A
		//      001:1 INS: C
		//      100:4 DEL
		//      101:5 SNP: + 1
		//      102:6 SNP: + 2
		//      103:7 SNP: + 3
		uint8_t alt_char_bin = ref_char_bin;
		if(change_type == 0){//SNP
			alt_char_bin += ((change_detail & 0x7) - 4);
			alt_char_bin %= 4;
			variant_supp.emplace_back();
			variant_supp.back().set(wb_ID, change_type, 0, 0, var_supp, Position_depth_item_Novel_CHANGE_DETAIL);
			variant_supp.back().alt_bin.clear();
			variant_supp.back().alt_bin.emplace_back(alt_char_bin);
		}else if(change_type == 1){//INS
			variant_supp.emplace_back();
			variant_supp.back().set(wb_ID, change_type, var_length, 0, var_supp, Position_depth_item_Novel_CHANGE_DETAIL);
			variant_supp.back().alt_bin.clear();
			variant_supp.back().alt_bin.emplace_back(ref_char_bin);
			for(uint i = 0; i < var_length; i++){
				uint8_t ins_char_bin = ((change_detail >> (24 - i*2)) & 0x3);
				variant_supp.back().alt_bin.emplace_back(ins_char_bin);
			}
			//alt_char_bin = change_detail;
		}else{//DEL
			variant_supp.emplace_back();
			variant_supp.back().set(wb_ID, change_type, var_length, 0, var_supp, Position_depth_item_Novel_CHANGE_DETAIL);
			variant_supp.back().alt_bin.clear();
			variant_supp.back().alt_bin.emplace_back(alt_char_bin);
		}
	}

	void store_in_vcf(FILE* out,
			std::vector<Window_t> &window_info,
			std::vector<VCF_item> &VCF_BUFFs,
			int MAX_LOW_depth, ALN_ONLINE::Simple_ref_handler *ref_h){
	  uint32_t total_depth = ref_supp_number_sum;
	  for(novel_var_supp_item & var_info : variant_supp)
	  {
	    if(var_info.var_supp == 1)         continue;
	    //simple GT
	    int GT_TYPE;// 0: 0/0; 1: 0/1; 2: 1/1
	    int PASS_TYPE; // 0: NOT_PASS;1: PASS; 2:LOW_COV
	    int32_t true_ref_supp = total_depth - var_info.var_supp;
		if(true_ref_supp < 0) true_ref_supp = 0;

	    get_GT_type_NOVEL(var_info.var_supp, true_ref_supp, GT_TYPE, PASS_TYPE, MAX_LOW_depth);
	    //uint32_t var_ID_this_wb, var_ID_pre_wb;
	    uint32_t var_pos = var_pos_in_wb + window_info[var_info.wb_ID].st_pos;
	    uint32_t chr_ID = window_info[var_info.wb_ID].chr_ID;

	    std::vector<char> &ref = var_info.ref_char; ref.clear();
	    std::vector<char> &alt = var_info.alt_char; alt.clear();

		if(var_info.var_type == 0){//SNP
			ref.emplace_back("ACGT"[ref_char_bin]);
			alt.emplace_back("ACGT"[var_info.alt_bin[0]]);
			var_info.var_len = 0;
		}else if(var_info.var_type == 1){//INS
			ref.emplace_back("ACGT"[ref_char_bin]);
			for(uint8_t c: var_info.alt_bin)
				alt.emplace_back("ACGT"[c]);
		}else{//DEL
			//var_pos -= 1;
			std::string ref_string;
			ref_h->load_ref_from_buff(chr_ID, var_pos, var_info.var_len + 1, ref_string);
			for(int i = 0; i < var_info.var_len + 1; i++){
				ref.emplace_back(ref_string[i]);
			}
			var_info.var_len = -var_info.var_len;
			alt.emplace_back(ref_string[0]);
		}
		ref.emplace_back(0);
		alt.emplace_back(0);
		VCF_BUFFs.emplace_back();
	    VCF_BUFFs.back().set(
	      chr_ID, var_pos,
		  &(ref[0]), &(alt[0]),
	      30, PASS_TYPE, GT_TYPE,
		  var_info.var_len,
	      true_ref_supp ,var_info.var_supp,
	      -1,
	      var_info.wb_ID,
		  var_info.signal_TYPE + 1
	    );
	  }
	}

	//	struct VCF_item{
	//		uint32_t chr_ID;
	//		uint32_t var_pos;
	//		const char *ref_p;
	//		const char *alt_p;
	//		int QUAL;
	//		int PASS_TYPE; // 0: NOT_PASS;1: PASS; 2:LOW_COV
	//		int GT_TYPE;// 0: 0/0; 1: 0/1; 2: 1/1
	//		int var_len;
	//		int ref_supp;
	//		int var_supp;
	//
	//		uint32_t var_ID;
	//		uint32_t var_wb;
	//	};


};

//this structure is used to store the read information before read assembly, those data is useful for explaining the assembly results
struct ASS_reads_info{
//	ASS_reads_info(int read_in_ref_offset_, int read_st_pos_, int read_ed_pos_, int read_list_index_, int soft_left_, int soft_right_, int number_NM_){
//		read_in_ref_offset = read_in_ref_offset_;
//		read_list_index = read_list_index_;
//		//read_st_pos = read_st_pos_;
//		//read_ed_pos = read_ed_pos_;
//		//soft_left = soft_left_;
//		//soft_right = soft_right_;
//		//number_NM = number_NM_;
//	}

	ASS_reads_info(int read_in_ref_offset_, int read_list_index_){
		read_in_ref_offset = read_in_ref_offset_;
		read_list_index = read_list_index_;
	}

	int read_list_index;
	int read_in_ref_offset;
	//int read_st_pos;//read_in_ref_offset - left_clip
	//int read_ed_pos;//read_in_ref_offset + middle_size
	//int soft_left;
	//int soft_right;
	//int number_NM;
	void print(FILE* log_file, int first_read_ID){
		fprintf(log_file, " id: %d(%d) ", read_list_index - first_read_ID, read_list_index);
		fprintf(log_file, " ref pos %d\t", read_in_ref_offset);
	}
};

struct Ass_Block{
	//input read info
	std::vector<ASS_reads_info> ass_read_list;
	//input read data
	std::vector<std::string> reads;
	//output: read result
	std::vector<AssemblyContig> contigs;

	void clear(){
		ass_read_list.clear();
		reads.clear();
	}

	void add_read_basic_info(int read_in_ref_offset_, int read_list_index_){
		ass_read_list.emplace_back(read_in_ref_offset_, read_list_index_);//true set data
	}

	void add_read_string_bin(uint8_t * read_str, int read_len, bool store_in_reverse){
		reads.emplace_back();
		std::string &s_store = reads.back(); //string to store: string pointer type
		s_store.resize(read_len);
		for(int i = 0; i < read_len; i++)
			s_store[i] = ("ACGTN"[read_str[i]]);//store binary to acgt
		if(store_in_reverse){
			char_seq_reverse(read_len, &(s_store[0]), true);
		}
	}

	void add_read_string(char * read_str, int read_len, bool store_in_reverse){
		reads.emplace_back(); std::string &s_store = reads.back(); //string to store: string pointer type
		s_store.resize(read_len);
		for(int i = 0; i < read_len; i++)
			s_store[i] = read_str[i];//store binary to ACGT
		if(store_in_reverse){
			char_seq_reverse(read_len, &(s_store[0]), true);
		}
	}

	void run_assembly(MainAssemblyHandler *am){
		std::swap(am->reads, reads);
		am->assembley();
		std::swap(am->contigs, contigs);
		std::swap(am->reads, reads);//swap back the read list
	}

};

struct Contig_Depth_Counter{
	//Suggest
	struct SUGGEST_POS_LIST_ITEM{
	  SUGGEST_POS_LIST_ITEM(int suggest_pos_){
	    suggest_pos = suggest_pos_;
	    read_count = 0;
	    low_wrong_base_read_number = 0;//
	    high_wrong_base_read_number = 0;
	  }

	  void add_read_start_pos(int wrong_base){
	    if(wrong_base < 2)
	      low_wrong_base_read_number ++;
	    else if(wrong_base < 200)
	      high_wrong_base_read_number ++;
	    else return;
	    read_count++;
	  }

	  void printf(FILE * log){
	    fprintf(log, "[SG: %d, read count: %d]\t", suggest_pos, read_count);
	  }

	  static inline int cmp_by_read_count(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.read_count > b.read_count; }
	  //static inline int cmp_by_read_position(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.ave_read_start < b.ave_read_start; }

	  int suggest_pos;
	  int read_count;
	  int low_wrong_base_read_number;
	  int high_wrong_base_read_number;

	};
	std::set<int> remove_read_set;
	std::map<int, int > suggest_st_pos_map;
	std::vector<uint16_t> contig_depth;
	std::vector<SUGGEST_POS_LIST_ITEM> SUGGEST_pos_list;
	int min_suggested_alignment_pos;
	int final_read_suggest_alignment_start_pos;

#define MAX_WRONG_BASE 8
	void get_contig_depth_and_suggented_contig_positions(AssemblyContig & contig, int ori_read_number,
			std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f ){
	const char *contig_seq = contig.seq.c_str();
	  int contig_seq_len = contig.seq.size();
	  //step1: get all the read positions at assembled CONTIG
	  //set read position for each actions
	  remove_read_set.clear();
	  for (auto &ca : contig.actions)
		if(ca.read_ID < ori_read_number)
		  ca.set_read_pos(read_list[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.wordLength, remove_read_set, ass_read_list[ca.read_ID].read_in_ref_offset);

	  //step2: find the suggestion alignment position of CONTIG at reference from each read position
	  suggest_st_pos_map.clear();

	  //get contig coverage
	  if((int)contig_depth.size() < contig_seq_len) contig_depth.resize(contig_seq_len);
	  memset(&(contig_depth[0]), 0, contig_seq_len*sizeof(uint16_t));

	  for (AssemblyReadAction &ca : contig.actions){
		ca.wrong_base = 999;
		if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
		//if(ass_read_list[ca.read_ID].signal_type != Read_type::SR) continue;
		if(remove_read_set.find(ca.read_ID) != remove_read_set.end()) continue;
		//set coverage:
		const char* read_seq = read_list[ca.read_ID].c_str();
		int st_pos_ref = ca.position_in_contig - contig.ass_begin_offset_in_contig - ca.position_read;
		int ed_pos_ref = st_pos_ref + read_list[ca.read_ID].size();
		int st_pos_read = 0; if(st_pos_ref < 0){  st_pos_read -= st_pos_ref;st_pos_ref = 0;}
		ed_pos_ref = MIN(contig_seq_len, ed_pos_ref);
		ca.wrong_base = 0;
		for(int i = st_pos_ref; i < ed_pos_ref && ca.wrong_base <= MAX_WRONG_BASE; i++, st_pos_read++){
		  if(read_seq[st_pos_read] != 'N' && contig_seq[i] != read_seq[st_pos_read])
			ca.wrong_base++;
		}
		if(ca.wrong_base <= MAX_WRONG_BASE){
		  st_pos_read = 0;
		  for(int i = st_pos_ref; i < ed_pos_ref; i++, st_pos_read++){
			if(contig_seq[i] == read_seq[st_pos_read]) contig_depth[i] ++;
		  }
		  std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
		  if(it != suggest_st_pos_map.end())  it->second++;
		  else suggest_st_pos_map[ca.suggest_contig_offset_in_ref] = 1;
		}else{
		  std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
		  if(it != suggest_st_pos_map.end())  it->second--;
		}
	  }
	  //if no result is there, simple selected the first read position
	  if(suggest_st_pos_map.empty())  suggest_st_pos_map[contig.actions[0].suggest_contig_offset_in_ref] = 1;

	  //step3: find the max suggest position
	  //get max suggest alignment position
	  //try simple merge:
	  int max_suggention = 0; int max_count = 0;
	  for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++)
		if(max_count < it->second) {max_suggention = it->first; max_count = it->second;}

	  //step4: find the suggestion list
	  //if a suggestion covers most of read, [unique_alignment_pos] is true and use it as alignment suggestion, otherwise store all suggestions in a list [suggent_pos_list]
	  // when no result or result > 1: [unique_alignment_pos] = false
	  //suggest_coverage only use when suggent_pos_list.size() > 1;
	  SUGGEST_pos_list.clear();

	  for(std::map<int, int >::value_type &sug : suggest_st_pos_map){
		if(sug.second >= 2){
		  SUGGEST_pos_list.emplace_back(sug.first);
		  SUGGEST_POS_LIST_ITEM & sug = SUGGEST_pos_list.back();
		  int c_suggest_pos = sug.suggest_pos;
		  for(auto &ca : contig.actions){
			if(ca.read_ID < ori_read_number && ca.suggest_contig_offset_in_ref == c_suggest_pos)
			  sug.add_read_start_pos(ca.wrong_base);
		  }
		  if(sug.low_wrong_base_read_number < 1 || sug.high_wrong_base_read_number > sug.low_wrong_base_read_number * 4)//only 0 support reads or less than 20% right supported reads
			SUGGEST_pos_list.erase(SUGGEST_pos_list.end() - 1);
		}
	  }
	  min_suggested_alignment_pos = MAX_int32t;
	  for(SUGGEST_POS_LIST_ITEM & s: SUGGEST_pos_list){
		  min_suggested_alignment_pos = MIN(s.suggest_pos, min_suggested_alignment_pos);
	  }
	  //print suggestion alignment range information
	  if(false){
		//all possible ranges
		fprintf(log_f, "All suggestions:\t");
		for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++) fprintf(log_f, "[SUG: %d CN: %d]\t", it->first, it->second);
		fprintf(log_f, "\n Used suggestions: \t");
		for(SUGGEST_POS_LIST_ITEM & sug : SUGGEST_pos_list) sug.printf(log_f);
		//possible true ranges
		if(SUGGEST_pos_list.empty())      {fprintf(log_f, "\nNO suggestion\n"); }
		else if(SUGGEST_pos_list.size() == 1)   {fprintf(log_f, "\nUNIQUE: [MAX_SUG: %d]\n", max_suggention);}
		else                  {fprintf(log_f, "\nMULTY\n"); }

		if(true){
		  for (auto &ca : contig.actions)
			if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
			{
			  ass_read_list[ca.read_ID].print(log_f, 0);
			  ca.print(log_f);
			}
		  fprintf(log_f, "\n");
		}
	  }
	}

};

void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
extern uint8_t charToDna5n[];
struct Contig_String_aligner{
private:
	int min_var_position;
	int max_var_position;
public:
	static void test(){
		Contig_String_aligner ca;
		ca.init();
		//char ref[500] =   "AAAGAATTATTAACTATAGTTAATAGAAACCAAGTGGAAATAATGTGAGCCAAAAGGAGCATTTTATCATAAGAATGTAAGGGTATGTCCTGAAACACAAAGCAAGGAACATAATTTGGTATGTTCCAGAACATACCAAGTAGGTTCCAGAACAAGTAGGTTCTGGAACCAGAAACTTGTCTCTAACTCCATCTTTCATCTCTATTTCTCCGTGAA";
		//char query[500] = "AAAGAATTATTAACTATAGTTAATAGAAACCAAGTGGAAATAATGTGAGCCAAAAGGAGCATTTTATCATAAGAATGTAAGGGTATGTCCTGAAACACAAAGCAAGGAACATAATTTGG   GTT        T                  CAAGTAGGTTCTGGAACCAGAAACTTGTCTCTAACTCCATCTTTCATCTCTATTTCTCCGTGAA";
		char ref[500] =     "TAATTCTTCTTTCACAGACAGAAGATGGATGTACAAAAAGTCAAGGTGGAGGAGACCTCTATCTTCGTAGGATAGAGGATAGTCCTTGAATGGGCGAGCATTGAATATTTGTGGATGGCTTAGATAAGTGAGTGTATTCTTCC";
		char query[500] =   "TAATTCTTCTTTCACAGACAGAAGATGGATGTACAAAAAGTCAAGGTGGAGGAGACCTCTATCTTCGTAGGAAGTACAAGGACAAACTGAATAGTCCTTGAATGGGCGAGCATTGAATATTTGTGGATGGCTTAGATAAGTGAGTGTATTCTTCC";

		std::vector<uint8_t> ref_bin;
		std::vector<uint8_t> query_bin;
		//store bin contig
		int ref_seq_len = strlen(ref);
		ref_bin.resize(ref_seq_len);
		for (int i = 0; i < ref_seq_len; ++i)
			ref_bin[i] = charToDna5n[(uint8_t)ref[i]];

		int q_seq_len = strlen(query);
		uint16_t * contig_depth = (uint16_t *)xcalloc(q_seq_len, 2);
		query_bin.resize(q_seq_len);
		for (int i = 0; i < q_seq_len; ++i)
			query_bin[i] = charToDna5n[(uint8_t)query[i]];

		ca.setRef(&(ref_bin[0]), ref_seq_len, 0, 0);
		ca.align_non_splice(&(query_bin[0]), q_seq_len, 0, ref_seq_len - 1);

		ca.printf_alignment_detail(stderr, 0, contig_depth);

	}

	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option_INDEL();
		ksw_gen_mat_D();
		min_var_position = 150;
		max_var_position = 300;
	}

	void setRef(uint8_t * ref_string, int ref_len, 	int chr_ID_, int global_ref_pos_){ tseq = ref_string; tlen = ref_len; chr_ID = chr_ID_;global_ref_pos = global_ref_pos_;}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}

	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		//ref_end_pos += 50;
		last_ref_st_pos = ref_st_pos;
		last_ref_end_pos = ref_end_pos;

		if(false){
			//debug code:
			uint32_t debug_qlen = qlen;
			uint8_t* debug_qseq = qseq;
			for(uint32_t i = 0; i < debug_qlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_qseq[i]]);
			fprintf(stderr, "\n");
			uint32_t debug_tlen = ref_end_pos - ref_st_pos;
			uint8_t* debug_tseq =  tseq + ref_st_pos;
			for(uint32_t i = 0; i < debug_tlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_tseq[i]]);
			fprintf(stderr, "\n");
		}
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	int adjustCIGAR(){
		uint32_t n_cigar = ez.n_cigar;
		int adj_size = cigar_adjust(&n_cigar, ez.cigar, false, 15);
		ez.n_cigar = n_cigar;
		return adj_size;
	}

	void printf_alignment_detail(FILE * output, int suggest_st_pos, uint16_t * contig_depth){//, int ref_pos_enough_match_base){
		log_output = output;
		printCIGAR();
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0) return;
		if(true) printQuerySeq(suggest_st_pos);//print contig sequence
		if(true) printTargetSeq();
		int contig_coverage = 0;
		if(true){ contig_coverage = print_X_E_sequence(suggest_st_pos);}//, ref_pos_enough_match_base); } //print X/= sequence
		if(true){ print_coverage(suggest_st_pos, contig_coverage, contig_depth); } 		//print coverage sequence
		print_SV_canditate(suggest_st_pos);
 		fprintf(output, "\n");
	}

	void new_candidate(int var_pos, uint8_t ref_char_bin, std::map<uint16_t, Position_depth_item_Novel> &candidate_var_positions){
			//get candidate
				auto it = candidate_var_positions.find(var_pos);
				if(it == candidate_var_positions.end()){
					Position_depth_item_Novel new_depth_item;
					new_depth_item.init(var_pos, ref_char_bin);
					candidate_var_positions[var_pos] = new_depth_item;
				}
	}

	bool get_canditate_Vars_and_ref_depth(int ref_st_pos, std::map<uint16_t, Position_depth_item_Novel> &candidate_var_positions, uint16_t * contig_depth, int wb_id,
			uint32_t * signal_read_depth_buff){
		bool next_wb_need_assembly_again = false;
		uint32_t* bam_cigar = ez.cigar;
		int n_cigar = ez.n_cigar;
		int target_index = ref_st_pos;
		int query_index = 0;
		int max_INDEL_length = 100; //todo:: used as PARA:
		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:
				for(int i = 0; i < c_length; i++){
					if(tseq[target_index + i] != qseq[query_index + i]){
						int snp_position = target_index + i;
						int snp_contig_position = query_index + i;
						int snp_depth = contig_depth[snp_contig_position];
						if(snp_position >= min_var_position && snp_position < max_var_position && snp_depth > 1){//SNP
							new_candidate(snp_position, tseq[target_index + i], candidate_var_positions);
							candidate_var_positions[snp_position].add_candidate_by_assembly(wb_id, Position_depth_item_Novel_SNP, qseq + query_index + i, 0, snp_depth);
						}else if(snp_position >= max_var_position && snp_depth > 1){
							next_wb_need_assembly_again = true;
						}
					}
				}
				//set depth:
				for(int i = 0; i < c_length; i++)
					signal_read_depth_buff[target_index + i] = MAX(signal_read_depth_buff[target_index + i], contig_depth[query_index + i]);
				//add index:
				target_index += c_length; query_index += c_length;	break;//M
			case 1:	//insertion
				//store the data
				{
					if(c_length <= max_INDEL_length){
						int ins_position = target_index - 1;
						int ins_contig_position = query_index - 1;
						int ins_depth = contig_depth[ins_contig_position];
						if(cigar_ID != 0 && cigar_ID != n_cigar - 1 && ins_position >= min_var_position && ins_position < max_var_position && ins_depth > 1){
							new_candidate(ins_position, tseq[ins_position], candidate_var_positions);
							candidate_var_positions[ins_position].add_candidate_by_assembly(wb_id, Position_depth_item_Novel_INS, qseq + query_index, c_length, ins_depth);
						}else if(ins_position >= max_var_position && ins_depth > 1){
							next_wb_need_assembly_again = true;
						}
					}
				}
				query_index += c_length;
				break;//I, int chr_ID; int
			case 2:	 //Deletions
				//store the data
				{
					if(c_length <= max_INDEL_length){
						int del_position = target_index - 1;
						int del_contig_position = query_index - 1;
						int del_depth = contig_depth[del_contig_position];
						if(cigar_ID != 0 && cigar_ID != n_cigar - 1 && del_position >= min_var_position && del_position < max_var_position && del_depth > 1){
							new_candidate(del_position, tseq[del_position], candidate_var_positions);
							candidate_var_positions[del_position].add_candidate_by_assembly(wb_id, Position_depth_item_Novel_DEL, NULL, c_length, del_depth);
						}else if(del_position >= max_var_position && del_depth > 1){
							next_wb_need_assembly_again = true;
						}
					}
				}
				//store the depth data
				for(int i = 0; i < c_length; i++)
					signal_read_depth_buff[target_index + i] = MAX(signal_read_depth_buff[target_index + i], contig_depth[query_index + i]);
				target_index += c_length;
				break;
			case 3:	target_index += c_length; query_index += c_length; break;//N, print N
			case 4:	target_index += c_length; break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}
		return next_wb_need_assembly_again;
	}

	inline uint8_t getTseq(int i){ return tseq[i];}

	uint8_t* get_ref(){return tseq;}
	int get_tlen(){return tlen;}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

	int get_cigar(uint32_t** bam_cigar){
		*bam_cigar = ez.cigar;
		return ez.n_cigar;
	}

	//part7: candidate SVs
	int last_ref_st_pos;  //the reference start and end position of last alignment
	int last_ref_end_pos;

private:
	//part1 : basic informations
	int chr_ID;
	int global_ref_pos;
	//part2: reference and target
	uint8_t *tseq; int tlen;
	uint8_t *qseq; uint32_t qlen; //query temp: used for analysis
	//part3: alignment result
	ksw_extz_t ez;
	//part4: buffs
	void *km;
	//part5: options
	int8_t mata_D[25]; int8_t match_D; int8_t mismatch_D;
	int8_t gap_open_D; int8_t gap_ex_D; int8_t gap_open2_D; int8_t gap_ex2_D;
	uint16_t zdrop_D; int bandwith; //for DNA zdrop = 400, 200 for RNA
	int flag;
	//part6: logs
	FILE * log_output;

	void copy_option_SV(){
			match_D = 2;
			mismatch_D= 10;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 0;
			zdrop_D= gap_open2_D + 600; //for DNA zdrop = 400, 200 for RNA
			bandwith = 600;//zdrop_D;
			flag = 0;
	}

	void copy_option_INDEL(){
			match_D = 3;
			mismatch_D= 5;
			gap_open_D= 4;
			gap_ex_D = 2;
			gap_open2_D= 100;
			gap_ex2_D = 1;
			zdrop_D= gap_open2_D + 500; //for DNA zdrop = 400, 200 for RNA
			bandwith = 2000;//zdrop_D;
			flag = 0;
	}

	void ksw_gen_mat_D(){
		int8_t l,k,m;
		for (l = k = 0; l < 4; ++l) {
			for (m = 0; m < 4; ++m) { mata_D[k] = l == m ? match_D : -(mismatch_D);	/* weight_match : -weight_mismatch */ k++; }
			mata_D[k] = 0; // ambiguous base
			k++;
		}
		for (m = 0; m < 5; ++m) { mata_D[k] = 0; k++; }
	}

	void printCIGAR(){
		fprintf(log_output, "\n");
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0){
			for(int i = 0; i < (int)qlen; i++) {fprintf(log_output, "%c", "ACGT"[qseq[i]]);}
			fprintf(log_output, "\n ???NO alignment result\n");
		}else{
			fprintf(log_output, "CIGAR: ");
			uint32_t* bam_cigar = ez.cigar;
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				fprintf(log_output, "%d", cigar_len);
				char cigar_char = '?';
				switch(type){
				case 0:	cigar_char = 'M'; break;
				case 1:	cigar_char = 'I'; break;
				case 2:	cigar_char = 'D'; break;
				case 3:	cigar_char = 'N'; break;
				case 4:	cigar_char = 'S'; break;
				}
				fprintf(log_output, "%c", cigar_char);
			}
			fprintf(log_output, "\n");
		}
	}

	void printQuerySeq(int suggest_st_pos){
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		int output_index = 0;
		int seq_i = 0;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", "ACGT"[qseq[seq_i]]); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
	}

	void printTargetSeq(){
		for(int seq_i = 0; seq_i < tlen; seq_i++){
			fprintf(log_output, "%c", "ACGT"[tseq[seq_i]]);
		}
		fprintf(log_output, "\n");
	}

	int print_X_E_sequence(int suggest_st_pos){//, int ref_pos_enough_match_base){
		int contig_coverage = 0;
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++)
				{
					if(qseq[seq_i] == tseq[output_index]){
						contig_coverage++;
						//if(output_index < ref_pos_enough_match_base) fprintf(log_output, "M");//1 == M; 2 == X; 0 == -;
						//else
						fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
					}
					else{
						fprintf(log_output, "X");//1 == M; 2 == X; 0 == -;
					}
					output_index++;
				}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
		return contig_coverage;
	}

	void print_SV_canditate(int suggest_st_pos){
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		bool have_long_SV = false;
		int minIndelLen = 1; //for Indels
		//int minSVlen = 50;//for SVs

		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	break;//M
			case 1:	case 2: case 3:	case 4:
				if(cigar_len >= minIndelLen) have_long_SV = true;
				break;//I, print nothing
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}

		int ref_index = suggest_st_pos;
		if(have_long_SV){
			fprintf(log_output, "Variant detected:");
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:	ref_index += cigar_len; break;//M
				case 1:	if(cigar_len >= minIndelLen) {fprintf(log_output, "Insertion: "); fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } seq_i += cigar_len; break;//I, int chr_ID; int
				case 2:	if(cigar_len >= minIndelLen) {fprintf(log_output, "Deletion: ");  fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } ref_index += cigar_len; break;//D, print '-'
				case 3:	ref_index += cigar_len; break;//N, print N
				case 4:	ref_index += cigar_len; break;//S, print -
				default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
		}
		fprintf(log_output, "\n");
	}

	void print_coverage(int suggest_st_pos, int contig_coverage, uint16_t * contig_depth){
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", contig_depth[seq_i] + '#'); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "0"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\t");
		fprintf(log_output, "\nCigar sequence: ");
		for(int i = 0; i < cigar_len; i++){
			fprintf(log_output, "%d%c", (bam_cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(bam_cigar[i] & BAM_CIGAR_MASK)]);
		}
		fprintf(log_output, "contig_coverage: [%d bp]\t", contig_coverage);

		fprintf(log_output, "\t");
		fprintf(log_output, "\n");
	}
};

struct Single_WB_signal_handler_Novel{
	//basic
	uint32_t wb_id;
	uint32_t wb_basic;
	//stored data
	std::vector<Single_read_aln_rst_unpack> gap_g;//GAP DATA
	std::vector<Single_read_aln_rst_unpack> gap_um;//UNMAPPED DATA
	bool already_assembly;
	//change detail map for each WB
	std::vector<uint32_t> parse_change_detail_buff;
	std::map<uint32_t, uint32_t> wb_change_detail;
	//read mapping position list
	bool is_mapping_position_sorted;
	std::vector<uint16_t> read_mapping_position;

	bool without_any_signal(){ return (gap_g.empty() && gap_um.empty()); }
	void clear();
	void sort_mapping_position();
	void combine_mapping_position_with_next(Single_WB_signal_handler_Novel * next,int step_length){
		int adjust_size = (step_length << 1);
		for(uint16_t p: next->read_mapping_position)
			read_mapping_position.emplace_back(adjust_size + p);
	}
	void change_detail_parse(std::vector<uint16_t> &ori, std::vector<uint32_t> &rst, int r_read_len);
	void add_signal(Single_read_aln_rst_unpack & r);
};

#define MIN_CHANGE_DETAIL_DEPTH_NOVEL_VAR 3
//candidate
class Candidate_result_Novel{
public:
	//parameters
	int max_depth_in_region_bg;
	int max_depth_in_region_ed;
	int wb_step_length;
	//WB info
	ALN_ONLINE::MM_idx_loader *wb_idx_p;
	ALN_ONLINE::window_block_info * pre_wb_info;
	//ref
	ALN_ONLINE::Simple_ref_handler *ref_h;//reference index
	std::vector<uint8_t> ref_string_buff;//450 bp

	//statistic basic value
	//min read depth for candidate var:
	uint32_t min_read_number_UM_ASS;
	uint32_t min_read_number_GAP_UM_ASS;
	//the input:
	Single_WB_signal_handler_Novel *pre_wb_signals;
	Single_WB_signal_handler_Novel *next_wb_signals;

	bool is_assembly;
	bool with_contig;

	//the output
	//P1: gap read depth
	uint32_t max_depth;
	uint32_t* ref_depth_buff;
	int signal_read_depth_buff_size;
	//P2: combined change details
	std::map<uint32_t, uint32_t> combined_change_detail;
	//P3: assemblies
	MainAssemblyHandler am;
	Ass_Block *ass_b;
	//P3.1 contig aligner
	Contig_String_aligner *c_a;
	Contig_Depth_Counter d_c;
	std::vector<uint8_t> bin_contig;
	//P4: the middle results of variants calling: it store the candidate positions and candidate vars
	std::map<uint16_t, Position_depth_item_Novel> candidate_var_positions;

public:
	//REF + ALT
	uint32_t get_total_depth_at_position(uint32_t novel_pos){
		return ref_depth_buff[novel_pos];
	}
	//ONLY REF
	uint32_t get_ref_depth_at_position(uint32_t novel_pos){
		if(ref_depth_buff[novel_pos] == 0 || is_assembly){
			return ref_depth_buff[novel_pos];
		}
		//BY ALN, but the depth is not 0
		int ref_depth = 0;
		int wb_step_length = 150;
		for(int i = 0; i < 2; i++){
			std::vector<Single_read_aln_rst_unpack> * l = (i == 0)?(&pre_wb_signals->gap_g):(&next_wb_signals->gap_g);
			int pos_adjust = (i == 0)?0:wb_step_length;
			for(Single_read_aln_rst_unpack & r: *l){
				uint32_t read_st = r.read_hamming_mapping_position + pos_adjust;
				uint32_t read_ed = read_st + r.read_len;
				if(novel_pos < read_st || novel_pos >= read_ed)
					continue;
				int indel_adjust = 0;
				bool with_var = false;
				for(uint16_t c:r.change_detail){
					uint32_t var_pos = read_st + (c >> 3) + indel_adjust;
					uint16_t change_type = c & 0x7;
					if(change_type <= 4){//INDEL
						if(var_pos - 1 == novel_pos)
							with_var = true;
					}else //SNP
						if(var_pos == novel_pos)
							with_var = true;
					if(change_type == 4) 	 indel_adjust ++;//DEL
					else if(change_type < 4) indel_adjust --;//INS
				}
				if(with_var == false){
					ref_depth++;
				}
			}
		}
		return ref_depth;
	}

	uint32_t get_max_gap_read_depth(){
		return max_depth;
	}

	std::map<uint16_t, Position_depth_item_Novel> &get_candidate_var_positions(){
		return candidate_var_positions;
	}

	void init(ALN_ONLINE::MM_idx_loader *wb_idx_p, ALN_ONLINE::Simple_ref_handler *ref_h){
		signal_read_depth_buff_size = 450; //todo::
		max_depth_in_region_bg = 150;
		max_depth_in_region_ed = 299;
		wb_step_length = 150;
		ref_depth_buff = (uint32_t *)xmalloc((signal_read_depth_buff_size + wb_step_length) * sizeof(uint32_t));
		min_read_number_UM_ASS = 6;//20% of read depth
		min_read_number_GAP_UM_ASS = 8;//30% of read depth
		this->wb_idx_p = wb_idx_p;
		this->ref_h = ref_h;
		c_a = (Contig_String_aligner *)xcalloc(1, sizeof(Contig_String_aligner));
		c_a->init();//todo:: set special parameters
		//c_a.setZdrop(50, 500);
		ass_b = (Ass_Block *)xcalloc(1, sizeof(Ass_Block));
	}

	void clear(){

	}

private:
	//get the max depth of gap read in the overlap region of two nearby WB
//	S1：添加read的信号，每个read添加两个信号，起始与结束分别添加一次，之后按照位置排序。
//	S2：同时考虑两个相邻的wb，合计长度为 450 bp；其中的 150bp的深度是可以直接计算的，周围的两个150bp不能直接计算。
//	S3：当前深度值记录为0。遍历拍过序的 read信号列表，遇见“起始”信号，当前深度值加1；遇见结束信号，深度值减一。
//	S4：统计中间的150bp的最高深度记录。
//return the position list
	void get_ref_depth_in_wb_using_gap_alignmnet(uint32_t * signal_read_depth_buff){
		pre_wb_signals->combine_mapping_position_with_next(next_wb_signals, wb_step_length);
		pre_wb_signals->sort_mapping_position();
		//next_wb_signals->sort_mapping_position();
		//get the depth:
		uint32_t cur_depth = 0;
		uint32_t position_bg = 0;// 0~449
		uint32_t position_ed = 0;// 0~449

		for(uint16_t  r: pre_wb_signals->read_mapping_position){
			position_ed = (r>>1);
			position_ed = MIN(position_ed, 450);
			//store depth list
			for(uint32_t i = position_bg; i < position_ed; i++){
				signal_read_depth_buff[i] = cur_depth;
			}

			//reset depth
			if(r & 0x1)	cur_depth--; //is end
			else		cur_depth++;  //is begin
			xassert(cur_depth >= 0, "");
			//reset the begin position
			position_bg = position_ed;
		}

		max_depth = 0;
		for(int i = max_depth_in_region_bg; i < max_depth_in_region_ed; i++)
			max_depth = MAX(max_depth, signal_read_depth_buff[i]);

		//debug code:
		if(false && max_depth != 0){
			fprintf(stderr, "pre_wb_signals->wb_id %d @wb_basic %d , next_wb_signals->wb_id %d @wb_basic %d \n", pre_wb_signals->wb_id, pre_wb_signals->wb_basic,
					next_wb_signals->wb_id, next_wb_signals->wb_basic);
			for(int i = 0; i < 450; i++){
				if(i % 150 == 0)
					fprintf(stderr, "\n");
				fprintf(stderr, "%d\t", signal_read_depth_buff[i]);
			}
			fprintf(stderr, "\n");
		}
	}

	// 函数主要涉及相邻 两个wb的坐标转换问题。
	//10 bit: pos in WB; 3 lower bit: change detail
	void combine_change_detail_list_of_nearby_wb(){
		combined_change_detail.clear();
		//copy the w_follow
		combined_change_detail = pre_wb_signals->wb_change_detail;
		uint32_t next_wb_adjust_number =  (150 << 3);//change:
		//combine the front
		for(auto it = next_wb_signals->wb_change_detail.begin(); it != next_wb_signals->wb_change_detail.end(); it++){
			uint32_t change_detailit = it->first + next_wb_adjust_number;
			auto find_it = combined_change_detail.find(change_detailit);
			if(find_it != combined_change_detail.end()){
				find_it->second += it->second;
			}else{
				combined_change_detail[change_detailit] = it->second;
			}
		}
	}

	//return true: need ASS;
	//return false: just call this way
	bool get_signals_using_change_detail(){
		 int position_bg = wb_step_length;
		 int position_ed = position_bg + wb_step_length;
		 bool with_INDEL_signal = false;
		for(auto combined_change_detail_it = combined_change_detail.begin(); combined_change_detail_it != combined_change_detail.end(); combined_change_detail_it++)
		{
			int var_depth = combined_change_detail_it->second;
			//condition 1:
			if(var_depth < MIN_CHANGE_DETAIL_DEPTH_NOVEL_VAR)
				continue;
			//condition 2:
			int var_pos = ((combined_change_detail_it->first >> 3) & 0x1fff);
			uint32_t change_type = (combined_change_detail_it->first >> 30);
			if(change_type != 0){ //INS or DEL:
				var_pos -= 1;
			}
			if(var_pos < position_bg || var_pos >= position_ed)
				continue;
			//condition 3:
			if(change_type != 0){ //INS or DEL:
				with_INDEL_signal = true;
			}
			auto it = candidate_var_positions.find(var_pos);
			if(it == candidate_var_positions.end()){
				Position_depth_item_Novel new_depth_item;
				uint8_t ref_char_bin = ref_string_buff[var_pos];
				new_depth_item.init(var_pos, ref_char_bin);
				candidate_var_positions[var_pos] = new_depth_item;
			}
			candidate_var_positions[var_pos].add_candidate_by_change_detail(
					pre_wb_signals->wb_id, (combined_change_detail_it->first), var_depth);
		}
		if(combined_change_detail.size() > 1 && with_INDEL_signal)
			return true;
		return false;
	}

public:
	void running_ASS(FILE * log_f, bool show_assembly){
		is_assembly = true;
		//clear:
		ass_b->clear();//clear the ass block
		with_contig = false;
		candidate_var_positions.clear();
		//clear depth counter:
		memset(ref_depth_buff, 0, sizeof(uint32_t) * signal_read_depth_buff_size);
		//read number check: skip the WB when the read number is too High
		uint64_t total_UM_read_number = pre_wb_signals->gap_g.size() + pre_wb_signals->gap_um.size() + next_wb_signals->gap_g.size() + next_wb_signals->gap_um.size();
		if(total_UM_read_number > 300){
			fprintf(log_f, "Skip assembling of this WB because the read number is too high. [ %d ]\n", total_UM_read_number);
			return;
		}
		//load reads signals
		//NOTE: all reads is stored forward, that is: with the same direction of reference
		for(int i = 0; i < 2; i++){
			Single_WB_signal_handler_Novel * h = (i == 0)? pre_wb_signals: next_wb_signals;
			int position_adjust = (i == 0)? 0: 150;
			for(int j = 0; j < 2; j++){
				std::vector<Single_read_aln_rst_unpack> & l = (j == 0)? (h->gap_g): (h->gap_um);
				for(Single_read_aln_rst_unpack & r: l){//load read in the pre_wb_signals
					r.generate_read_string(wb_idx_p, r.read_len);
					ass_b->add_read_basic_info(r.read_hamming_mapping_position + position_adjust, ass_b->ass_read_list.size());
					ass_b->add_read_string_bin(&(r.read_string[0]), r.read_string.size(), false);
					if(show_assembly)//debug
						r.print_detail(log_f);
				}
			}
		}
		uint32_t ori_read_number = ass_b->ass_read_list.size();
		//fprintf(log_f, "Using ASS\n");
		ass_b->run_assembly(&am);
		//get contigs: ONLY the first contig is used
		if(ass_b->contigs.empty() == false){
			with_contig = true;
			AssemblyContig contig = ass_b->contigs[0];
			//Skip contigs in repeat region
			if(contig.ending_reason[0] == 0 && contig.ending_reason[1] == 0){
				//get depth of contig using assembly actions
				d_c.get_contig_depth_and_suggented_contig_positions(contig, ori_read_number, ass_b->reads, ass_b->ass_read_list, stderr);
				//debug:
				if(show_assembly)
					contig.debug_print(log_f);
				//get_ref_string(ref_string_buff);
				c_a->setRef(&(ref_string_buff[0]), ref_string_buff.size(), pre_wb_info->chrID, pre_wb_info->region_st);
				store_bin_contig(contig.seq, bin_contig);
				int read_suggest_alignment_start_pos = d_c.min_suggested_alignment_pos;//the suggested alignment position of contig in reference
				//given more 30 bp of ref for better alignment
				read_suggest_alignment_start_pos -= 30; read_suggest_alignment_start_pos = MAX(read_suggest_alignment_start_pos, 0);
				int read_suggest_alignment_end_pos = ref_string_buff.size();
				c_a->align_non_splice(&(bin_contig[0]), bin_contig.size(), read_suggest_alignment_start_pos, read_suggest_alignment_end_pos);
				//c_a.printf_alignment_detail(stderr, read_suggest_alignment_start_pos, NULL);
				read_suggest_alignment_start_pos += c_a->adjustCIGAR();
				//show the alignment results
				if(show_assembly)
					c_a->printf_alignment_detail(log_f, read_suggest_alignment_start_pos, &(d_c.contig_depth[0]));
				//get_canditate_Vars_and_ref_depth
				bool next_wb_need_assembly_again = c_a->get_canditate_Vars_and_ref_depth(read_suggest_alignment_start_pos, candidate_var_positions, &(d_c.contig_depth[0]), pre_wb_signals->wb_id, ref_depth_buff);
				//note:: when next_wb_need_assembly_again, set already_assembly to false
				d_c.final_read_suggest_alignment_start_pos = read_suggest_alignment_start_pos;
				if(next_wb_need_assembly_again)
					next_wb_signals->already_assembly = false;
			}
		}
	}

	void running_ALN(FILE * log_f, bool show_assembly, bool & running_assembly){
		//get change detail in this block
		is_assembly = false;
		with_contig = false;
		//note: the change detail information may be used in next WB, therefore the value need calculate even
		candidate_var_positions.clear();
		//clear depth counter:
		memset(ref_depth_buff, 0, sizeof(uint32_t) * signal_read_depth_buff_size);
		combine_change_detail_list_of_nearby_wb();
		get_ref_depth_in_wb_using_gap_alignmnet(ref_depth_buff);
		running_assembly = get_signals_using_change_detail();
		if(running_assembly){
			uint total_ass_read_number = pre_wb_signals->gap_um.size() + next_wb_signals->gap_um.size() +  pre_wb_signals->gap_g.size() + next_wb_signals->gap_g.size();
			if(total_ass_read_number < min_read_number_GAP_UM_ASS){
				running_assembly = false;
			}
		}
		if(false && !running_assembly){// debug code
			for(int i = 0; i < 2; i++){
				Single_WB_signal_handler_Novel * h = (i == 0)? pre_wb_signals: next_wb_signals;
				for(int j = 0; j < 2; j++){
					std::vector<Single_read_aln_rst_unpack> & l = (j == 0)? (h->gap_g): (h->gap_um);
					for(Single_read_aln_rst_unpack & r: l){//load read in the pre_wb_signals
						r.generate_read_string(wb_idx_p, r.read_len);
						if(true)//debug
							r.print_detail(log_f);
					}
				}
			}
		}
	}

	//return whether using Assembly
	bool call_novel_vars_single(Single_WB_signal_handler_Novel * pre, Single_WB_signal_handler_Novel *next, FILE * log_f, bool show_assembly){
		//clear candidate
		candidate_var_positions.clear();
		//skip when without any signals
		if(pre->without_any_signal() && next->without_any_signal()){
			//clear depth counter:
			memset(ref_depth_buff, 0, sizeof(uint32_t) * signal_read_depth_buff_size);
			return false;
		}
		pre_wb_signals = pre;
		next_wb_signals = next;

		//S1: get depth string of mapped gap reads
		//S2: basic:
		//set WB info
		pre_wb_info = wb_idx_p->wb_info + pre->wb_id;
		//load reference
		ref_h->load_ref_bin_from_buff(
				pre_wb_info->chrID, pre_wb_info->region_st,
				pre_wb_info->region_length + wb_step_length + wb_step_length,//600 bp:
				ref_string_buff);

		//S2.1 statistics T2: by depth of change detail
		if(false && pre->wb_basic == 546808){
			fprintf(log_f, " ");
		}
		if(show_assembly){
			fprintf(log_f, "pre->wb_id %d pre GAP %ld pre UM %ld\t",pre->wb_id, pre->gap_g.size(), pre->gap_um.size());
			fprintf(log_f, "next->wb_id %d next GAP %ld next UM %ld\n",next->wb_id, next->gap_g.size(), next->gap_um.size());
		}
		//S3: store ALT candidate positions and ALT signals
		//S3.1 statistics T1: by depth of gap and total number of UM reads
		//bool max_depth_is_over_cutoff = (max_depth >= min_read_depth_for_candidate_var);
		//for conditions:
		//C1 is mean the NEXT has enough read to assembly
		//and C2 is means the pre or next have data (at least 2) but not enough to assembly
		//meanwhile "next+pre" has enough reads to assembly
		//pre 0 + next 5: fail ;pre 2 + next 5: pass ;pre 0 + next 10: pass ; pre 10 next 0: fail
		bool max_UM_number_is_over_cutoff_C1 = (pre_wb_signals->gap_um.size() >= (uint64_t)min_read_number_UM_ASS && pre_wb_signals->already_assembly == false);
		bool max_UM_number_is_over_cutoff_C2 = (next_wb_signals->gap_um.size() >= (uint64_t)min_read_number_UM_ASS && next_wb_signals->already_assembly == false);
		bool max_UM_number_is_over_cutoff_C3 =
				(pre_wb_signals->gap_um.size() + next_wb_signals->gap_um.size() >= (uint64_t)min_read_number_UM_ASS) &&
				(pre_wb_signals->gap_um.size() > 2) &&
				(next_wb_signals->gap_um.size() > 2);
		bool running_assembly = max_UM_number_is_over_cutoff_C1 || max_UM_number_is_over_cutoff_C2 || max_UM_number_is_over_cutoff_C3;
		//either call var directly or using assembly; NOT both.
		is_assembly = false;
		with_contig = false;
		if(!running_assembly){//using change detail
			//get change detail in this block
			//note: the change detail information may be used in next WB, therefore the value need calculate even
			running_ALN(log_f, show_assembly, running_assembly);
		}
		if(running_assembly){
			next_wb_signals->already_assembly = true;
			running_ASS(log_f, show_assembly);
		}
		if(show_assembly)
			fprintf(log_f, "END: pre->wb_id %d \n",pre->wb_id);
		return running_assembly;
	}

	//Suggest
	struct SUGGEST_POS_LIST_ITEM{
		SUGGEST_POS_LIST_ITEM(int suggest_pos_){
			suggest_pos = suggest_pos_;
			read_count = 0;
			low_wrong_base_read_number = 0;//
			high_wrong_base_read_number = 0;
			read_start_pos_sum = 0;
			ave_read_start = 0;
		}

		void add_read_start_pos(int read_start_pos, int wrong_base){
			if(wrong_base < 2)
				low_wrong_base_read_number ++;
			else if(wrong_base < 200)
				high_wrong_base_read_number ++;
			else return;

			read_count++;
			read_start_pos_sum += read_start_pos;
		}

		void count_ave_read_start(){ ave_read_start = (float)read_start_pos_sum/((float)(read_count) + 0.1); }

		void printf(FILE * log){
			if(ave_read_start < -5000 || ave_read_start > 5000)
				fprintf(log, "Fatal ERROR, ave_read_start wrong");
			fprintf(log, "[SG: %d, read count: %d, ave_read_start: %f]\t", suggest_pos, read_count, ave_read_start);
		}

		static inline int cmp_by_read_count(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.read_count > b.read_count; }
		static inline int cmp_by_read_position(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.ave_read_start < b.ave_read_start; }

		int suggest_pos;
		int read_count;
		int low_wrong_base_read_number;
		int high_wrong_base_read_number;
		int read_start_pos_sum;
		float ave_read_start;

	};

};

//handler all signals for Novel variants in nearby 4 WBs
class Nearby_WB_Signal_handler_Novel{
	//statistics
	//stored signal data
	std::vector<Single_WB_signal_handler_Novel> cur_4;//signals of current 4 WBs
	std::vector<Single_WB_signal_handler_Novel> old_4;//signals of previous 4 WBs, some informations may be useful
	FILE* log_f;
	bool show_assembly;

	//statistics
	//
	uint32_t total_wb_number;
	uint32_t total_ASS_number;

public:
	//results
	//4 results: that is: (1) WB -1 and 0; (2) WB 0 and 1; (3) WB 1 and 2; (4) WB 2 and 3;
	std::vector<Candidate_result_Novel> result_4;

	void init(ALN_ONLINE::MM_idx_loader *wb_idx, ALN_ONLINE::Simple_ref_handler *ref_h, FILE*vcf_out_log, bool show_assembly);
	void destroy(){
		//todo::
	}
	void clear(uint32_t wb_basic);
	inline void add_signal(Single_read_aln_rst_unpack & r, int window_basic_bg){
		cur_4[r.wb_ID - window_basic_bg].add_signal(r);
	}
	void store_cur_to_old(){ std::swap(cur_4, old_4); }
	void generate_novel_vars_candidate_in_4_nearby_WB();//calling variants in for nearby WB
	void show_statistics(){
		fprintf(log_f, "Nearby_WB_Signal_handler_Novel STATISTICS: total_wb_number %d total_ASS_number %d\n",total_wb_number, total_ASS_number);
	}
};

}


#endif /* SRC_BWT_ONLINE_NOVEL_VAR_CALLER_HPP_ */
