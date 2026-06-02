/*
 * variant_count.hpp
 *
 *  Created on: Jun 8, 2022
 *      Author: zaq
 */

#ifndef VARIANT_CALLER_MAIN_HPP_
#define VARIANT_CALLER_MAIN_HPP_
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include "../CPPLIB/tools.hpp"
#include "../BWT_idx/var_map.hpp"
#include "../index_building/haplotype_online.hpp"
#include "known_var_hap_counter.hpp"
extern "C" {
#include "../clib/utils.h"
}

namespace Aln_online{

//show the support is from which haplotype??
struct Var_supp_haplotype_item{
	uint32_t wb_ID;
	int var_ID; // if var ID is -1, it support REFerence
	int hap_ID; //which hap is it from
	int supp;	//depth in the haplotype
	Var_supp_haplotype_item(uint32_t wb_ID_, int var_ID_, int hap_ID_, int supp_){
		wb_ID = wb_ID_;
		var_ID = var_ID_,
		hap_ID = hap_ID_;
		supp = supp_;
	}
};

class Position_depth_item{
private:
	struct var_supp_item{
		uint32_t wb_ID;
		uint32_t var_ID;
		uint32_t depth;
		//uint8_t /*BOOL*/var_from_pre_wb;
		var_supp_item(uint32_t wb_ID_, uint32_t var_ID_, uint32_t depth_){
			wb_ID = wb_ID_;
			var_ID = var_ID_;
			depth = depth_;
		}

		void print(FILE * log_f){
			fprintf(log_f,
					",wb_ID %d var_ID %d, depth %d \n",
					wb_ID, var_ID, depth);//, var_from_pre_wb);
		}
	};
	int32_t ref_supp;

public:
	std::vector<var_supp_item> variant_supp;
	//select the max haplotype set for each high support vars
	void get_max_depth_hap(uint min_var_depth_select, std::set<int> &select_hap_idx_set, FILE * log_f){
		if(false)
			fprintf(log_f, "POS: ref_supp %d, already_combined %d\n" , ref_supp, already_combined);
		for(auto & v: variant_supp){
			if(v.depth >= min_var_depth_select){
				//select max_depth hap as the main HAP
				for(auto & vh:vh_l){
					if(vh.var_ID == (int)v.var_ID){
						 select_hap_idx_set.emplace(vh.hap_ID);
					}
				}
			}
			if(false) v.print(log_f);
		}
	}

	std::vector<Var_supp_haplotype_item> vh_l;
	bool already_combined;

	Position_depth_item(){ ref_supp = 0; vh_l.clear(); already_combined = false;}

	void update_ref_supp(uint32_t window_id, int hap_ID_, int supp_){
		ref_supp += supp_;
		vh_l.emplace_back(window_id, -1, hap_ID_, supp_);
	}

	void update_ref_supp_novel(int supp_){
		ref_supp += supp_;
	}

	void update_variant_supp(uint32_t window_id, uint32_t variant_ID, int hap_ID_, int supp_) {
		bool find_flag = false;
		for (auto &c_var : variant_supp)  // Loop over all the maps
			if(c_var.var_ID == variant_ID){
				c_var.depth += supp_;
				find_flag = true;
				break;
			}
		if (!find_flag) {
			variant_supp.emplace_back(window_id, variant_ID, supp_);
			//variant_supp.back().var_from_pre_wb = false;
			//variant_supp.back().var_ID_pre_wb = MAX_uint32_t;
		}
		vh_l.emplace_back(window_id, variant_ID, hap_ID_, supp_);
	}

	void combine(Position_depth_item & to_combine, std::vector<Window_t> &window_info){
		//combine REF depth
		ref_supp += to_combine.ref_supp;
		//combine var list: support by this and both
		for(auto & var: variant_supp){
			for(auto & to_combine_var: to_combine.variant_supp){
				if(window_info[var.wb_ID].var_list[var.var_ID].isSAME_ref_alt(window_info[to_combine_var.wb_ID].var_list[to_combine_var.var_ID])){
					var.depth += to_combine_var.depth;
					to_combine_var.depth = var.depth;
					break;
				}
			}
		}
		//combine var list: support only by to_combine
		for(auto & to_combine_var: to_combine.variant_supp){
			bool already_combined = false;
			for(auto & var: variant_supp){
				if(window_info[var.wb_ID].var_list[var.var_ID].isSAME_ref_alt(window_info[to_combine_var.wb_ID].var_list[to_combine_var.var_ID])){
					already_combined = true;
					break;
				}
			}
			if(!already_combined)
				variant_supp.emplace_back(to_combine_var);
		}
	}

	void store_in_vcf(FILE* out, std::vector<Window_t> &window_info, std::vector<VCF_item> &VCF_BUFFs,int MAX_LOW_depth);
};

class Window_block_counter_Known_Var{
	Hap_counter BLANK_window_process_block;

	//pos depth
	std::map<uint32_t, Position_depth_item> BLANK_pos_depth;

	//buffs
	std::map<std::tuple<uint32_t, uint8_t>, std::tuple<uint32_t, uint32_t>> variant_count;
	std::map<uint32_t, int> candidate_pos;

	std::map<uint32_t, Position_depth_item> pos_depth[4];
	std::map<uint32_t, Position_depth_item> pre_wb_pos_depth;
	Hap_counter pre_window_process_block;
	Hap_counter window_process_block[4];

public:
	ALN_ONLINE::Simple_ref_handler *ref_h;

	FILE* log_f;

	std::map<uint32_t, Position_depth_item> *get_pos_depth_by_id(int id){
		if(id == -1) return &(pre_wb_pos_depth);
		else return &(pos_depth[id]);
	}

	Hap_counter *get_process_block_by_id(int id){
		if(id == -1) return &(pre_window_process_block);
		else return &(window_process_block[id]);
	}

	void init(std::vector<Window_t> &window_info, FILE* log_f_);
	void destroy(){
		//todo::
	}

	void set_new_wb(std::vector<Window_t> &window_info, uint32_t window_basic){
		for(int i = 0; i < 4; i++)	window_process_block[i].init( window_info[window_basic + i].total_hap_length, window_basic + i, log_f);
	}
	void add_signal(int32_t window_count, int32_t pos, int read_len){
		window_process_block[window_count].add_signal(pos, read_len);
	}
	void generate_known_var_candidate(std::vector<Window_t> &window_info);

	void final_process_after_var_calling(){
		//final:
		std::swap(pre_wb_pos_depth, pos_depth[3]);
		std::swap(pre_window_process_block, window_process_block[3]);
	}

	void show_statistics(){

	}

private:
	void position_depth_counting(uint32_t l, std::vector<Window_t> &window_info);
	void haplotype_counter_split_and_realignment(Hap_counter &hap_count, std::vector<Window_t> &window_info);
} ;

struct VC_PARA{
	//basic option
	int MAX_WB_NUM;
	int MIN_WB_NUM;
	char * window_block_info_fn;
	char * mapping_rst_fn;
	char * index_path;

	int get_option(int argc, char *argv[]){
		FILE * log_f = stderr;
		options_list l;
		l.show_command(log_f, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  variant_calling  [Options] <INDEX_PATH> <WB_index.map> <Map_Resultst.bin>\n");
	    l.add_title_string("  Basic:   \n");
	    l.add_title_string("    <WB_index.map>  FILE   	The map file\n");
	    l.add_title_string("    <Map_Resultst.bin>  FILE    mapping result file of 'bwt_aln' \n");
		//thread number
		//l.add_option("MAX_MAPPING_RST_NUM",  'R', "DEBUG: max load number of mapping results", true, MAX_int32t); l.set_arg_pointer_back((void *)&MAX_MAPPING_RST_NUM);
	    l.add_option("MIN_WB_NUM", 			'M', "DEBUG: the first wb in used", true, 0); l.set_arg_pointer_back((void *)&MIN_WB_NUM);
	    l.add_option("MAX_WB_NUM", 			'W', "DEBUG: the last  wb in used", true, MAX_int32t); l.set_arg_pointer_back((void *)&MAX_WB_NUM);

		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(log_f);
			return 1;
		}
		l.show_c_value(log_f);
		if (argc - optind < 3)
			return l.output_usage();
		index_path = strdup(argv[optind]);
		window_block_info_fn = strdup(argv[optind + 1]);
		mapping_rst_fn = strdup(argv[optind + 2]);
		return 0;
	}
};

struct MAP_rst_48{
	uint8_t r[6];
	void set(int32_t chr_ID, int32_t POS){
		uint64_t store = ((uint64_t)chr_ID << 22) + (POS & 0x3fffff);
		r[5] = store >> 40;
		r[4] = store >> 32;
		r[3] = store >> 24;
		r[2] = store >> 16;
		r[1] = store >> 8;
		r[0] = store;
	}
	void restore(int32_t &chr_ID, int32_t &POS){
		uint64_t store = 0;

		store += (r[5]); store <<= 8;
		store += (r[4]); store <<= 8;
		store += (r[3]); store <<= 8;
		store += (r[2]); store <<= 8;
		store += (r[1]); store <<= 8;
		store += (r[0]);
		chr_ID = store >> 22;
		POS = store & 0x3fffff;
	}
};

#define WINDOW_HAPS_BLOCK_NUM 1000000 //1M
class Variant_Caller_main {
private:
	VC_PARA para;
	std::vector<Window_t> window_info;
	MAP_rst_48 * mapping_rst_load_buffer;
	FILE *mapping_rst_fp;

	void show_simple_vcf_header(FILE* output){
		fprintf(output,
				"##fileformat=VCFv4.2\n"
				"##contig=<ID=chr1,length=248956422>\n"
				"##contig=<ID=chr2,length=242193529>\n"
				"##contig=<ID=chr3,length=198295559>\n"
				"##contig=<ID=chr4,length=190214555>\n"
				"##contig=<ID=chr5,length=181538259>\n"
				"##contig=<ID=chr6,length=170805979>\n"
				"##contig=<ID=chr7,length=159345973>\n"
				"##contig=<ID=chr8,length=145138636>\n"
				"##contig=<ID=chr9,length=138394717>\n"
				"##contig=<ID=chr10,length=133797422>\n"
				"##contig=<ID=chr11,length=135086622>\n"
				"##contig=<ID=chr12,length=133275309>\n"
				"##contig=<ID=chr13,length=114364328>\n"
				"##contig=<ID=chr14,length=107043718>\n"
				"##contig=<ID=chr15,length=101991189>\n"
				"##contig=<ID=chr16,length=90338345>\n"
				"##contig=<ID=chr17,length=83257441>\n"
				"##contig=<ID=chr18,length=80373285>\n"
				"##contig=<ID=chr19,length=58617616>\n"
				"##contig=<ID=chr20,length=64444167>\n"
				"##contig=<ID=chr21,length=46709983>\n"
				"##contig=<ID=chr22,length=50818468>\n"
				"##contig=<ID=chrX,length=156040895>\n"
				"##contig=<ID=chrY,length=57227415>\n"
				"##contig=<ID=chrM,length=16569>\n"
				"##FILTER=<ID=LOW_COV,Description=\"Low coverage\">\n"
				"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant TYPE\">\n"
				"##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total read depth summed across all datasets\">\n"
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus Genotype across all datasets with called genotype\">\n"
				"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INTEGRATION\n"
				"");
	}

public:
	int run(int argc, char *argv[]) {
		FILE * out_f = stdout;
		FILE * log_f = stderr;

		para.get_option(argc, argv);
		Window_t::load_variant_map(para.window_block_info_fn, para.MIN_WB_NUM, para.MAX_WB_NUM,window_info);
		//load bin reference
	    ALN_ONLINE::Simple_ref_handler ref_h;
		ref_h.load_bin_ref(para.index_path);

		mapping_rst_fp = xopen(para.mapping_rst_fn, "rb+");
		fseek(mapping_rst_fp, 0, SEEK_END);// non-portable
		uint64_t total_read_numer = ftell(mapping_rst_fp)/6;
		uint64_t ave_read_depth = total_read_numer*150/2800000000;
		fprintf(log_f, "tatal_read_number %ld depth %ld\n", total_read_numer, ave_read_depth);
		rewind(mapping_rst_fp);

		mapping_rst_load_buffer = (MAP_rst_48*) calloc(WINDOW_HAPS_BLOCK_NUM, sizeof(MAP_rst_48));
		//output header
		show_simple_vcf_header(out_f);
		uint32_t buff_idx = 0;
		int window_basic = -10;
		int32_t window_id = -1, pos;
		Window_block_counter_Known_Var *wb_c = (Window_block_counter_Known_Var*) calloc(1, sizeof(Window_block_counter_Known_Var));
		wb_c->init(window_info, log_f);

		wb_c->ref_h = &ref_h;
		while (true) {
			//re-load data
			uint64_t true_load_rst_num = fread(mapping_rst_load_buffer, sizeof(MAP_rst_48), WINDOW_HAPS_BLOCK_NUM, mapping_rst_fp);
			if(true_load_rst_num == 0){ break;} //final
			if(true && window_id > para.MAX_WB_NUM) break;
			buff_idx = 0;
			while (buff_idx < true_load_rst_num) {
				mapping_rst_load_buffer[buff_idx++].restore(window_id, pos);
				if(window_id < para.MIN_WB_NUM){
					continue;
				}
				//if(rand() % 100 > 60) continue; //debug code: 30X
				//fprintf(log_f, "window_id %d, pos %d\n ", window_id, pos);
				if(true && window_id > para.MAX_WB_NUM) break;
				//new a window block counter
				if (window_id > window_basic + 3) {
					//get the read number
					wb_c->generate_known_var_candidate(window_info);
					//init 4 window block per time, not only one, because wb ID is not sorted
					window_basic = ((window_id >> 2) << 2);
					//if(window_basic == 6460)fprintf(log_f, " ");
					wb_c->set_new_wb(window_info, window_basic);
				}
				if(pos >= (int)window_info[window_id].total_hap_length)
					pos = 0x3fffff - pos + window_info[window_id].total_hap_length - 300;

				int32_t st_pos = pos;
				int32_t ed_pos = st_pos + 150;
				//window_info[window_id].with_in_repeat_check(st_pos, ed_pos);
				if(false && window_id == 6430 && st_pos > 1051 && st_pos < 1351){
					fprintf(log_f, "st_pos %d, ed_pos %d\n", st_pos, ed_pos);
				}
				//if(ed_pos > st_pos + 16){//skip too short reads
				wb_c->add_signal(window_id - window_basic, pos, ed_pos - st_pos);
				//}
//				if(wb_c->window_process_block[2].hap_count_list[1254] != 0)
//					fprintf(log_f, " ");
			}
		}
		//final process
		if(window_basic )
		wb_c->generate_known_var_candidate(window_info);
		fclose(mapping_rst_fp);
		// parse variant
		return 0;
	}
};
}


#endif /* VARIANT_CALLER_MAIN_HPP_ */
