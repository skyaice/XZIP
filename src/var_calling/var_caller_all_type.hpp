/*
 * var_caller_all_type.hpp
 *
 *  Created on: 2023年2月13日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_VAR_CALLING_VAR_CALLER_ALL_TYPE_HPP_
#define SRC_BWT_ONLINE_VAR_CALLING_VAR_CALLER_ALL_TYPE_HPP_

#include "known_var_candidate_generator.hpp"
#include "novel_var_candidate_generator.hpp"

namespace Aln_online{

class Var_generater_all_type{

	FILE* out_f; FILE* log_f;

	ALN_ONLINE::Simple_ref_handler *ref_h;

	int ave_read_depth;
	int MIN_depth_of_low_ALT_check;
	int MIN_depth_of_low_REF_check;
	int MAX_depth_of_low_check;//when the support reads of vars is below 25% of average depth (12 for 50X; 7 for 30X), the reads support the VARs will be realigned
	int MAX_depth_of_low_check_LV2;//when the support reads of vars is below 25% of average depth (12 for 50X; 7 for 30X), the reads support the VARs will be realigned
	int MAX_LOW_depth;

	std::vector<VCF_item> KNWON_VCF_BUFFs;
	std::vector<VCF_item> NOVEL_VCF_BUFFs;

	int known_var_ID;
	int novel_var_ID;
	char *var_name_buff;
	char *sample_name_string;

	struct Compare_main_hap{
		int id;
		std::string str;
		Compare_main_hap(int id_, std::string &str_){
			id = id_;
			std::swap(str, str_);
		}
	};
public:
	void init(int ave_read_depth, FILE* out_f, FILE* log_f, ALN_ONLINE::Simple_ref_handler *ref_h, char *sample_name_string){
		this->ave_read_depth = ave_read_depth;
		this->MIN_depth_of_low_ALT_check = 2;
		this->MIN_depth_of_low_REF_check = 2;
		this->MAX_depth_of_low_check = (float)ave_read_depth*0.3;//when the support reads of vars is below 25% of average depth (15 for 50X; 10 for 30X), the reads support the VARs will be realigned
		this->MAX_depth_of_low_check_LV2 = (float)ave_read_depth*0.4;
		this->MAX_LOW_depth = (float)ave_read_depth*0.08;
		this->MAX_LOW_depth = MAX(3, this->MAX_LOW_depth);
		this->out_f = out_f;
		this->log_f = log_f;
		this->ref_h = ref_h;

		known_var_ID = 0;
		novel_var_ID = 0;
		var_name_buff = (char *)xmalloc(1024*sizeof(char));
		this->sample_name_string = sample_name_string;
	}

	void destroy(){
		//todo::
	}

	//this function only combine the same VAR in novel list to the known VAR list
	//when the novel variant is UNIQ/TRUE, is is kept in original list
	void combine_novel_to_known(Position_depth_item &known, Position_depth_item_Novel & novel, uint novel_position, std::vector<Window_t> &window_info){
		for(int novel_idx = 0; novel_idx < (int)novel.variant_supp.size(); novel_idx++){
			bool is_combined = false;
			for(auto & known_var_info: known.variant_supp){
				bool is_same_var_novel_to_known = false;
				{
					Variant & known_var = window_info[known_var_info.wb_ID].var_list[known_var_info.var_ID];
					novel_var_supp_item & novel_var = novel.variant_supp[novel_idx];
					int novel_var_len = (novel_var.var_type == Position_depth_item_Novel_DEL)?(-novel_var.var_len):(novel_var.var_len);
					if(known_var.get_length() != novel_var_len){
						// do nothing
					}else if(novel_var_len == 0){//SNP
						if(known_var.get_alt()[0] == "ACGT"[novel_var.alt_bin[0]])
							is_same_var_novel_to_known = true;
						else{ /*do nothing*/}
					}else if(novel_var.var_type == Position_depth_item_Novel_DEL){//del
						is_same_var_novel_to_known = true;
					}else{//INS
						const char * known_alt = known_var.get_alt();
						is_same_var_novel_to_known = true;
						for(int i = 0; i < novel_var_len + 1 ; i++){
							if("ACGT"[novel_var.alt_bin[i]] != known_alt[i]){
								is_same_var_novel_to_known = false;
							}
						}
					}
				}
				if(is_same_var_novel_to_known){
					known_var_info.depth += novel.variant_supp[novel_idx].var_supp;
					//NOTE: after the combination, the depth of KNOWN REF -= the depth of novel VAR
					//this is because the "depth of KNOWN REF" is ONLY count for the REF (exclude the depth of VARs)
					known.update_ref_supp_novel( - novel.variant_supp[novel_idx].var_supp);
					is_combined = true;
					break;
				}
			}
			if(is_combined){
				novel.variant_supp.erase(novel.variant_supp.begin() + novel_idx);
				novel_idx --;
			}
		}
	}

public:

	struct Haplotype_with_data{
		Haplotype_with_data(){}
		std::string hap_string;
		Window_t *window = NULL;
		int hap_id;

		void set(Window_t *window, std::string &hap_str, int hap_id){
			this->window = window;
			std::swap(hap_str, this->hap_string);
			this->hap_id = hap_id;
		}

	};
	std::vector<Haplotype_with_data> h_with_data;
	void store_hap_with_data_known(Window_t &window, Hap_counter &hap_counter){
		std::string ref_string, hap_string_buff;
		ref_h->load_ref_from_buff(window.chr_ID, window.st_pos, 300, ref_string);
		//store the reference into the target list
		//h_with_data.emplace_back(); h_with_data.back().set(&window, ref_string, -2); ref_string = h_with_data.back().hap_string;//copy back
		std::vector<hap_t> &c_hap_list = window.hap_list;
		uint hap_size = c_hap_list.size();
		for(uint hap_ID = 0; hap_ID < hap_size; hap_ID++){
			if(hap_counter.is_hap_with_read(hap_ID)){
				//uint32_t pre_hap_covarge = counter_p.get_depth(hap_ID, c_hap_list[hap_ID].hap_length/2);
				//if(pre_hap_covarge >= 0){
				window.get_hap_string(hap_ID, ref_string, hap_string_buff);
				h_with_data.emplace_back();
				h_with_data.back().set(&window, hap_string_buff, (hap_ID == hap_size-1)?-2:hap_ID);
			}
		}
	}

	void store_hap_with_data_novel(Candidate_result_Novel * cur_novel, std::string &contig_core){
		std::string &ori_contig = cur_novel->ass_b->contigs[0].seq;
		uint32_t* bam_cigar;
		int n_cigar = cur_novel->c_a->get_cigar(&bam_cigar);
		int contig_st_pos_in_ref = cur_novel->d_c.final_read_suggest_alignment_start_pos;
		int target_index = contig_st_pos_in_ref;
		int query_index = 0;
		int min_var_position = 150;//todo
		int max_var_position = 300-1;//todo
		int INDEL_ADJUST = 0;

		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:
				target_index += c_length; query_index += c_length;	break;//M
			case 1:	//insertion
				if(target_index >= min_var_position && target_index < max_var_position)
					INDEL_ADJUST+=c_length;
				query_index += c_length; break;//I, int chr_ID; int
			case 2:	 //Deletions
				if(target_index >= min_var_position && target_index < max_var_position)
					INDEL_ADJUST -= c_length;
				target_index += c_length; break;
			case 3:	target_index += c_length; query_index += c_length; break;//N, print N
			case 4:	target_index += c_length; break;//S, print -
			default: break;
			}
		}
		int load_bg = min_var_position - contig_st_pos_in_ref; load_bg = MAX(load_bg, 0); load_bg = MIN((int)ori_contig.size() - 1, load_bg);
		int load_ed = max_var_position - contig_st_pos_in_ref + INDEL_ADJUST; load_ed = MIN((int)ori_contig.size(), load_ed);
		contig_core = ori_contig.substr(load_bg, load_ed - load_bg);
	}

	//other calling functions:: support ONLY when Novel Var were detected:
	bool check_novel_haplotype_same_as_known_haplotype(
			std::vector<Window_t> &window_info,
			uint32_t pre_window_ID,
			uint32_t cur_window_ID,
			Hap_counter &pre_counter,
			Hap_counter &cur_counter,
			Candidate_result_Novel * cur_novel){
		//pre_counter.get_depth(hap_ID, POS);
		//S1: generate all haplotype list:
		//clear the list
		h_with_data.clear();
		//store list for known haplotypes
		store_hap_with_data_known(window_info[pre_window_ID], pre_counter);
		store_hap_with_data_known(window_info[cur_window_ID], cur_counter);
		//store the list for the Novel
		if(!cur_novel->is_assembly)//assembly the data when ASS is not running
			cur_novel->running_ASS(log_f, true);
		if(cur_novel->with_contig){
			std::string contig_core;
			store_hap_with_data_novel(cur_novel, contig_core);
			fprintf(log_f, "Novel contig check BG: %s\n", contig_core.c_str());
			//check:
			for(Haplotype_with_data &h : h_with_data){
				fprintf(log_f, "Novel contig check: with: st_pos %d, %s , hap_id %d\n", h.window->st_pos, h.hap_string.c_str(), h.hap_id);
				int contig_realigned_pos = h.hap_string.find(contig_core);
				if(contig_realigned_pos != -1){
					fprintf(log_f, "Novel contig repeat check SUCCESSFUL\n");
					NOVEL_VCF_BUFFs.clear();
					return true;
				}
			}
		}
		return false;
	}


	bool vcf_generatig(
			std::map<uint32_t, Position_depth_item> &pre_var_depth,
			std::map<uint32_t, Position_depth_item> &cur_var_depth,
			std::vector<Window_t> &window_info,
			uint32_t pre_window_ID,
			uint32_t cur_window_ID,
			Hap_counter &pre_counter,
			Hap_counter &cur_counter,
			Candidate_result_Novel * cur_novel
		){
		//INIT:
		//for novel vars:
		std::map<uint16_t, Position_depth_item_Novel> & novel_var = cur_novel->get_candidate_var_positions();
		if(false && cur_window_ID == 546780){
			fprintf(stderr, " ");
		}
		//=========================================================================================================================
		//S1: Novel : try to combine with known variants
		//=========================================================================================================================
		if(!novel_var.empty()){
			//merge the SAME VARs in the novel VAR list and known VAR list
			//when a VAR is truly NOVEL, it kept in the original NOVEL VAR list
			//combine the novel and the known variants
			for(std::map<uint16_t, Position_depth_item_Novel>::iterator novel_it = novel_var.begin(); novel_it != novel_var.end(); novel_it++){
				//try to combine
				//S1: get the var position in the previous WB
				uint32_t novel_pos = novel_it->first;
				if(novel_pos < 150)
					continue;
				//try to combine
				//S1: get the var position in the current(next) WB
				uint32_t cur_pos = novel_pos - 150;
				//S2: find the var in pre/cur_var_depth
				auto cur_it = cur_var_depth.find(cur_pos);
				if(cur_it != cur_var_depth.end()){
					//try to combine
					combine_novel_to_known(cur_it->second, novel_it->second, cur_pos, window_info);
				}else{
					uint32_t pre_pos = novel_pos;
					auto pre_it = pre_var_depth.find(pre_pos);
					if(pre_it != pre_var_depth.end()){
						//try to combine
						combine_novel_to_known(pre_it->second, novel_it->second, pre_pos, window_info);
					}
				}
			}
		}
		//=========================================================================================================================
		//S2: Known VARs: generate
		//=========================================================================================================================
		//update the depth of known(CUR) using the novel depth
		if(!cur_var_depth.empty() && cur_novel->get_max_gap_read_depth() != 0){
			for(std::map<uint32_t, Position_depth_item>::iterator cur_it = cur_var_depth.begin(); cur_it != cur_var_depth.end(); cur_it++){
				uint32_t cur_pos = cur_it->first;
				uint32_t novel_pos = cur_pos + 150;
				uint32_t novel_hap_covarge = cur_novel->get_ref_depth_at_position(novel_pos);
				if(novel_hap_covarge >= MIN_CHANGE_DETAIL_DEPTH_NOVEL_VAR && novel_pos >= 150 && novel_pos < 300)
					cur_it->second.update_ref_supp_novel(novel_hap_covarge);
			}
		}

		//vcf buffs
		KNWON_VCF_BUFFs.clear();
		//search for the var depth belonged to cur and both
		if(!cur_var_depth.empty()){
			for(std::map<uint32_t, Position_depth_item>::iterator curr_it = cur_var_depth.begin(); curr_it != cur_var_depth.end(); curr_it++){
				uint32_t curr_pos = curr_it->first;
				if(curr_pos >= 150)
					continue;
				//try to combine
				//S1: get the var position in the previous WB
				uint32_t pre_pos = curr_pos + 150;
				//S2: find the var in pre_var_depth
				auto pre_it = pre_var_depth.find(pre_pos);
				if(pre_it != pre_var_depth.end()){
					//try to combine
					curr_it->second.combine(pre_it->second, window_info);
					pre_it->second.already_combined = true; //set the combined flag, but not truely erase it."//cur_var_depth.erase(cur_it);"
					curr_it->second.already_combined = true;
				}else{
					//when not find var, try to combine the depth of haplotypes in the PRE WB
					uint pre_hap_size = window_info[pre_window_ID].hap_list.size();
					for(uint hap_ID = 0; hap_ID < pre_hap_size; hap_ID++){
						uint32_t pre_hap_covarge = pre_counter.get_depth(hap_ID, pre_pos);
						if(pre_hap_covarge != 0)
							curr_it->second.update_ref_supp(pre_window_ID, hap_ID, pre_hap_covarge);
					}
				}
				//basic depth check
				//STORE VCFS
				curr_it->second.store_in_vcf(out_f, window_info, KNWON_VCF_BUFFs, MAX_LOW_depth);
			}
		}

		//update the depth of known(PRE) using the novel depth
		if(!pre_var_depth.empty() && cur_novel->get_max_gap_read_depth() != 0){
			for(std::map<uint32_t, Position_depth_item>::iterator pre_it = pre_var_depth.begin(); pre_it != pre_var_depth.end(); pre_it++){
				uint32_t pre_pos = pre_it->first;
				uint32_t novel_pos = pre_pos;
				uint32_t novel_hap_covarge = cur_novel->get_ref_depth_at_position(novel_pos);
				if(novel_hap_covarge >= MIN_CHANGE_DETAIL_DEPTH_NOVEL_VAR && novel_pos >= 150 && novel_pos < 300)
					pre_it->second.update_ref_supp_novel(novel_hap_covarge);
			}
		}

		//search for the var depth belonged to PRE only:
		if(!pre_var_depth.empty()){
			for(std::map<uint32_t, Position_depth_item>::iterator pre_it = pre_var_depth.begin(); pre_it != pre_var_depth.end(); pre_it++){
				if(pre_it->second.already_combined == true)
					continue;
				uint32_t pre_pos = pre_it->first;
				if(pre_pos < 150)
					continue;
				//get overall depth:
				uint32_t curr_pos = pre_pos - 150;
				//S2: find the var in cur_var_depth
				uint cur_hap_size = window_info[cur_window_ID].hap_list.size();
				for(uint hap_ID = 0; hap_ID < cur_hap_size; hap_ID++){
					uint32_t cur_hap_covarge = cur_counter.get_depth(hap_ID, curr_pos);
					if(cur_hap_covarge != 0)
						pre_it->second.update_ref_supp(cur_window_ID, hap_ID, cur_hap_covarge);
				}
				pre_it->second.store_in_vcf(out_f, window_info, KNWON_VCF_BUFFs, MAX_LOW_depth);
			}
		}

		//=========================================================================================================================
		//S3 KNWON: check repeat
		//=========================================================================================================================
		//Repeat and low quality VAR/REF support check;
		//some conditions trigger the check,
		//for example: (1)WB is in duplication region;
		//(2) tow variants is same, but write in different format
		//(3) there is too many vars in one WB
		//(4) the number of reads support ALT or REF is low, but not low enough to make sure what`s the true GT
		//(4.1) for example, when reads support REF is 0 or 1., the GT is 1/1; but when the REF is 2 or 3 while the ALT is 20,
		//(4.2) we should check the 2 reads that support REF
		if(!KNWON_VCF_BUFFs.empty()){
			//low quality VARs existence check
			std::vector<VCF_item *> low_qualty_var_list; low_qualty_var_list.clear();
			std::vector<VCF_item *> low_qualty_REF_list; low_qualty_REF_list.clear();
			if(KNWON_VCF_BUFFs.size() >= 4){
				fprintf(log_f, "TOO MANY VARS %ld ST\n", KNWON_VCF_BUFFs.size());
				//output all results
				for(VCF_item & vi:KNWON_VCF_BUFFs){
					vi.VCF_dump(log_f);
				}
				fprintf(log_f, "TOO MANY VARS %ld END \n", KNWON_VCF_BUFFs.size());
			}
		//	if(window_info[pre_window_ID].repeat_region_len > 0 || window_info[cur_window_ID].repeat_region_len > 0){
			bool repeat_1 = window_info[pre_window_ID].repeat_region_len > 30 && window_info[pre_window_ID].repeat_region_len < 200;
			bool repeat_2 = window_info[cur_window_ID].repeat_region_len > 30 && window_info[cur_window_ID].repeat_region_len < 200;
			bool with_repeat_var = with_DUPLICATE_VARs(KNWON_VCF_BUFFs);
			if(with_repeat_var){// when with repeat var, skip check for ALT
				for(VCF_item & vi:KNWON_VCF_BUFFs){
					if(vi.ref_supp >= MIN_depth_of_low_REF_check && vi.ref_supp <= MAX_depth_of_low_check_LV2)
						low_qualty_REF_list.emplace_back(&vi);
				}
			}
			else if(repeat_2 || repeat_1){ //when with repeat , check ALL INDELs
				for(VCF_item & vi:KNWON_VCF_BUFFs){
					if(vi.var_len != 0){
						if(vi.var_supp >= MIN_depth_of_low_ALT_check)
							low_qualty_var_list.emplace_back(&vi);
						if(vi.ref_supp >= MIN_depth_of_low_REF_check)
							low_qualty_REF_list.emplace_back(&vi);
					}else{
						if(vi.var_supp >= MIN_depth_of_low_ALT_check && vi.var_supp <= MAX_depth_of_low_check_LV2)
							low_qualty_var_list.emplace_back(&vi);
						if(vi.ref_supp >= MIN_depth_of_low_REF_check && vi.ref_supp <= MAX_depth_of_low_check_LV2)
							low_qualty_REF_list.emplace_back(&vi);
					}
				}
			}
			else if(KNWON_VCF_BUFFs.size() >= 4){
				for(VCF_item & vi:KNWON_VCF_BUFFs){
					if(vi.var_supp >= MIN_depth_of_low_ALT_check && vi.var_supp <= MAX_depth_of_low_check_LV2)
						low_qualty_var_list.emplace_back(&vi);
					if(vi.ref_supp >= MIN_depth_of_low_REF_check && vi.ref_supp <= MAX_depth_of_low_check_LV2)
						low_qualty_REF_list.emplace_back(&vi);
				}
			}else{
				for(VCF_item & vi:KNWON_VCF_BUFFs){
					if(vi.var_len != 0){
						if(vi.var_supp >= MIN_depth_of_low_ALT_check && vi.var_supp <= MAX_depth_of_low_check)
							low_qualty_var_list.emplace_back(&vi);
						if(vi.ref_supp >= MIN_depth_of_low_REF_check && vi.ref_supp <= MAX_depth_of_low_check_LV2)
							low_qualty_REF_list.emplace_back(&vi);
					}else{
						if(vi.var_supp >= MIN_depth_of_low_ALT_check && vi.var_supp <= MAX_depth_of_low_check)
							low_qualty_var_list.emplace_back(&vi);
						if(vi.ref_supp >= MIN_depth_of_low_REF_check && vi.ref_supp <= MAX_depth_of_low_check)
							low_qualty_REF_list.emplace_back(&vi);
					}
				}
			}

			//low quality VARs realignment
			if(!low_qualty_var_list.empty()){
				std::map<uint32_t, Position_depth_item> * pos_depth_map_p[2] = {&pre_var_depth, &cur_var_depth};
				Hap_counter *counter_p[2] = {&pre_counter, &cur_counter};
				low_quality_var_realignment( pos_depth_map_p,
						window_info, pre_window_ID, cur_window_ID,
						counter_p, low_qualty_var_list);
			}

			//low quality REFs realignment
			if(!low_qualty_REF_list.empty()){
				std::map<uint32_t, Position_depth_item> * pos_depth_map_p[2] = {&pre_var_depth, &cur_var_depth};
				Hap_counter *counter_p[2] = {&pre_counter, &cur_counter};
				low_quality_ref_realignment( pos_depth_map_p,
						window_info, pre_window_ID, cur_window_ID,
						counter_p, low_qualty_REF_list);
			}
			//duplicated VARs check:
			//like:
			//chr1  59347102  ID_KNOWN  C CTTTTT  30  PASS  TYPE=INS  GT:DP 0/1:39,21
			//chr1  59347103  ID_KNOWN  T TTTTTT  30  PASS  TYPE=INS  GT:DP 0/1:26,35
			//OR:
			//chr1  61554944  ID_KNOWN  G GCTC  30  PASS  TYPE=INS  GT:DP 0/1:21,45
			//chr1  61554947  ID_KNOWN  C CCTC  30  PASS  TYPE=INS  GT:DP 0/1:49,21
			if(with_repeat_var)
				DUPLICATE_VARs_check(KNWON_VCF_BUFFs);
		}
		bool with_var = (NOVEL_VCF_BUFFs.size() > 0 || KNWON_VCF_BUFFs.size() > 0);
		if(false && (with_var)){
			fprintf(out_f, "####CHECH BG\n");
		}

		//=========================================================================================================================
		//S4 NOVEL VAR: check repeat:
		//=========================================================================================================================
		//NOVEL VAR: check whether with NOVEL vars that CANN`T be combined with known VARs (that is: truly NOVEL)
		NOVEL_VCF_BUFFs.clear();
		bool with_novel_var_that_cannt_be_combined = false;
		for(std::map<uint16_t, Position_depth_item_Novel>::iterator novel_it = novel_var.begin(); novel_it != novel_var.end(); novel_it++){
			if(novel_it->second.variant_supp.empty() == false)
				with_novel_var_that_cannt_be_combined = true;
		}
		if(with_novel_var_that_cannt_be_combined){
			bool novel_is_same_as_hap = check_novel_haplotype_same_as_known_haplotype(window_info, pre_window_ID, cur_window_ID, pre_counter, cur_counter, cur_novel);
			//condition: neither the contig is not same as any other haplotype nor the Known haplotype didnot call any results.
			bool novel_is_true = (!novel_is_same_as_hap) || (KNWON_VCF_BUFFs.empty());
			//condtion: when the novel is true and with novel var called by Assembly.
			if(novel_is_true && cur_novel->is_assembly && cur_novel->with_contig && !novel_var.empty()){
				//update the depth of novel using the known depth
				for(std::map<uint16_t, Position_depth_item_Novel>::iterator novel_it = novel_var.begin(); novel_it != novel_var.end(); novel_it++){
					if(novel_it->second.variant_supp.empty())
						continue;
					uint32_t novel_var_position = novel_it->first;
					{ //adding the depth of PRE
						int pre_pos = novel_var_position;
						if(pre_pos >= 0 && pre_pos < 300){
							uint pre_hap_size = window_info[pre_window_ID].hap_list.size();
							for(uint hap_ID = 0; hap_ID < pre_hap_size; hap_ID++){
								uint32_t pre_hap_covarge = pre_counter.get_depth(hap_ID, pre_pos);
								if(pre_hap_covarge != 0)
									novel_it->second.add_ref_supp(pre_hap_covarge);
							}
						}
					}
					{//adding the depth of CUR
						int cur_pos = novel_var_position - 150;
						if(cur_pos >= 0 && cur_pos < 300){
							uint cur_hap_size = window_info[cur_window_ID].hap_list.size();
							for(uint hap_ID = 0; hap_ID < cur_hap_size; hap_ID++){
								uint32_t cur_hap_covarge = cur_counter.get_depth(hap_ID, cur_pos);
								if(cur_hap_covarge != 0)
									novel_it->second.add_ref_supp(cur_hap_covarge);
							}
						}
					}
					{//adding the depth of NOVEL
						uint32_t total_novel_hap_covarge = cur_novel->get_total_depth_at_position(novel_var_position);
						if(total_novel_hap_covarge != 0)
							novel_it->second.add_ref_supp(total_novel_hap_covarge);
					}
					//STORE VCFS
					novel_it->second.store_in_vcf(out_f, window_info, NOVEL_VCF_BUFFs, MAX_LOW_depth, ref_h);
				}
			}
			//remove VARs that already included in Known list;
			//sometimes the haplotype was not same, but the var was same.
			for(VCF_item & vi_n:NOVEL_VCF_BUFFs){
				bool is_same_with_known_var = false;
				for(VCF_item & vi_k:KNWON_VCF_BUFFs){
					if(vi_n.var_is_same(vi_k))
						is_same_with_known_var = true;
				}
				if(is_same_with_known_var){
					vi_n.var_supp = 0;
					if(true) { fprintf(log_f, "Known var is removed because it is same as known vars\n"); vi_n.VCF_dump(log_f);}
				}
			}
		}

		//=========================================================================================================================
		//S5 OUTPUT
		//=========================================================================================================================
		//output all results: Novel
		for(VCF_item & vi:NOVEL_VCF_BUFFs){
			if(vi.var_supp != 0){
				//debug code
				if(true) {fprintf(log_f, "DEBUG\t"); vi.VCF_dump(log_f);}
				sprintf(var_name_buff, "%s_novel_%d",sample_name_string, novel_var_ID++);
				vi.VCF_dump_core(out_f, var_name_buff);
			}
		}
		//output all results: Known
		for(VCF_item & vi:KNWON_VCF_BUFFs){
			if(vi.var_supp != 0){
				if(true){ fprintf(log_f, "DEBUG\t"); vi.VCF_dump(log_f); }
				sprintf(var_name_buff, "%s_known_%d", sample_name_string, known_var_ID++);
				vi.VCF_dump_core(out_f, var_name_buff);
			}
		}

		if(false && (with_var)){
			fprintf(out_f, "####CHECH ED\n");
		}
		//=========================================================================================================================
		// return
		//=========================================================================================================================
		return with_var;
	}

private:
	void load_hap_list_with_vars(Window_t & window_p, std::string &ref_string, std::string &hap_string_buff, std::vector<Compare_main_hap> &compare_main_hap_list){
		ref_h->load_ref_from_buff(window_p.chr_ID, window_p.st_pos, 300, ref_string);
		//store the reference into the target list
		compare_main_hap_list.clear(); compare_main_hap_list.emplace_back(-2, ref_string); ref_string = compare_main_hap_list.back().str;
		std::vector<hap_t> &c_hap_list = window_p.hap_list;
		uint hap_size = c_hap_list.size();
		for(uint hap_ID = 0; hap_ID < hap_size; hap_ID++){
			//uint32_t pre_hap_covarge = counter_p.get_depth(hap_ID, c_hap_list[hap_ID].hap_length/2);
			//if(pre_hap_covarge >= 0){
			window_p.get_hap_string(hap_ID, ref_string, hap_string_buff);
			compare_main_hap_list.emplace_back(hap_ID, hap_string_buff);
			//}
		}
	}

	int search_var(Window_t & w, int search_pos, Variant & old_var){
		int search_var_ID = -1;
		for(uint var_ID = 0; var_ID < w.var_list.size(); var_ID++){
			auto & search_var = w.var_list[var_ID];
			if(search_var.ref_pos != search_pos)	continue;
			if(search_var.isSAME_ref_alt(old_var)) { search_var_ID = var_ID; break; }
		}
		return search_var_ID;
	}

	bool is_read_cover_vars(Window_t *w, int hap_ID, int varID, int read_st_pos, int read_ed_pos){
		uint var_st, var_ed;
		var_covered_by_region(w, hap_ID, read_st_pos, read_ed_pos, var_st, var_ed);
		return read_cover_vars_core(w, hap_ID, var_st, var_ed,  varID);
	}

	void low_quality_var_realignment(
			std::map<uint32_t, Position_depth_item> * pos_depth_map_p[2],
			std::vector<Window_t> &window_info,
			uint32_t pre_window_ID,
			uint32_t cur_window_ID,
			Hap_counter *counter_p[2],
			std::vector<VCF_item *> &low_qualty_var_list
			){
		//loading reference
		Window_t *window_p[2] = { &(window_info[pre_window_ID]), &(window_info[cur_window_ID]) };
		uint32_t wb_ID[2] = {pre_window_ID, cur_window_ID};
		int varID[2];
		std::vector<Compare_main_hap> compare_main_hap_list[2];
		std::string ref_string[2];
		std::string hap_string_buff;
		std::string read_buff;
		//buffs
		std::vector<int> read_pos_list;
		//loading haplotype list
		for(int i = 0; i < 2; i++)
			load_hap_list_with_vars(*(window_p[i]), ref_string[i], hap_string_buff, compare_main_hap_list[i]);
		//check the reads supp the low quality, whether they support the reference or the main haplotype?
		for(VCF_item * vi:low_qualty_var_list){
			//get support haplotype:
			int belong_to_other_read_num = 0;
			int not_belong_to_other_read_num = 0;
			//get VAR ID:
			if(vi->var_wb == wb_ID[0]){
				varID[0] = vi->var_ID;
				Variant &c_var = window_info[vi->var_wb].var_list[vi->var_ID];
				varID[1] =  search_var(*(window_p[1]), c_var.ref_pos - 150, c_var);
			}else{
				varID[1] = vi->var_ID;
				Variant &c_var = window_info[vi->var_wb].var_list[vi->var_ID];
				varID[0] =  search_var(*(window_p[0]), c_var.ref_pos + 150, c_var);
			}

			//FOR each wrong VARs:
			for(int i = 0; i < 2; i++){
				//get the support haplotypes of low cover variants: the low quality vars may get support from both of the WB; i == 0: the first (pre)WB; i == 1: the second (NEXT) hap
				Position_depth_item * pos_depth = NULL;
				{
					//In condition(1), the pos of first WB is needed (pos += 150), other conditions, get pos_var directly
					uint32_t var_pos_wb = window_info[vi->var_wb].var_list[vi->var_ID].ref_pos;
					if(i == 0 && vi->var_wb != wb_ID[i]) var_pos_wb += 150;
					//get the Position_depth_item for this WB and position.
					{ auto it = pos_depth_map_p[i]->find(var_pos_wb); pos_depth = &(it->second); if(it == pos_depth_map_p[i]->end()) continue; }
				}
				//search support hap and read list using [Var_supp_haplotype_item] and var_ID
				int LAST_search_hap_idx_BUFF = -1;
				for(Var_supp_haplotype_item & var_hap:pos_depth->vh_l){
					if(var_hap.var_ID != varID[i]) continue;
					if(var_hap.wb_ID != wb_ID[i]) continue;
					//when find, get the haplotype string using hapID
					window_p[i]->get_hap_string(var_hap.hap_ID, ref_string[i], hap_string_buff);
					if(false) fprintf(log_f, "ALL_hap: %s\nREF: %s\n", hap_string_buff.c_str(),  ref_string[i].c_str());
					//get all support reads that aligned to this haplotype, all reads will be re-alignment
					counter_p[i]->get_hap_read_list(var_hap.hap_ID, read_pos_list);
					for(uint read_id = 0; read_id < read_pos_list.size(); read_id++){
						int read_st_pos  = read_pos_list[read_id];
						//S1: //read cover var check:
						if(false == is_read_cover_vars(window_p[i], var_hap.hap_ID, varID[i], read_st_pos, read_st_pos + 150))
							continue;
						//S2: get read string: //get read string using read_st_pos
						read_buff = hap_string_buff.substr(read_st_pos, 150);//todo::read length
						if(false) fprintf(log_f, "read @ read_st_pos %d: %s\n",read_st_pos, read_buff.c_str());
						//S3: search reads in other haplotypes
						int search_hap_idx = search_read_in_compare_main_hap_list(compare_main_hap_list[i], read_buff, window_p[i],  varID[i], false, LAST_search_hap_idx_BUFF);
						//check read whether support the reference??
						if(search_hap_idx == -1)	not_belong_to_other_read_num++; //if(false) fprintf(log_f, "Read of ID %d,with original pos %d, wb_ID %d, cannot be realigned\n", read_id, read_st_pos, wb_ID[i]);
						else						{belong_to_other_read_num++; fprintf(log_f, "@%d ", search_hap_idx);}
					}
				}
			}
			//re-genotyping
			if(true){
				fprintf(log_f, "LOW VAR re-alignment: belong_to_var: %d, belong_to_other: %d original: %d\t", not_belong_to_other_read_num, belong_to_other_read_num, vi->var_supp);
				vi->var_supp = not_belong_to_other_read_num;
				vi->VCF_dump(log_f);
				get_GT_type_KNOWN(vi->var_supp, vi->ref_supp, vi->GT_TYPE , vi->PASS_TYPE, MAX_LOW_depth);
			}
		}
	}

	//given a hap ID, a st and ed pos in the hap, return which var is covered by the region
		//input: skip
		//output: the start index of var in haplotype.
		//If the region covered the 1st,2ed,3rd of vars of this hap ,return 0,2. if only cover 1st,return 0,0; if only cover nothing,return MAX_INT32,0;
		void var_covered_by_region(Window_t * wb, int hap_id, int st, int ed, uint &var_st, uint &var_ed){
			var_st = MAX_int32t; var_ed = 0;
			std::vector<uint16_t> &var_l = wb->hap_list[hap_id].variant_list;
			int c_INDLE_OFFSET = 0;
			int global_INDLE_OFFSET = 0;
			for(uint i = 0; i < var_l.size(); i++){
				Variant &c_var = wb->var_list[var_l[i]];
				c_INDLE_OFFSET = c_var.get_length();///####
				uint32_t var_pos = c_var.ref_pos + global_INDLE_OFFSET;
				if(var_pos + ABS(c_INDLE_OFFSET) >= st && var_pos < ed + ABS(c_INDLE_OFFSET)){
					var_st = MIN(var_st, i);
					var_ed = MAX(var_ed, i);
				}
				global_INDLE_OFFSET += c_INDLE_OFFSET;///####
			}
		}

	void low_quality_ref_realignment(
			std::map<uint32_t, Position_depth_item> * pos_depth_map_p[2],
			std::vector<Window_t> &window_info,
			uint32_t pre_window_ID,
			uint32_t cur_window_ID,
			Hap_counter *counter_p[2],
			std::vector<VCF_item *> &low_qualty_var_list){
		//return;
		//loading reference
		Window_t *window_p[2] = { &(window_info[pre_window_ID]), &(window_info[cur_window_ID]) };
		uint32_t wb_ID[2] = {pre_window_ID, cur_window_ID};
		int varID[2];
		std::vector<Compare_main_hap> compare_main_hap_list[2];
		std::string ref_string[2];
		std::string hap_string_buff;
		std::string read_buff;
		//buffs
		std::vector<int> read_pos_list;

		//loading haplotype list
		for(int i = 0; i < 2; i++)
			load_hap_list_with_vars(*(window_p[i]), ref_string[i], hap_string_buff, compare_main_hap_list[i]);

		//check the reads supp the low quality, whether they support the reference or the main haplotype?
		for(VCF_item * vi:low_qualty_var_list){
			//get support haplotype:
			int belong_to_other_read_num = 0;
			int not_belong_to_other_read_num = 0;
			//get VAR ID:
			if(vi->var_wb == wb_ID[0]){
				varID[0] = vi->var_ID;
				Variant &c_var = window_info[vi->var_wb].var_list[vi->var_ID];
				varID[1] =  search_var(*(window_p[1]), c_var.ref_pos - 150, c_var);
			}else{
				varID[1] = vi->var_ID;
				Variant &c_var = window_info[vi->var_wb].var_list[vi->var_ID];
				varID[0] =  search_var(*(window_p[0]), c_var.ref_pos + 150, c_var);
			}

			Position_depth_item * pos_depth[2] = {0};
			for(int i = 0; i < 2; i++){
				//get the support haplotypes of low cover variants: the low quality vars may get support from both of the WB; i == 0: the first (pre)WB; i == 1: the second (NEXT) hap
				//In condition(1), the pos of first WB is needed (pos += 150), other conditions, get pos_var directly
				uint32_t var_pos_wb = window_info[vi->var_wb].var_list[vi->var_ID].ref_pos;
				if(i == 0 && vi->var_wb != wb_ID[i]) var_pos_wb += 150;
				//get the Position_depth_item for this WB and position.
				auto it = pos_depth_map_p[i]->find(var_pos_wb); pos_depth[i] = &(it->second);
				if(it == pos_depth_map_p[i]->end())
					pos_depth[i] = NULL;
			}
			//bool both_pos_depth_exist = (pos_depth[0] != NULL && pos_depth[1] != NULL);

			//FOR each wrong VARs:
			for(int i = 0; i < 2; i++){
				if(pos_depth[i] == NULL) continue;
				//search support hap and read list using [Var_supp_haplotype_item] and var_ID
				int LAST_search_hap_idx_BUFF = -1;
				for(Var_supp_haplotype_item & var_hap:pos_depth[i]->vh_l){
					int c_i = i;
					if(var_hap.wb_ID != wb_ID[i]){
						if(pos_depth[i]->already_combined) 	continue;
						else if(var_hap.wb_ID == wb_ID[1 - c_i])
							c_i = 1 - c_i;
					}
					if(var_hap.var_ID == varID[c_i]) continue;//#3
					//when find, get the haplotype string using hapID
					window_p[c_i]->get_hap_string(var_hap.hap_ID, ref_string[c_i], hap_string_buff);
					if(false) fprintf(log_f, "ALL_hap: %s\nREF: %s\n", hap_string_buff.c_str(),  ref_string[c_i].c_str());
					//get all support reads that aligned to this haplotype, all reads will be re-alignment
					counter_p[c_i]->get_hap_read_list(var_hap.hap_ID, read_pos_list);
					for(uint read_id = 0; read_id < read_pos_list.size(); read_id++){
						int read_st_pos  = read_pos_list[read_id];
						//skip read near the edge of WB
						if(vi->var_len < 0){
							if(read_st_pos + 150 - vi->var_len >= (int)hap_string_buff.size()) continue; //M2
						}

						//S0 + S1: read cover the position, but not cover the VAR, that mean it may support the REF.
						//S0: POS check, ONLY USED IN REF:  //#4, that means that the read cover the position of VARs, if not cover the position, skip
						{
							uint var_pos = window_p[c_i]->var_list[varID[c_i]].ref_pos;
							if(!REF_pos_covered_by_ALT_region(window_p[c_i], var_hap.hap_ID, read_st_pos, read_st_pos +150, var_pos))
								continue;
						}
						//S1: //read cover var check: that means that the read NOT cover the of VARs, if cover the VAR, skip
						if(true == is_read_cover_vars(window_p[c_i], var_hap.hap_ID, varID[c_i], read_st_pos, read_st_pos + 150)) //#2
							continue;
						//S2: get read string: //get read string using read_st_pos
						read_buff = hap_string_buff.substr(read_st_pos, 150);//todo::read length
						if(false) fprintf(log_f, "read @ read_st_pos %d: %s\n",read_st_pos, read_buff.c_str());
						//S3: search reads in other haplotypes
						int search_hap_idx = search_read_in_compare_main_hap_list(compare_main_hap_list[c_i], read_buff, window_p[c_i],  varID[c_i], true, LAST_search_hap_idx_BUFF);//#1
						//check read whether support the reference??
						if(search_hap_idx == -1)	not_belong_to_other_read_num++; //if(false) fprintf(log_f, "Read of ID %d,with original pos %d, wb_ID %d, cannot be realigned\n", read_id, read_st_pos, wb_ID[c_i]);
						else						belong_to_other_read_num++;
					}
				}
			}

			//re-genotyping
			if(true){
				fprintf(log_f, "LOW REF re-alignment: belong_to_REF: %d, belong_to_other: %d original: %d\t", not_belong_to_other_read_num, belong_to_other_read_num, vi->ref_supp);
				vi->ref_supp = not_belong_to_other_read_num;
				vi->VCF_dump(log_f);
				get_GT_type_KNOWN(vi->var_supp, vi->ref_supp, vi->GT_TYPE , vi->PASS_TYPE, MAX_LOW_depth);
			}
		}
	}
	bool REF_pos_covered_by_ALT_region(Window_t * wb, int hap_id, int st, int ed, int REF_POS){
		std::vector<uint16_t> &var_l = wb->hap_list[hap_id].variant_list;
		int global_INDLE_OFFSET = 0;
		for(uint i = 0; i < var_l.size(); i++){
			Variant &c_var = wb->var_list[var_l[i]];
			uint32_t var_pos = c_var.ref_pos + global_INDLE_OFFSET;
			if(var_pos > st){ break;}
			global_INDLE_OFFSET += c_var.get_length();///####
		}
		REF_POS += global_INDLE_OFFSET;
		return (st <= REF_POS && ed > REF_POS);
	}


	int search_read_in_compare_main_hap_list(std::vector<Compare_main_hap> & compare_main_hap_list, std::string & read_buff, Window_t *window_p, int var_ID, bool REF_check, int &LAST_search_hap_idx_BUFF){
		int search_hap_idx = -1;
		for(int hap_idx = -1;hap_idx < (int)compare_main_hap_list.size();hap_idx++){
			int c_hap_idx = hap_idx;
			if(hap_idx == -1){
				if(LAST_search_hap_idx_BUFF != -1 && LAST_search_hap_idx_BUFF <  (int)compare_main_hap_list.size())
					c_hap_idx = LAST_search_hap_idx_BUFF;
				else
					 continue;
			}
			int new_realigned_pos = compare_main_hap_list[c_hap_idx].str.find(read_buff);
			int new_realigned_hap_ID = compare_main_hap_list[c_hap_idx].id;
			if(new_realigned_pos != -1){
				//when the read was realigned to reference
				if(new_realigned_hap_ID == -2){
					//if(false) fprintf(log_f, "Read of ID %d,with original pos %d, realigned to wb_ID %d, @hap_idx.id %d @pos %d\n", read_id, read_st_pos, wb_ID[i], new_realigned_hap_ID, new_realigned_pos);
					if(!REF_check) search_hap_idx = -2;
					continue;
				}
				else{//when the read was realigned to the haplotype NOT REF
					//check whether the newly realigned pos is covered by the original vars
					uint var_st, var_ed;
					var_covered_by_region(window_p, new_realigned_hap_ID, new_realigned_pos, new_realigned_pos + read_buff.size(), var_st, var_ed);
					bool read_cover_var = read_cover_vars_core(window_p, new_realigned_hap_ID, var_st, var_ed, var_ID);
					if(REF_check && read_cover_var == false) continue;//IN REF check, when new realigned reads NOT cover its original support var, skip
					if(!REF_check && read_cover_var == true) continue;//when new realigned reads cover its original support var, skip
					search_hap_idx = c_hap_idx;
					LAST_search_hap_idx_BUFF = search_hap_idx;
					break;
				}
			}

			if(false)
				fprintf(log_f, "hap_idx %d, find_pos_ref %d hap %s\n",c_hap_idx, new_realigned_pos, compare_main_hap_list[c_hap_idx].str.c_str());
		}
		return search_hap_idx;
	}

	//it return whether a var(in var ID of WB) is covered by a region in a haplotype
	bool read_cover_vars_core(Window_t * wb, int hap_id, uint var_st, uint var_ed, int search_var_ID){
		for(uint i = var_st; i <= var_ed; i++){
			if(wb->hap_list[hap_id].variant_list[i] == search_var_ID)
				return true;
		}
		return false;
	}

	//check whether or not there are two vars that is exact the SAME
	//when two vars is SAME, return true, when no var is same, return false
	bool with_DUPLICATE_VARs(std::vector<VCF_item> &VCF_BUFFs){
		if(VCF_BUFFs.size() <= 1)
			return false;
		uint vcf_size = VCF_BUFFs.size();
		for(uint v1_id = 0; v1_id < vcf_size; v1_id ++){
			int v1_len = VCF_BUFFs[v1_id].var_len;
			if(v1_len == 0)
				continue;
			for(uint v2_id = v1_id + 1; v2_id < vcf_size; v2_id ++){
				if(VCF_BUFFs[v2_id].var_len != v1_len)
					continue;
				int pos_diff = ABS_U(VCF_BUFFs[v1_id].var_pos, VCF_BUFFs[v2_id].var_pos);
				if(pos_diff > 20)
					continue;
				if(VCF_BUFFs[v1_id].var_supp == 0|| VCF_BUFFs[v2_id].var_supp == 0)
					continue;
				//IS DUP VARs?
				VCF_item & v1 = VCF_BUFFs[v1_id];
				VCF_item & v2 = VCF_BUFFs[v2_id];
				int load_ref_st = MIN(v1.var_pos, v2.var_pos);
				int load_ref_ed = MAX(v1.var_pos, v2.var_pos) + 30;
				std::string ref_string;
				ref_h->load_ref_from_buff(v1.chr_ID, load_ref_st, load_ref_ed -load_ref_st, ref_string);
				std::string v1_string = ref_string;
				std::string v2_string = ref_string;
				v1_string.erase(v1.var_pos - load_ref_st, strlen(v1.ref_p));
				v1_string.insert(v1.var_pos - load_ref_st, v1.alt_p);
				v2_string.erase(v2.var_pos - load_ref_st, strlen(v2.ref_p));
				v2_string.insert(v2.var_pos - load_ref_st, v2.alt_p);
				if(v2_string == v1_string)
					return true;
			}
		}
		return false;
	}

	void DUPLICATE_VARs_check(std::vector<VCF_item> &VCF_BUFFs){
		uint vcf_size = VCF_BUFFs.size();
		for(uint v1_id = 0; v1_id < vcf_size; v1_id ++){
			int v1_len = VCF_BUFFs[v1_id].var_len;
			if(v1_len == 0)
				continue;
			for(uint v2_id = v1_id + 1; v2_id < vcf_size; v2_id ++){
				if(VCF_BUFFs[v2_id].var_len != v1_len)
					continue;
				int pos_diff = ABS_U(VCF_BUFFs[v1_id].var_pos, VCF_BUFFs[v2_id].var_pos);
				if(pos_diff > 20)
					continue;
				if(VCF_BUFFs[v1_id].var_supp == 0|| VCF_BUFFs[v2_id].var_supp == 0)
					continue;
				//IS DUP VARs?
				VCF_item & v1 = VCF_BUFFs[v1_id];
				VCF_item & v2 = VCF_BUFFs[v2_id];
				int load_ref_st = MIN(v1.var_pos, v2.var_pos);
				int load_ref_ed = MAX(v1.var_pos, v2.var_pos) + 30;
				std::string ref_string;
				ref_h->load_ref_from_buff(v1.chr_ID, load_ref_st, load_ref_ed -load_ref_st, ref_string);
				std::string v1_string = ref_string;
				std::string v2_string = ref_string;
				fprintf(log_f, "DUP VARs\n");
				fprintf(log_f, "REF: [%d, %d] %s\n", load_ref_st, load_ref_ed, ref_string.c_str());
				v1.VCF_dump(log_f);
				v2.VCF_dump(log_f);
				//var check:
				//v1:
				v1_string.erase(v1.var_pos - load_ref_st, strlen(v1.ref_p));
				v1_string.insert(v1.var_pos - load_ref_st, v1.alt_p);
				v2_string.erase(v2.var_pos - load_ref_st, strlen(v2.ref_p));
				v2_string.insert(v2.var_pos - load_ref_st, v2.alt_p);

				fprintf(log_f, "v1_string %s\n", v1_string.c_str());		//REF
				fprintf(log_f, "v2_string %s\n", v2_string.c_str());		//ALT
				if(v2_string == v1_string){
					fprintf(log_f, "String is SAME\n");
					//remove duplicate VARs
					v1.ref_supp -= v2.var_supp; v1.ref_supp = MAX(0,v1.ref_supp);
					v2.ref_supp -= v1.var_supp; v2.ref_supp = MAX(0,v2.ref_supp);
					v1.ref_supp = MIN(v1.ref_supp, v2.ref_supp);
					v1.var_supp += v2.var_supp;
					v2.var_supp = 0;
					get_GT_type_KNOWN(v1.var_supp, v1.ref_supp, v1.GT_TYPE , v1.PASS_TYPE, MAX_LOW_depth);
					get_GT_type_KNOWN(v2.var_supp, v2.ref_supp, v2.GT_TYPE , v2.PASS_TYPE, MAX_LOW_depth);
				}
				else
					fprintf(log_f, "String is NOT SAME\n");
				fprintf(log_f, "DUP VARs END\n");
			}
		}
	}
};


}

#endif /* SRC_BWT_ONLINE_VAR_CALLING_VAR_CALLER_ALL_TYPE_HPP_ */
