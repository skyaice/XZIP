#include "known_var_candidate_generator.hpp"
#include "GT_TYPE.hpp"

using namespace Aln_online;

void Position_depth_item::store_in_vcf(FILE* out, std::vector<Window_t> &window_info, std::vector<VCF_item> &VCF_BUFFs, int MAX_LOW_depth){
	uint32_t total_depth = ref_supp;
	for(var_supp_item & var_info : variant_supp)
		total_depth += var_info.depth;

	for(var_supp_item & var_info : variant_supp)
	{
		if(var_info.var_ID == MAX_uint32_t)		continue;
		if(var_info.depth == 1)					continue;
		//simple GT
		int GT_TYPE;// 0: 0/0; 1: 0/1; 2: 1/1
		int PASS_TYPE; // 0: NOT_PASS;1: PASS; 2:LOW_COV
		int32_t true_ref_supp = total_depth - var_info.depth;
		if(true_ref_supp < 0) true_ref_supp = 0;

		get_GT_type_KNOWN(var_info.depth, true_ref_supp, GT_TYPE, PASS_TYPE, MAX_LOW_depth);
		//uint32_t var_ID_this_wb, var_ID_pre_wb;
		Variant * c_var = &(window_info[var_info.wb_ID].var_list[var_info.var_ID]);
		uint32_t var_pos = c_var->ref_pos + window_info[var_info.wb_ID].st_pos;
		uint32_t chr_ID = window_info[var_info.wb_ID].chr_ID;

		if(c_var->isNoVar())
			continue;
		VCF_BUFFs.emplace_back();
		VCF_BUFFs.back().set(
			chr_ID, var_pos,
			c_var->get_ref(), c_var->get_alt(),
			30, PASS_TYPE, GT_TYPE,
			c_var->get_length(),
			true_ref_supp ,var_info.depth,
			var_info.var_ID,
			var_info.wb_ID,
			0
		);
	}
}

void Window_block_counter_Known_Var::init(std::vector<Window_t> &window_info, FILE* log_f_){
	set_new_wb(window_info, 0);
	BLANK_window_process_block.init(window_info[0].total_hap_length, 0, log_f);
	BLANK_window_process_block.realign_hap_list_final.emplace_back();
	BLANK_window_process_block.realign_hap_list_final.back().resize(1000, 0);
	log_f = log_f_;
}

void Window_block_counter_Known_Var::generate_known_var_candidate(std::vector<Window_t> &window_info){
	for(int l=0; l<4; l++)
		position_depth_counting(l, window_info);
	//generate vcfs
	uint32_t st_wb_ID = window_process_block[0].get_window_id();
	if(false){
		uint32_t debug_var_pos = 82034692;
		uint32_t debug_wb_ID = debug_var_pos/150;
		//log_f = stderr;
		if(st_wb_ID == debug_wb_ID/4*4){
			uint32_t debug_wb_st = 	debug_wb_ID*150+1;
			int var_pos_wb = debug_var_pos - debug_wb_st;
			fprintf(log_f, "debug_wb_ID %d, debug_wb_st %d, var_pos_wb %d\n", debug_wb_ID, debug_wb_st, var_pos_wb);
			fprintf(log_f, "window_info[debug_wb_ID].repeat_region_len %d\n", window_info[debug_wb_ID].repeat_region_len);
			//test:
			uint32_t c_st_wb_ID = debug_wb_ID;
			std::string hap_string;
			std::string ref;
			ref_h->load_ref_from_buff(window_info[c_st_wb_ID].chr_ID, window_info[c_st_wb_ID].st_pos, 300, ref);
			fprintf(log_f, "ref %s\n", ref.c_str());
			for(uint i = 0; i < window_info[c_st_wb_ID].hap_list.size(); i++){
				window_info[c_st_wb_ID].get_hap_string(i, ref, hap_string);
				fprintf(log_f, "i:%d  %s\n", i, hap_string.c_str());
			}
			window_process_block[debug_wb_ID - st_wb_ID].print_realign_hap_list_final(window_info[debug_wb_ID].hap_list.size(), log_f);		fprintf(log_f, "\n");
			window_info[debug_wb_ID].print(log_f);
			window_process_block[debug_wb_ID - st_wb_ID].print_realign_hap_list_ori(window_info[debug_wb_ID].hap_list.size(), log_f);			fprintf(log_f, "\n");
			for(int l=0; l<4; l++)
				position_depth_counting(l, window_info);

		}

		debug_wb_ID = (debug_var_pos-150)/150;
		if(st_wb_ID == debug_wb_ID/4*4){
			uint32_t debug_wb_st = 	debug_wb_ID*150+1;
			int var_pos_wb = debug_var_pos - debug_wb_st;
			fprintf(log_f, "debug_wb_ID %d, debug_wb_st %d, var_pos_wb %d\n", debug_wb_ID, debug_wb_st, var_pos_wb);
			fprintf(log_f, "window_info[debug_wb_ID].repeat_region_len %d\n", window_info[debug_wb_ID].repeat_region_len);
			//test:
			uint32_t c_st_wb_ID = debug_wb_ID;
			std::string hap_string;
			std::string ref;
			ref_h->load_ref_from_buff(window_info[c_st_wb_ID].chr_ID, window_info[c_st_wb_ID].st_pos, 300, ref);
			fprintf(log_f, "ref %s\n", ref.c_str());
			for(uint i = 0; i < window_info[c_st_wb_ID].hap_list.size(); i++){
				window_info[c_st_wb_ID].get_hap_string(i, ref, hap_string);
				fprintf(log_f, "i:%d  %s\n", i, hap_string.c_str());
			}
			window_process_block[debug_wb_ID - st_wb_ID].print_realign_hap_list_final(window_info[debug_wb_ID].hap_list.size(), log_f);		fprintf(log_f, "\n");
			window_info[debug_wb_ID].print(log_f);
			window_process_block[debug_wb_ID - st_wb_ID].print_realign_hap_list_ori(window_info[debug_wb_ID].hap_list.size(), log_f);			fprintf(log_f, "\n");
			for(int l=0; l<4; l++)
				position_depth_counting(l, window_info);
		}
	}

//	if(pre_wb_ID + 1 == st_wb_ID){
//		var_generater.vcf_generatig(pre_wb_pos_depth, pos_depth[0], 	window_info, pre_wb_ID, st_wb_ID + 0, 	pre_window_process_block, 	window_process_block[0]);
//	}else{
//		var_generater.vcf_generatig(pre_wb_pos_depth, BLANK_pos_depth,window_info, pre_wb_ID, 0, 				pre_window_process_block, 	BLANK_window_process_block);
//		var_generater.vcf_generatig(BLANK_pos_depth, pos_depth[0], 	window_info, 0, 		st_wb_ID + 0, 	BLANK_window_process_block, window_process_block[0]);
//	}
//	var_generater.vcf_generatig(pos_depth[0], 	pos_depth[1], window_info, st_wb_ID + 0, st_wb_ID + 1, window_process_block[0], window_process_block[1]);
//	var_generater.vcf_generatig(pos_depth[1], 	pos_depth[2], window_info, st_wb_ID + 1, st_wb_ID + 2, window_process_block[1], window_process_block[2]);//##
//	var_generater.vcf_generatig(pos_depth[2], 	pos_depth[3], window_info, st_wb_ID + 2, st_wb_ID + 3, window_process_block[2], window_process_block[3]);

}

void Window_block_counter_Known_Var::position_depth_counting(uint32_t l, std::vector<Window_t> &window_info){
	haplotype_counter_split_and_realignment(window_process_block[l], window_info);
	Hap_counter &chc = window_process_block[l];
	std::map<uint32_t, Position_depth_item> &cpd = pos_depth[l];
	cpd.clear();
	variant_count.clear();
	candidate_pos.clear();
	uint32_t window_id = chc.get_window_id();
	Window_t &c_window = window_info[window_id];
	//get candidate position
	for(uint32_t hap_id = 0; hap_id < c_window.hap_list.size(); hap_id++){
		for(uint32_t variant_ID : c_window.hap_list[hap_id].variant_list){
			uint32_t var_ref_pos = c_window.var_list[variant_ID].ref_pos;
			int32_t var_len = c_window.var_list[variant_ID].get_length();
			int32_t REF_ALT_SAME_len = c_window.var_list[variant_ID].get_REF_ALT_SAME_len();
//			if(var_len > 0 && REF_ALT_SAME_len < var_len){
//				REF_ALT_SAME_len = 0;
//			}

			uint16_t covarge;
			if(var_len == 0){
				if(var_ref_pos >= 300)	covarge = 0;
				else					covarge = chc.get_depth(hap_id, var_ref_pos);
			}
			else if(var_len < 0){//deletion:
				if(var_ref_pos - var_len + 1 >=300){
					covarge = chc.get_depth(hap_id, var_ref_pos);
				}else{
					uint16_t covarge_left = chc.get_depth(hap_id, var_ref_pos);
					uint16_t right_pos = var_ref_pos - var_len + 1 + REF_ALT_SAME_len;
					uint16_t covarge_right = (right_pos >= 300)?0:chc.get_depth(hap_id, right_pos);
					covarge = MIN(covarge_left, covarge_right);
				}
			}else{
				if(var_ref_pos + 1 >=300){
					covarge = chc.get_depth(hap_id, var_ref_pos);
				}else{
					uint16_t covarge_left = chc.get_depth(hap_id, var_ref_pos);
					uint16_t right_pos = var_ref_pos + 1 + REF_ALT_SAME_len;
					uint16_t covarge_right = (right_pos >= 300)?0:chc.get_depth(hap_id, right_pos);
					covarge = MIN(covarge_left, covarge_right);
				}
			}
			if(covarge == 0)
				continue;
			candidate_pos[var_ref_pos] = REF_ALT_SAME_len;
			variant_count[std::make_tuple(var_ref_pos,hap_id)]=std::make_tuple(variant_ID, covarge);
		}
	}

	//get vcf
	for(auto it = candidate_pos.begin(); it!= candidate_pos.end(); it++){
		uint32_t c_pos = it->first;
		int REF_ALT_SAME_len = it->second;
		Position_depth_item cdi;
		for(uint32_t hap_id = 0; hap_id < c_window.hap_list.size(); hap_id++){
			std::map<std::tuple<uint32_t, uint8_t>, std::tuple<uint32_t, uint32_t>>::iterator it = variant_count.find(std::make_tuple(c_pos, hap_id));
			if(it==variant_count.end()){
				int ref_depth = chc.get_depth(hap_id, c_pos);
				if(REF_ALT_SAME_len != 0){
					uint16_t right_pos = c_pos + REF_ALT_SAME_len;
					uint16_t covarge_right = (right_pos >= 300)?0:chc.get_depth(hap_id, right_pos);
					ref_depth = MIN(ref_depth, covarge_right);
				}
				if(ref_depth != 0)
					cdi.update_ref_supp(window_id, hap_id, ref_depth);
			}
			else{
				uint32_t variant_ID = std::get<0>(it->second);
				uint32_t alt_depth = std::get<1>(it->second);
				cdi.update_variant_supp(window_id, variant_ID, hap_id, alt_depth);
			}
		}
		cpd[c_pos] = (cdi);
	}
}

void Window_block_counter_Known_Var::haplotype_counter_split_and_realignment(Hap_counter &hap_count, std::vector<Window_t> &window_info){
	if(false){
		fprintf(log_f, "post proces window %d \n", hap_count.get_window_id());
	}
	Window_t &window = window_info[hap_count.get_window_id()];
	uint32_t hap_num = window.hap_list.size();
	//resize
	if(hap_count.realign_hap_list.size() < hap_num){
		hap_count.realign_hap_list.resize(hap_num);
		hap_count.realign_hap_list_final.resize(hap_num);
	}

	std::vector<std::vector<uint16_t>> &realign_hap_list = hap_count.realign_hap_list;
	std::vector<std::vector<uint16_t>> &realign_hap_list_return = hap_count.realign_hap_list_final;

	int pre_bias = 0;
	int hap_index = 0;
	int bias = 0;
	//hap_count.print_count_list();

	// copy each counter in global string into each haplotype
	for(uint32_t hap_id = 0;hap_id < window.hap_list.size(); hap_id++){
		//if(hap_id == 33){
		//	fprintf(log_f, " ");
		//}
		hap_t & hap = window.hap_list[hap_id];
		uint32_t hap_length = hap.hap_length;
		realign_hap_list[hap_id].resize(hap_length);
		bias = (hap_length-(8-pre_bias))%8;
		uint32_t cp_to_idx = 0;
		if(pre_bias==0)
			for(; cp_to_idx + 7 < hap_length; cp_to_idx += 8, hap_index++)
				hap_count.cp_from_hap_count_list(&(realign_hap_list[hap_id][cp_to_idx]), hap_index, 0, 8);
		else if(pre_bias != 0){
			hap_count.cp_from_hap_count_list(&(realign_hap_list[hap_id][0]), hap_index, pre_bias, 8);
			cp_to_idx = 8 - pre_bias;
			hap_index += 1;
			for(; cp_to_idx + 7 < hap_length; cp_to_idx += 8, hap_index++)
				hap_count.cp_from_hap_count_list(&(realign_hap_list[hap_id][(cp_to_idx)]), hap_index, 0, 8);
		}
		if(bias!=0)
			hap_count.cp_from_hap_count_list(&(realign_hap_list[hap_id][cp_to_idx]), hap_index, 0, bias);
		pre_bias = bias;
		if(false){
			fprintf(log_f, "part 1 start \n");
			for(uint t=0; t<hap_length;t++)
				fprintf(log_f, "realign i %d : %d \n", t, realign_hap_list[hap_id][t]);
			fprintf(log_f, " part 1 over \n");
		}
	}
	// according to variant, realign count list
	for(uint32_t hap_id = 0;hap_id < window.hap_list.size(); hap_id++){
		hap_t & hap = window.hap_list[hap_id];
		realign_hap_list_return[hap_id].resize(300);
		int alt_st = 0;
		int ref_st = 0;
		uint16_t *cp_to = &(realign_hap_list_return[hap_id][0]);
		if(hap.variant_list.empty() && hap.hap_length != 300){
			memset(cp_to,0,600);
			continue;
		}
		uint16_t *cp_from = &(realign_hap_list[hap_id][0]);
		for(uint32_t var_ID: hap.variant_list){
			Variant &var = window.var_list[var_ID];
			int var_len = var.get_length();
			if(var_len == 0)
				continue;
			int32_t cp_length = var.ref_pos - ref_st + 1;
			if(cp_length < 0)
				continue;
			memcpy(cp_to + ref_st, cp_from + alt_st, (cp_length << 1));
			ref_st += cp_length;
			alt_st += cp_length;
			if(var_len < 0){ //del
				int16_t del_depth = cp_to[ref_st - 1];
				if(del_depth > 0){
					for(int i = 0; i < -var_len; i++){
						cp_to[ref_st + i] = del_depth;
					}
					//memset(cp_to + ref_st, del_depth, ((-var_len)<<1));
					//fprintf(log_f, " ");
				}else{
					memset(cp_to + ref_st, del_depth, ((-var_len)<<1));
				}
				ref_st -= var_len;
			}
			else if(var_len > 0) // ins
				alt_st += var_len;
		}
		if(300 > ref_st)
			memcpy(cp_to + ref_st, cp_from + alt_st, (300 - ref_st)<<1);
		//calling variants
		if(false){
			fprintf(log_f, " realign \n ");
			for(int i=0; i<300;i++)
				fprintf(log_f, " %d : %d\n", i, realign_hap_list_return[hap_id][i]);
			fprintf(log_f, "over \n");
		}
	}
}
