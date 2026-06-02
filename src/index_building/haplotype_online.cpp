/*
 * haplotype.cpp
 *
 *  Created on: 2022年4月14日
 *      Author: fenghe
 */

#include "haplotype_online.hpp"
//#include "SV_handler.hpp"
using namespace ALN_ONLINE;
//return pos
void VCF_Reader::load_line(int & chr_ID, int & pos){
	chr_ID = -1; pos = -1;
	if(is_header_line){
		is_header_line = false;
		//header check:
		fprintf(stderr, "FILTER ITEM: %d" " %s\n", filter_line_number, item_value[filter_line_number].c_str());
		fprintf(stderr, "INFO_line ITEM: %d" " %s\n", INFO_line_idx, item_value[INFO_line_idx].c_str());
		fprintf(stderr, "SAMPLE_BG ITEM: %d" " %s\n", genotype_begin_index, item_value[genotype_begin_index].c_str());
		fprintf(stderr, "ITEM before SAMPLE_BG: %d" " %s\n", genotype_begin_index, item_value[genotype_begin_index - 1].c_str());
		fprintf(stderr, "vcf_ref_line: %d" " %s\n", vcf_ref_line, item_value[vcf_ref_line].c_str());
		fprintf(stderr, "vcf_alt_line: %d" " %s\n", vcf_alt_line, item_value[vcf_alt_line].c_str());
		//check code:
		sample_numer = item_value.size() - genotype_begin_index;
		sample_numer = std::min(sample_numer, max_load_sample);
		fprintf(stderr, "Total sample number: %d\n", sample_numer);
		return;
		//check code: for VENN value
		//AF_analysis_ITEM::get_VENN_value(item_value, true, 0);
	}else{
		//filter line
		if(max_loading_line != -1 && ++total_line_numer >= max_loading_line + 2)
			{reach_EOF = true; return;}
		std::map<std::string, uint32_t>::iterator filter_map_it = filter_map.find(item_value[filter_line_number]);
		if(filter_map_it != filter_map.end())	filter_map_it->second++;
		else filter_map[item_value[filter_line_number]] = 1;
		if(false && item_value[filter_line_number].compare("PASS") != 0) {return;} //NOT RUNNING filter
		chr_ID = get_chr_ID_from_string(item_value[0]);
		pos = atol(item_value[1].c_str());
	}
}

//return true load var number
void VCF_Reader::load_vcf_in_region(int chr_ID, int region_st, int region_ed){
	var_st_idx = 0;
	var_ed_idx = 0;
	//delete VAR from list when POS < region_st
	uint erase_number = 0;
	for(;erase_number < region_GT.size() && (region_GT[erase_number].chrID < chr_ID || region_GT[erase_number].POS < region_st); erase_number++);
	region_GT.erase(region_GT.begin(), region_GT.begin() + erase_number);
	if(!region_GT.empty() && (region_GT[0].chrID > chr_ID || region_GT[0].POS >= region_ed)) return;
	//load new vars into the queue
	while(true){
		if(NULL == gzgets(vcf_file, analysis_line, MAX_LINE_LENGTH))
			{reach_EOF = true; break;}
		else
			split_string(item_value, temp, analysis_line, "\t");
		if(analysis_line[0] == '#' && analysis_line[1] == '#') continue;

		int line_chrID, line_pos;
		load_line(line_chrID, line_pos);

		if(reach_EOF) break;
		if(line_chrID == -1) continue;
		if(line_chrID < chr_ID || line_pos < region_st)
			continue;
		else
			region_GT.emplace_back(item_value, line_chrID, line_pos);

		//fprintf(stderr, "item_value.size() @ %ld, item_value.back() %s \n",item_value.size(), item_value.back().c_str());

		if(line_chrID > chr_ID || line_pos >= region_ed)
			break;
	}
	region_GT_true_size = region_GT.size() + (reach_EOF?0:-1);

	var_st_idx = 0;
	var_ed_idx = region_GT_true_size;
}

void Hap_seq_handler::get_reference(){
	true_region_load_len = 0;
	if(ref_pointer != NULL) free(ref_pointer);
	ref_pointer = NULL;

	int st = load_region_st;
	int ed = (load_region_ed == 0)?0:load_region_ed - 1;
	char reg[1024];
	if(chr_ID < 0){ return;}
	if(chr_ID < 22)			sprintf(reg, "chr%d:%d-%d", chr_ID + 1, st, ed);
	else if(chr_ID == 22) 	sprintf(reg, "chrX:%d-%d", st, ed);
	else if(chr_ID == 23) 	sprintf(reg, "chrY:%d-%d", st, ed);
	else if(chr_ID == 24) 	sprintf(reg, "chrM:%d-%d", st, ed);
	else return;

	ref_pointer = fai_fetch(fai, reg, &true_region_load_len);
}

void Hap_seq_handler::get_haplotype_seq(std::vector<GT_analysis> &region_GT, int h, int region_GT_st, int region_GT_ed){
	var_number = 0;
	HAP_str.clear();
	for(int i = region_GT_st; i < region_GT_ed; i++){
		GT_analysis & region_var = region_GT[i];
		if(region_var.GT[h] == 1){
			var_number++;
			if(var_number == 1){
				HAP_str.clear();
				HAP_str.append(ref_pointer);
				modified_info_init();
				ALT_store_offset = 0;
			}
			//modify the ALT string
			if(region_var.REF.size() == region_var.ALT.size()){//SNV
				int alt_store_pos = region_var.POS - load_region_st + ALT_store_offset;
				if(alt_store_pos < 0 || alt_store_pos >= (long int)HAP_str.size()){
					continue;
				}
				//check:
				if(HAP_str[alt_store_pos] != region_var.REF[0]){
					if(alt_store_pos < (long int)modify_info.size() && modify_info[alt_store_pos].is_modified != true){
						with_wrong_check = true;
						fprintf(stderr, "check error @ h: %d\n", h);
					}
				}
				//change alt string
				for(uint i  = 0; i< region_var.REF.size(); i++){
					if(alt_store_pos + i < HAP_str.size())
						HAP_str[alt_store_pos + i] = region_var.ALT[i];
					if(alt_store_pos - i >= 0 && alt_store_pos + i < (long int)modify_info.size()){
						modify_info[alt_store_pos + i].is_modified =  true;
					}
				}
				//change info string
				for(int i  = 1; i < kmer_size; i++){
					if(alt_store_pos - i >= 0 && alt_store_pos - i < (long int)modify_info.size()){
						modify_info[alt_store_pos - i].is_modified =  true;
					}
				}
			}else{//INDEL
				int alt_store_pos = region_var.POS - load_region_st + ALT_store_offset;
				if(alt_store_pos < 0 || alt_store_pos >= (long int)HAP_str.size()){ continue; }
				ALT_store_offset += region_var.ALT.size() - region_var.REF.size();
				//check:
				if(HAP_str[alt_store_pos] != region_var.REF[0]){
					if(alt_store_pos < (long int)modify_info.size() && modify_info[alt_store_pos].is_modified != true){
						with_wrong_check = true;
						fprintf(stderr, "check error @ h: %d\n", h);
					}
				}
				//erase original ref
				HAP_str.erase(alt_store_pos, region_var.REF.size());
				HAP_str.insert(alt_store_pos, region_var.ALT);
				//change modify string
				if(alt_store_pos < (long int)modify_info.size()){
					modify_info.resize(HAP_str.size());
				}
				//change info string
				for(uint i = 0; i < region_var.ALT.size(); i++){
					if(alt_store_pos + i >= 0 && alt_store_pos + i < (long int)modify_info.size())
						modify_info[alt_store_pos + i].is_modified =  true;
				}
				//WARNING:: only search back kmer_size -1 positions, NOT SAME AS the condition of SNPs
				for(int i  = 1; i < kmer_size - 1; i++){
					if(alt_store_pos - i >= 0 && alt_store_pos - i < (long int)modify_info.size())
						modify_info[alt_store_pos - i].is_modified =  true;
				}
			}
		}
	}
}

void Hap_seq_handler::store_alt_string(std::unordered_set<std::string> & to){
	if(var_number == 0) return;
	//get and store the ALT string
	uint modify_info_size = modify_info.size();
	for(uint i = 0; i < modify_info_size; i++){
		if(modify_info[i].is_modified == true){
			uint kmer_str_st = i;
			for(;i < modify_info_size; i++){
				if(modify_info[i].is_modified !=  true)
					break;
			}
			std::string kmer_str;
			//basic info
			kmer_str.append(std::to_string(chr_ID));
			kmer_str.append("_");
			kmer_str.append(std::to_string(modify_info[kmer_str_st].original_pos));
			kmer_str.append("_");
			kmer_str.append(std::to_string(modify_info[i - 1].original_pos));
			kmer_str.append("_");
			//alt string
			kmer_str.append(HAP_str.begin() + kmer_str_st, HAP_str.begin() + i + kmer_size - 1);
			//start and end info
			to.insert(kmer_str);
		}
	}
}

void HAP_string_single_VCF_builder::store_match(std::vector<uint8_t> &hap_seq,
		int length) {
	while (length > 127) {
		hap_seq.emplace_back(0xff);
		length -= 127;
	}
	if (length > 0)
		hap_seq.emplace_back((length << 1) + 1);
}

void HAP_string_single_VCF_builder::store_SNP(std::vector<uint8_t> &hap_seq,
		char alt_char) {
	switch (alt_char) {
	case 'A':
	case 'a':
		hap_seq.emplace_back(60 << 2); //60 << 2
		break;
	case 'C':
	case 'c':
		hap_seq.emplace_back(61 << 2);
		break;
	case 'G':
	case 'g':
		hap_seq.emplace_back(62 << 2);
		break;
	case 'T':
	case 't':
		hap_seq.emplace_back(63 << 2);
		break;
	default:
		hap_seq.emplace_back(59 << 2);
		break; //59 means error
	}
}

void HAP_string_single_VCF_builder::store_del(std::vector<uint8_t> &hap_seq,
		int length) {
	if (length > 0 && length < 50)
		hap_seq.emplace_back((length << 2) + 0);
	if (length >= 50) {
		//only handle SNP and INDELS,when it is SV, donothing
		return;
		hap_seq.emplace_back((63 << 2) + 2); //0xfe
		xassert(length < 0xffffff, "SV too long");
		hap_seq.emplace_back(length >> 16);
		hap_seq.emplace_back((length >> 8) & 0xffff);
		hap_seq.emplace_back((length) & 0xff);
	}
}

void HAP_string_single_VCF_builder::store_ins(std::vector<uint8_t> &hap_seq,
		int length, std::string &alt) {
	if (length > 0 && length < 50) {
		hap_seq.emplace_back((length << 2) + 2);
		const char *ALT_str = alt.c_str() + 1;
		for (int i = 0; i < length; i++)
			hap_seq.emplace_back(ALT_str[i]);
	}
	if (length >= 50) {
		//only handle SNP and INDELS,when it is SV, do nothing
		return;
		hap_seq.emplace_back((62 << 2) + 2); //0xfa
		xassert(length < 0xffffff, "SV too long");
		hap_seq.emplace_back(length >> 16);
		hap_seq.emplace_back((length >> 8) & 0xffff);
		hap_seq.emplace_back((length) & 0xff);
		const char *ALT_str = alt.c_str() + 1;
		for (int i = 0; i < length; i++)
			hap_seq.emplace_back(ALT_str[i]);
	}
}

/**
 * input: a vcf reads, and all variance in a window
 * output:
 * (1) an minimizer list for alt-string (todo::)
 * (2) a compact data struct for hap-seq
 * * ***/
std::vector<uint8_t>& HAP_string_single_VCF_builder::build_alt_string_idx(
		VCF_Reader &vcf_reader, Hap_seq_handler &hh, int load_region_st_,
		int load_region_ed_, int step_length_, int kmer_size_) {
	//final:
	std::vector<uint8_t> &window_block = window_block_buff;
	window_block.clear();
	//block count, 16 bit, occupied by 0 at first
	window_block.emplace_back(0);
	window_block.emplace_back(0);
	if (vcf_reader.is_empty()) {
		return window_block;
	}
	//S1: loading reference
	hh.set_ref(vcf_reader.get_chr_ID(), load_region_st_, load_region_ed_,
			step_length_, kmer_size_);
	//S2.0 init
	unique_alt_string_set.clear(); //unused by now
	unique_hap_seq_set.clear();
	hap_seq_list.clear();
	//store ref string to hap_seq_list
	ref_hap_seq_string.clear();
	ref_hap_seq_string.append(hh.ref_pointer);
	unique_hap_seq_set.insert(ref_hap_seq_string);
	//S2: generate alt string for each haplotype
	int GT_size = vcf_reader.get_sample_number() * 2;
	for (int h = 0; h < GT_size; h++) {
		//generate haplotype seq for this hap
		hh.get_haplotype_seq(vcf_reader.region_GT, h, vcf_reader.var_st_idx,
				vcf_reader.var_ed_idx);
		if (hh.var_number == 0)
			continue; //reference
		//store hap_string
		bool is_new_hap = unique_hap_seq_set.insert(hh.HAP_str).second;
		if (is_new_hap) {
			//if(load_region_st_ == 75511350)	fprintf(stderr, "new @ %d \n", h);
			//store in compact struct
			//store new haplotype
			hap_seq_list.emplace_back();
			std::vector<uint8_t> &hap_seq = hap_seq_list.back();
			//if(load_region_st_ == 75511350)	fprintf(stderr, "hap_seq_list.size() @ %ld \n", hap_seq_list.size());
			int old_ref_pos = load_region_st_;
			std::vector<GT_analysis> &region_GT = vcf_reader.region_GT;
			for (int i = vcf_reader.var_st_idx; i < vcf_reader.var_ed_idx; i++){
				GT_analysis &region_var = region_GT[i];
				if (region_var.GT[h] == 1) {
					if( region_var.POS < load_region_st_ || region_var.POS >= load_region_ed_)
						continue;
					int offset_change = region_var.ALT.size()
							- region_var.REF.size();
					if(offset_change <= -50 || offset_change >= 50){
						//skip SVs DO NOTHING
					}
					else if (offset_change == 0) {				//SNV
						int match_len = region_var.POS - old_ref_pos;
						xassert(match_len < 300, "hap match length over 300 bp\n");//debug::todo::
						store_match(hap_seq, match_len);
						old_ref_pos = region_var.POS + 1;
						//store SNPs
						store_SNP(hap_seq, region_var.ALT[0]);
					} else {				//DEL + INS
						int match_len = region_var.POS - old_ref_pos + 1;
						store_match(hap_seq, match_len);
						if (offset_change < 0) {				//DEL
							old_ref_pos = region_var.POS - offset_change + 1;
							store_del(hap_seq, -offset_change);
						} else {				//INS
							old_ref_pos = region_var.POS + 1;
							//store INS:
							store_ins(hap_seq, offset_change, region_var.ALT);
						}
					}
				}
			}
			//store the final part
			int match_len = load_region_ed_ - old_ref_pos;
			store_match(hap_seq, match_len);
		}
	}
	//if(load_region_st_ == 75511350)	fprintf(stderr, "ALL hap_seq_list.size() @ %ld \n", hap_seq_list.size());
	//seq list
	for (auto &hap_seq_bin : hap_seq_list) {
		uint32_t hap_seq_size = hap_seq_bin.size();
		xassert(hap_seq_size < 0x7fffff, "");
		if (hap_seq_size <= 127)
			window_block.emplace_back(hap_seq_size << 1);
		else {
			window_block.emplace_back(((hap_seq_size >> 16) << 1) + 1);
			window_block.emplace_back((hap_seq_size >> 8) & 0xffff);
			window_block.emplace_back((hap_seq_size) & 0xff);
		}
		//if(load_region_st_ == 75511350)	fprintf(stderr, "new @ %ld \n", window_block.size());
	}
	uint32_t block_count = window_block.size() - 2;
	xassert(block_count < 0xffff, "");
	//if(load_region_st_ == 75511350)	fprintf(stderr, "HEAD \n");
	window_block[0] = (block_count >> 8);
	window_block[1] = (block_count & 0xff);
	//store seq detail:
	for (auto &hap_seq_bin : hap_seq_list){
		window_block.insert(window_block.end(), hap_seq_bin.begin(),
					hap_seq_bin.end());
		//if(load_region_st_ == 75511350)	fprintf(stderr, "new @ %ld \n", window_block.size());
	}
	//fprintf(stderr, "ALL @ %ld \n", window_block.size());
	return window_block;
}

void HAP_string_single_VCF_builder::dump_final_window_block(char *wb_dump_fn) {
	if (wb_dump_fn != NULL) {
		FILE *dump_f = xopen(wb_dump_fn, "wb");	//equal to : fopen(wb_dump_fn, wb_dump_fn);
		CPP_vector_dump_bin(dump_f, &(final_window_block_info[0]),
				final_window_block_info.size() * sizeof(window_block_info));
		CPP_vector_dump_bin(dump_f, &(final_window_block_buff[0]),
				final_window_block_buff.size() * sizeof(uint8_t));
		fclose(dump_f);
	}
	fprintf(stderr, "final_window_block_info.size : %ld \t final_window_block_buff.size %ld \n", final_window_block_info.size(), final_window_block_buff.size());
}

void Minimizer_index_builder::minimizer_generater_hap_seq(std::vector<Hap_seq_modify_info> &modify_info, std::vector<char> &hap_seq,
		uint64_t window_ID, uint64_t hap_ID, std::vector<mm128_t> &mm_v, bool is_from_SV){
	xassert(hap_seq.size() <= 0xfffff, ""); // at most 1M
	xassert(window_ID <= 0x1ffffff, ""); // at most 32M
	xassert(hap_ID <= 0x7ffff, ""); // at most 0.5M
	uint modify_info_size = modify_info.size();
	for(uint i = 0; i < modify_info_size; i++){
		if(modify_info[i].is_modified == true){
			uint alt_str_st = i;
			for(;i < modify_info_size; i++){
				if(modify_info[i].is_modified !=  true)
					break;
			}
			//basic info
			alt_string.clear();
			alt_string.append(hap_seq.begin() + alt_str_st, hap_seq.begin() + i + kmer_size_alt_string - 1);
			//basic check:
			bool used_in_this_window = false;
			if(is_from_SV){
				used_in_this_window = true;
			}else{
				if(modify_info[alt_str_st].original_pos >= 75 && modify_info[alt_str_st].original_pos < 225){//when in sweet region
					used_in_this_window = true;
				}
				else if(modify_info[alt_str_st].original_pos < 75){//in the too early region
					alt_string_with_pos.clear();
					alt_string_with_pos.append(std::to_string(modify_info[alt_str_st].original_pos + 150));
					alt_string_with_pos.append("_");
					alt_string_with_pos.append(alt_string);
					if(alt_string_uniq_set_old.find(alt_string_with_pos) == alt_string_uniq_set_old.end()){//search in the alt_string_uniq_set of previous window block
						used_in_this_window = true;
					}
				}else if(modify_info[alt_str_st].original_pos >= 225){//in the too old region
					if(i == modify_info_size - 1){//reach the end of the window block, keep the alt string
						used_in_this_window = true;
					}
				}
			}
			if(!used_in_this_window)
				continue;
			alt_string_with_pos.clear();
			alt_string_with_pos.append(std::to_string(modify_info[alt_str_st].original_pos));
			alt_string_with_pos.append("_");
			alt_string_with_pos.append(alt_string);
			//alt string
			//start and end info
			if(alt_string_uniq_set.insert(alt_string_with_pos).second){// true insert:
				if(false)//debug code
					fprintf(stderr, "Alt string：%s\n", alt_string_with_pos.c_str());
				//generate minimizer
				uint64_t mm_st_idx = mm_v.size();
				mg.mm_sketch(&(alt_string[0]), alt_string.size(), mm_v);
				//modified the minimizer
				//From original:
//						X：[56bit hash value] + [8bit(useless)]
//						Y：[32bit rid] + [31bit offset] + [1 bit Reverse or Forward]
				//new modified:
				// X：[56bit hash value] + [7bit(useless)] + [1 bit Reverse or Forward]
				// Y：[25bit window ID @ MAX 32M space] + [19bit: hap ID @ MAX 0.5M space] + [20bit: hap offset @ MAX 1M]
				// set window id/ hap ID and offset and direction for minimizer
				uint64_t mm_ed_idx = mm_v.size();
				for(uint32_t mm_i = mm_st_idx; mm_i < mm_ed_idx; mm_i++){
					mm128_t & c_mm = mm_v[mm_i];
					c_mm.x &= 0xffffffffffffff00;
					c_mm.x += c_mm.y & 0x1;
					uint32_t hap_offset = ((c_mm.y & 0xffffffff) >> 1) + alt_str_st;
					xassert(hap_offset <= 0xfffff, ""); // at most 1M
					c_mm.y = (window_ID << 39) + (hap_ID << 20) + (hap_offset);
				}
			}
		}
	}
}

void Minimizer_index_builder::minimizer_generater_unitig(std::string &UNITIG_str, uint64_t unitig_ID, std::vector<mm128_t> &mm_v){
		xassert(UNITIG_str.size() <= 0xffffffff, ""); // at most 4G
		xassert(unitig_ID <= 0xffffffff, ""); // at most 4G
		uint64_t mm_st_idx = mm_v.size();
		mg.mm_sketch(&(UNITIG_str[0]), UNITIG_str.size(), mm_v);
		uint64_t mm_ed_idx = mm_v.size();
		//From original:
		//X：[56bit hash value] + [8bit(useless)]
		//Y：[32bit rid] + [31bit offset] + [1 bit Reverse or Forward]
		//new modified:
		// X：[56bit hash value] + [7bit(useless)] + [1 bit Reverse or Forward]
		// Y：[25bit window ID @ MAX 32M space] + [19bit: hap ID @ MAX 0.5M space] + [20bit: hap offset @ MAX 1M]
		// set window id/ hap ID and offset and direction for minimizer
		for(uint32_t mm_i = mm_st_idx; mm_i < mm_ed_idx; mm_i++){
			mm128_t & c_mm = mm_v[mm_i];
			c_mm.x &= 0xffffffffffffff00;
			c_mm.x += c_mm.y & 0x1;
			uint32_t hap_offset = ((c_mm.y & 0xffffffff) >> 1);
			c_mm.y = (unitig_ID << 32) + (hap_offset);
		}
	}

	void Minimizer_index_builder::building_mm_idx_and_dump(FILE * mm_v_f, FILE * mm_idx_f, std::vector<mm128_t> &sort_mm_v, int kmer_size_minimizer){
		uint32_t max_mm_buff_size = 0xffff;//64K
		//open the mm_v dump file
		xassert(kmer_size_minimizer > 13 && kmer_size_minimizer <= 28, "");
		xassert(sort_mm_v.size() < 0xffffffff, " Number of Minimizer is over 4G");
		int index_length = 26;//26 bit
		int minimizer_length = kmer_size_minimizer * 2; //at most 56 bit
		int mini_left_length = minimizer_length - index_length; // at most 30 bp
		int mini_left_mask = (0x1 << mini_left_length) - 1;//0x3fffffff when left == 30; and 0x3ffff when left = 18;
		mm128_t * mmv_bg = &(sort_mm_v[0]);
		mm128_t * mmv_ed = &(sort_mm_v[0]) + sort_mm_v.size();
		std::vector<mm_96_t> mm_96_buff;//
		//init index:
		std::vector<uint32_t> mm_96_idx;
		uint32_t index_size = (0x1 << index_length) + 1;
		mm_96_idx.resize(index_size); //64m * 4 = 256M
		for(uint32_t i = 0; i < index_size; i++){
			mm_96_idx[i] = MAX_uint32_t;
		}

		for(mm128_t * c_mm = mmv_bg; c_mm < mmv_ed; c_mm++){
			//handler X part
			uint64_t hash_key = (c_mm->x >> 8);
			uint64_t idx_bucket = (hash_key >> mini_left_length);
			uint64_t left_key = (hash_key & mini_left_mask);
			uint64_t is_forward = (c_mm->x & 0x1);
			uint32_t X = (left_key << 2) + is_forward;
			uint32_t Y = (c_mm->y >> 32);
			uint32_t Z = (c_mm->y & 0xffffffff);
			//store the mm_96
			mm_96_buff.emplace_back(X, Y, Z);
			if(mm_96_buff.size() == max_mm_buff_size || c_mm == mmv_ed - 1){
				fwrite(&(mm_96_buff[0]), sizeof(mm_96_t),mm_96_buff.size(), mm_v_f);
				mm_96_buff.clear();
			}

			//store index
			if(mm_96_idx[idx_bucket] == MAX_uint32_t){
				mm_96_idx[idx_bucket] = c_mm - mmv_bg;
			}
		}
		//store the final one:
		mm_96_idx[index_size - 1] = mmv_ed - mmv_bg;
		//check all -1 value
		for(int32_t i = index_size - 1; i >= 0; i--){
			if(mm_96_idx[i] == MAX_uint32_t){
				mm_96_idx[i] = mm_96_idx[i + 1];
			}
		}
		//write the index file
		fwrite(&(mm_96_idx[0]), sizeof(uint32_t),mm_96_idx.size(), mm_idx_f);
	}

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8)
//
//#define sort_key_64(x) (x)
//KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

void Minimizer_index_builder::mm_sort(std::vector<mm128_t> &mm_v){
	radix_sort_128x(&(mm_v[0]), &(mm_v[0]) + mm_v.size());
}

#define MIN_sweet_region 0
#define MAX_sweet_region 300

std::vector<char> & Window_block_handler::unfold_hap_seq(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, bool generate_modify_string, bool generate_var_list, int kmer_size_alt_string){
	//get hap seq:
	std::vector<char> & rst_haplotype = rst_haplotype_BUFF;
	rst_haplotype.clear();
	//modify string init:
	if(generate_modify_string){
		kmer_size_alt_string = MAX(kmer_size_alt_string, 22);
		modify_info.clear();
		int modify_info_size = ref.size();
		for(int i = 0; i < modify_info_size; i++){
			modify_info.emplace_back(i);
		}
	}
	int32_t c_ref_pos = 0;
	int32_t c_alt_pos = 0;
	if(generate_var_list){
		var_l.clear();
	}

	//generate hap strings
	for(uint32_t i = 0; i < hap_seq_len;){
		if((bin_hap_seq[i] & 1) == 1){//match
			int match_len = 0;
			for(;i < hap_seq_len && (bin_hap_seq[i] & 1) == 1; i++){ match_len += (bin_hap_seq[i] >> 1); }
			//debug code:
			match_len = (c_ref_pos + match_len > (int)ref.size())?((int)ref.size() - c_ref_pos):match_len;
			int32_t max_load_pos = std::min(c_ref_pos + match_len, (int)ref.size());
			if(c_ref_pos < max_load_pos)
				rst_haplotype.insert(rst_haplotype.end(), ref.begin() + c_ref_pos, ref.begin() + max_load_pos);
			c_ref_pos += match_len;
			c_alt_pos += match_len;
		}
		else if((bin_hap_seq[i] & 2) == 0){//SNP or DEL
			uint8_t real_data = (bin_hap_seq[i] >> 2);
			if(real_data <= 63 && real_data >= 60){//SNP
				char SNV_char = "ACGT"[(real_data) - 60];
				i++;
				rst_haplotype.emplace_back(SNV_char);
				if(generate_modify_string){
					for(int mod_kmer_i = 0; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
						if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
							modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
						}
					}
				}
				if(generate_var_list){
					xassert(c_alt_pos >= 0, "");
					std::string REF; REF.insert(REF.begin(),ref.begin() + c_ref_pos, ref.begin() + c_ref_pos + 1);
					std::string ALT; ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + 1);
					var_l.emplace_back(c_ref_pos, REF, ALT);
					if(c_ref_pos < MIN_sweet_region || c_ref_pos >= MAX_sweet_region ){
						rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size());
						rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin(), var_l.back().REF.end());
						var_l.erase(var_l.end() - 1);
					}
				}
				c_ref_pos += 1;
				c_alt_pos += 1;
			}else{//short del
				int32_t del_length = real_data;
				xassert(del_length > 0 && del_length < 50,"");
				i++;
				if(generate_modify_string){
					for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
						if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
							modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
						}
					}
				}
				if(generate_var_list){
					std::string REF;  std::string ALT;
					if(c_ref_pos + del_length >= ref.size())
						 del_length = ref.size() - c_ref_pos;
					if(del_length >= 0){
						if(c_ref_pos >= 1){
							REF.insert(REF.begin(),ref.begin() + c_ref_pos - 1, ref.begin() + c_ref_pos + del_length);
							ALT.insert(ALT.begin(),ref.begin() + c_ref_pos - 1, ref.begin() + c_ref_pos);
						}else{//deletion reach the edge of windows block
							ALT = "N";
							REF = "N" + REF;
							REF.insert(REF.begin(),ref.begin() + c_ref_pos, ref.begin() + c_ref_pos + del_length);
						}
						var_l.emplace_back(c_ref_pos - 1, REF, ALT);
						if(c_ref_pos < MIN_sweet_region || c_ref_pos - 1 >= MAX_sweet_region ){
							//rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size());
							rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin() + 1, var_l.back().REF.end());
							var_l.erase(var_l.end() - 1);
							//c_ref_pos -= del_length;
							c_alt_pos += del_length;//
						}
					}
				}
				c_ref_pos += del_length;
			}
		}else{//short INS or SV
			uint8_t real_data = (bin_hap_seq[i] >> 2);
			if((real_data) == 63){//long DEL
				uint32_t del_length = ((bin_hap_seq[i+1]) << 16) + (bin_hap_seq[i + 2] << 8) + (bin_hap_seq[i+3]);
				i+= 4;
				//skip the long DEL in this function
				if(false){
					if(generate_modify_string){
						for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
							if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
							}
						}
					}
					c_ref_pos += del_length;//todo:: consider the condition when the deletion is too long
				}
			}else{//INS
				bool is_long_ins = ((real_data) == 62);
				uint32_t ins_length = 0;
				if(is_long_ins){
					ins_length = ((bin_hap_seq[i+1]) << 16) + (bin_hap_seq[i + 2] << 8) + (bin_hap_seq[i+3]);
					i += 4 + ins_length;
				}
				else{
					ins_length = real_data;
					i += 1 + ins_length;
				}
				//skip the long INS in this function
				if(!is_long_ins){
					uint64_t be_ins_pos = rst_haplotype.size();
					rst_haplotype.resize(rst_haplotype.size() + ins_length);
					for(uint32_t i_char_idx = 0; i_char_idx < ins_length; i_char_idx++){
						rst_haplotype[be_ins_pos + i_char_idx] = bin_hap_seq[i - ins_length + i_char_idx];
					}
					if(generate_modify_string){
						modify_info.resize(rst_haplotype.size());
						for(uint32_t mod_kmer_i = 0; mod_kmer_i < ins_length; mod_kmer_i++){
							if(c_alt_pos + mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos + mod_kmer_i].is_modified =  true;
							}
						}
						for(int mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
							if(c_alt_pos - mod_kmer_i >= 0 && c_alt_pos - mod_kmer_i < (long int)modify_info.size()){
								modify_info[c_alt_pos - mod_kmer_i].is_modified =  true;
							}
						}
					}

					if(generate_var_list){
						std::string REF;  std::string ALT;
						if(c_ref_pos >= 1){
							REF.insert(REF.begin(),rst_haplotype.begin() + c_alt_pos - 1, rst_haplotype.begin() + c_alt_pos);
							ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos - 1, rst_haplotype.begin() + c_alt_pos + ins_length);
						}else{//deletion reach the edge of windows block
							REF = "N";
							ALT = "N" + ALT;
							ALT.insert(ALT.begin(),rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + ins_length);
						}
						var_l.emplace_back(c_ref_pos - 1, REF, ALT);

						if(c_ref_pos < MIN_sweet_region || c_ref_pos - 1 >= MAX_sweet_region ){
							rst_haplotype.erase(rst_haplotype.begin() + c_alt_pos, rst_haplotype.begin() + c_alt_pos + var_l.back().ALT.size() - 1);
							//rst_haplotype.insert(rst_haplotype.begin() + c_alt_pos, var_l.back().REF.begin(), var_l.back().REF.end());
							var_l.erase(var_l.end() - 1);
							c_alt_pos -= ins_length;
						}
					}
					c_alt_pos += ins_length;
				}
			}
		}
	}
	if(generate_modify_string){
		modify_info.resize(rst_haplotype.size());
	}
	return rst_haplotype;
}

std::vector<char> & Window_block_handler::unfold_hap_seq_SV(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, window_block_info & c_wb, bool generate_modify_string, uint32_t kmer_size_alt_string){
	//get hap seq:
	std::vector<char> & rst_haplotype = rst_haplotype_BUFF;
	rst_haplotype.clear();
/*
	int bin_hap_se_idx = 0;
	//int sv_type = SV_TYPE::flag_2_type(bin_hap_seq[bin_hap_se_idx]); bin_hap_se_idx += 1;//SV type
	bin_hap_se_idx += 2;//

	int sv_length = SV_STORE_BASIC::get_32bit(bin_hap_seq + bin_hap_se_idx); bin_hap_se_idx += 4;//SV length
	bin_hap_se_idx += 6;//pos of BP 1
	bin_hap_se_idx += 6;//pos of BP 2
	uint32_t hap_seq_length = SV_STORE_BASIC::get_32bit(bin_hap_seq + bin_hap_se_idx);

	bin_hap_se_idx += 4;//hap seq_length
	rst_haplotype.resize(hap_seq_length + 1);
	for(uint32_t i = 0; i < hap_seq_length; i++){
		rst_haplotype[i] = (char)bin_hap_seq[bin_hap_se_idx + i];
	}
	rst_haplotype[hap_seq_length] = 0;
	bin_hap_se_idx += hap_seq_length;

	if(generate_modify_string){
		kmer_size_alt_string = MAX(kmer_size_alt_string, 22);
		modify_info.clear();
		uint32_t modify_info_size = hap_seq_length;
		for(uint32_t i = 0; i < modify_info_size; i++){
			modify_info.emplace_back(i);
		}
		//reset
		bin_hap_se_idx = 0;
		bin_hap_se_idx += 1;//SV type
		uint32_t alt_base_bg = SV_STORE_BASIC::get_16bit(bin_hap_seq + bin_hap_se_idx);	bin_hap_se_idx += 2;//SV type
		bin_hap_se_idx += 4;//SV length

		uint16_t chr_ID1, chr_ID2;
		uint32_t POS1, POS2;
		SV_STORE_BASIC::get_chrID_POS_pair(bin_hap_seq + bin_hap_se_idx, chr_ID1, POS1); bin_hap_se_idx += 6;//pos of BP 1
		SV_STORE_BASIC::get_chrID_POS_pair(bin_hap_seq + bin_hap_se_idx, chr_ID2, POS2); bin_hap_se_idx += 6;//pos of BP 1
		if(sv_type == SV_TYPE::INS || sv_type == SV_TYPE::DEL){
			uint32_t modified_base_bg = alt_base_bg;
			uint32_t modified_base_ed = alt_base_bg + ((sv_type == SV_TYPE::DEL)?(1):(sv_length + 1));

			//fprintf(stderr, "modify_info.size() %d\n", modify_info.size());
			for(uint32_t mod_kmer_i = modified_base_bg; mod_kmer_i < modified_base_ed; mod_kmer_i++){
				if(mod_kmer_i < modify_info.size()){
					//fprintf(stderr, "%d\n", mod_kmer_i);
					modify_info[mod_kmer_i].is_modified =  true;
				}
			}
			for(uint32_t mod_kmer_i = 1; mod_kmer_i < kmer_size_alt_string; mod_kmer_i++){
				if(alt_base_bg >= mod_kmer_i && alt_base_bg - mod_kmer_i < modify_info.size()){
					//fprintf(stderr, "%d\n", alt_base_bg - mod_kmer_i);
					modify_info[alt_base_bg - mod_kmer_i].is_modified =  true;
				}
			}
		}
	}

*/
	return rst_haplotype;
}

void Window_block_handler::unfold_window_block_all(uint8_t * window_block_p){
	if(window_block_p_old == window_block_p)	return;
	else										window_block_p_old = window_block_p;

	block_all_hap_list.clear();
	uint32_t block_count = (window_block_p[0] << 8) + window_block_p[1];
	int c_hap_id = 0;
	uint8_t * hap_seq_load_p = window_block_p + 2;
	uint64_t hap_offset = 0;
	uint64_t hap_seq_len = 0;
	//get hap_seq_len and hap_offset
	while(true){
		if(hap_seq_load_p - window_block_p >= block_count + 2){ hap_offset = -1; break; }//out of boundary
		if(((hap_seq_load_p[0]) & 0x1) == 0){ 	hap_seq_len = ( hap_seq_load_p[0] >> 1); 															hap_seq_load_p ++;	} //one byte
		else{									hap_seq_len = ((hap_seq_load_p[0] >> 1) << 16) + (hap_seq_load_p[1] << 8) + (hap_seq_load_p[2]);	hap_seq_load_p += 3;} //three bytes
		block_all_hap_list.emplace_back();
		block_all_hap_list.back().p = window_block_p + 2 + block_count + hap_offset;
		block_all_hap_list.back().len = hap_seq_len;
		hap_offset += hap_seq_len;
		c_hap_id++;
	}
}

void Window_block_handler::unfold_window_block_i(uint8_t * window_block_p, uint64_t hap_i){
	block_i_rst_address.clear();
	if(hap_i < 0) {return;}
	uint32_t block_count = (window_block_p[0] << 8) + window_block_p[1];
	int c_hap_id = 0;
	uint8_t * hap_seq_load_p = window_block_p + 2;
	uint64_t hap_offset = 0;
	uint64_t hap_seq_len = 0;
	//get hap_seq_len and hap_offset
	while(true){
		if(hap_seq_load_p - window_block_p > block_count + 2){ hap_offset = -1; break; }//out of boundary
		if(((hap_seq_load_p[0]) & 0x1) == 0){ 	hap_seq_len = ( hap_seq_load_p[0] >> 1); 															hap_seq_load_p ++;	} //one byte
		else{									hap_seq_len = ((hap_seq_load_p[0] >> 1) << 16) + (hap_seq_load_p[1] << 8) + (hap_seq_load_p[2]);	hap_seq_load_p += 3;} //three bytes
		if(hap_i == c_hap_id) break;//find, break
		hap_offset += hap_seq_len;
		c_hap_id++;
	}
	if(hap_offset == -1){ return; }
	block_i_rst_address.p = window_block_p + 2 + block_count + hap_offset;
	block_i_rst_address.len = hap_seq_len;
}

#define BASE_PER_LINE 70
int hapseq2fa(int argc, char *argv[]){
	MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
	std::string INDEX_PATH = argv[optind];
	//std::string INDEX_PATH = "/media/fenghe/data2/4_20/index_dir";
	idx->load_all_index(INDEX_PATH.c_str(), NULL, true);
    //idx->load_all_index("/media/fenghe/data2/4_20/index_dir", NULL, true);
	hap_string_loader_single_thread hl_r1;
	Simple_ref_handler ref;
	//ref1.load_bin_ref("/media/fenghe/data2/4_20/index_dir");
	ref.load_bin_ref(INDEX_PATH.c_str());
	std::string ref_str;
    //count the wb data
	int new_line_count = 0;
	//---------------------------------------------------------------------------------------------
	//-------------------------------------PART1 HAPLOTYPES----------------------------------------
	//---------------------------------------------------------------------------------------------
	int64_t total_hap_length = 0; int64_t total_150_read_size = 0;
    for(uint64_t window_ID = 0; window_ID < idx->wb_info_size; window_ID++){
//    	if(idx->wb_info[i].flag != 1)
//    		continue;
    	ref.load_ref_from_buff(idx->wb_info[window_ID].chrID, idx->wb_info[window_ID].region_st, idx->wb_info[window_ID].region_length,  ref_str);
    	//std::vector<std::string> & hap_seq_v = hl_r1.get_string_list(window_ID, ref_str, idx);
    	String_list_and_var_list & hap_seq_v = hl_r1.get_string_list_and_var_list(window_ID, ref_str, idx);
		if(new_line_count % BASE_PER_LINE != 0)
			fprintf(stdout, "\n");
    	fprintf(stdout, ">window_%ld HAP_SEQ\n", window_ID);
		new_line_count = 0;
		if(hap_seq_v.hap_string_l.empty()){
			fprintf(stdout, "NNNNNNNNNNNNNNNN\n");
		}

		//check: get all kmers
		std::unordered_set<std::string> st;
		if(false){
			int total_length = 0;
			st.clear();
			total_length += ref_str.size();
			for(uint i = 0; i  + 150 < ref_str.size(); i++){
				std::string I;
				I.insert(I.begin(), ref_str.begin() + i,ref_str.begin() + i + 150);
				st.emplace(I);
			}

			for(std::string & h:hap_seq_v.hap_string_l){
				total_length += h.size();
				for(uint i = 0; i + 150 < h.size(); i++){
					std::string I;
					I.insert(I.begin(), h.begin() + i,h.begin() + i + 150);
					st.emplace(I);
				}
			}
			total_hap_length += (total_length -150); total_150_read_size += (st.size());
			fprintf(stderr, "Useful BWT rate: %f total_length %d ",(float)st.size()/total_length ,total_length);
			fprintf(stderr, "Useful BWT rate: %f total_length %ld @ window_ID %ld\n",(float)total_150_read_size/total_hap_length ,total_hap_length,window_ID);
		}

		for(auto & hap_str: hap_seq_v.hap_string_l){
			if(hap_str.empty()){
				fprintf(stderr, "BLANK string: %s\n", hap_str.c_str());
				idx->wb_info[window_ID].show_wb(stderr);
			}
//			for(uint i = 0; i < hap_str.size(); i++){
//				if(hap_str[i] != 'A' && hap_str[i] != 'C' && hap_str[i] != 'G'&& hap_str[i] != 'T'&& hap_str[i] != 'N'){
//					fprintf(stderr, "Error base at %d, %ld %s\n",i, global_count, hap_str.c_str());
//					hap_str[i] = 'N';
//				}
//			}
			for(uint i = 0; i < hap_str.size() ; i++){
				fprintf(stdout, "%c",hap_str[i]);
				new_line_count++;
				if(new_line_count % BASE_PER_LINE == 0){
					fprintf(stdout, "\n");
				}
			}
			//fprintf(stdout, "%s\n", S.c_str());
		}
    }
	if(new_line_count % BASE_PER_LINE != 0)
		fprintf(stdout, "\n");

	//---------------------------------------------------------------------------------------------
	//-------------------------------------PART2 REFERENCE-----------------------------------------
	//---------------------------------------------------------------------------------------------
	uint32_t chr_number = idx->wb_info[idx->wb_info_size - 1].chrID + 1;
	chr_number = MIN(ref.get_chr_N(), chr_number);
	for(uint32_t chr_id = 0; chr_id < chr_number; chr_id++){
		uint32_t chr_len = ref.get_chr_length(chr_id);
		fprintf(stdout, ">%s len_%d\n", ref.get_chr_name(chr_id).c_str(), chr_len);
		ref.load_ref_from_buff(chr_id, 1, chr_len,  ref_str);
		//check G string:
		uint G_string_size = 0;
		for(uint i = 0; i < chr_len + 1; i++){
			if(i < chr_len && ref_str[i] == 'G')
				G_string_size ++;
			else{
				if(G_string_size > 300)
					for(uint j = 0; j < G_string_size ; j++)
						ref_str[i - 1 - j] = 'N';
				G_string_size = 0;
			}
		}
		//output:
		new_line_count = 0;
		for(uint i = 0; i < chr_len ; i++){
			fprintf(stdout, "%c",ref_str[i]);
			new_line_count++;
			if(new_line_count % BASE_PER_LINE == 0){
				fprintf(stdout, "\n");
			}
		}
		if(new_line_count % BASE_PER_LINE != 0)
			fprintf(stdout, "\n");
	}
	return 0;
}

//int print_var_list(int argc, char *argv[]){
//	std::string INDEX_PATH = argv[optind];
//	MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
//    idx->load_all_index(INDEX_PATH.c_str(), NULL, true);
//
//    Simple_ref_handler ref1;
//	ref1.load_bin_ref(INDEX_PATH.c_str());
//	hap_string_loader_single_thread hl_r1;
//	std::string ref;
//	for(uint i = 0; i < idx->wb_info_size; i++){
//		window_block_info &c_wb_info = idx->wb_info[i];
//		if(c_wb_info.is_SV()){ break; }
//		fprintf(stderr, "ID %d\n", i);
//		c_wb_info.show_wb(stdout);
//		ref1.load_ref_from_buff(c_wb_info.chrID, c_wb_info.region_st, c_wb_info.region_length, ref);
//		String_list_and_var_list & A = hl_r1.get_string_list_and_var_list(i, ref, (idx));
//		A.print(stdout);
//	}
//	return 0;
//}

struct Repeat_region{ int st; int ed; Repeat_region(int st_, int ed_):st(st_),ed(ed_){} };
struct Simple_repeat_checker{
private:
	bool already_clear = true;
	bool with_data = false;
	uint8_t* repeat_region;
	void clear_repeat_region(uint repeat_region_size){
		memset(&(repeat_region[0]), 0, repeat_region_size);
		already_clear = true;
	}
public:
	std::vector<Repeat_region> r;
	void init(int repeat_region_size){
		repeat_region = (uint8_t*)xcalloc(repeat_region_size + 100,1);
		with_data = false;
		already_clear = true;
	}

	void simple_repeat_check(std::string &ref){
		r.clear();
		//return; //this function is abandoned
		const char * s = ref.c_str();
		uint max_check_size = 25;
		uint max_repeat_len = 0;
		for(uint i = 1; i < max_check_size + 1; i++){
			for(uint j = i; j < ref.size(); j++){
				if(s[j] == s[j - i] && s[j] != 'N'){
					max_repeat_len++;
				}else{
					if(max_repeat_len > 15){
						uint region_ed = j;
						uint region_bg = (region_ed-1)-((max_repeat_len) - 1);
						with_data = true;
						for(uint k = region_bg; k < region_ed; k++)
							repeat_region[k] = 1;
						if(false)
							fprintf(stderr, "max_repeat_len %d, i %d , [%d, %d) \n", max_repeat_len, i, region_bg, region_ed);
					}
					max_repeat_len = 0;
				}
			}
			if(max_repeat_len > 15){
				uint region_ed = ref.size();
				uint region_bg = (region_ed-1)-((max_repeat_len) - 1);
				with_data = true;
				for(uint k = region_bg; k < region_ed; k++)
					repeat_region[k] = 1;
				if(false)
					fprintf(stderr, "max_repeat_len %d, i %d ,[%d, %d) \n", max_repeat_len, i, region_bg, region_ed);
			}
			max_repeat_len = 0;
		}
		//finally:
		if(with_data){
			//get data
			if(false)
				fprintf(stderr, "REF：%s\n", ref.c_str());
			uint region_bg = 0;
			uint region_ed = 0;
			for(uint j = 0; j < ref.size(); j++){
				if(repeat_region[j] == 1){
					region_ed = j+1;
				}else{
					if(region_ed > region_bg){
						r.emplace_back(region_bg, region_ed);
						if(false)
							fprintf(stderr, "load region [%d, %d) \n", region_bg, region_ed);
					}
					region_bg = region_ed = j+1;
				}
			}
			//final one
			if(region_ed > region_bg){
				r.emplace_back(region_bg, region_ed);
				if(false)
					fprintf(stderr, "load region [%d, %d) \n", region_bg, region_ed);
			}
			clear_repeat_region(ref.size());
			with_data = false;
		}
	}
};

void remove_duplicated_VARs(std::vector<HAP_VAR_ITEM> & var_v){
	for(uint i = 1; i < var_v.size(); i++){
		//when both vars are in SAME position, try to combine
		if(var_v[i - 1].ref_pos == var_v[i].ref_pos){
			if(var_v[i - 1].try_to_combine(var_v[i])){
				//successfully combine
				var_v.erase(var_v.begin() + i);
				i--;
				if(var_v[i].is_no_var()){
					var_v.erase(var_v.begin() + i);
				}
			}
		}
		//when both first var is DEL and the second var is INDEL, and they are in similar position, try to combine
		else if(var_v[i - 1].get_len() < 0 && var_v[i - 1].ref_pos - var_v[i - 1].get_len() >= var_v[i].ref_pos && var_v[i].get_len() != 0 && var_v[i - 1].get_len() != 0 ){
			if(var_v[i - 1].try_to_combine(var_v[i])){
				//successfully combine
				var_v.erase(var_v.begin() + i);
				i--;
				if(var_v[i].is_no_var()){
					var_v.erase(var_v.begin() + i);
				}
			}
		}
	}
}

int get_var_ALT_REF_same_length(HAP_VAR_ITEM & var_v, std::string & REF){
//	std::string varREF = "CAAAAAAAAA";
//	std::string varALT = "C";
//	HAP_VAR_ITEM var_v(46, varREF, varALT);
//	std::string REF;
//	REF = "CTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGA";

//	fprintf(stderr, "%s \n", REF.c_str());
//	fprintf(stderr, "%s ", var_v.REF.c_str());
//	fprintf(stderr, "%s ", var_v.ALT.c_str());
//	fprintf(stderr, "\n");

	int SAM_len = 0;
	//complex vars
	if(var_v.ALT.size() > 1 && var_v.REF.size() > 1){
		std::string ALT = REF;
		ALT.erase(ALT.begin() + var_v.ref_pos, ALT.begin() + var_v.ref_pos + var_v.REF.size());
		ALT.insert(ALT.begin() + var_v.ref_pos, var_v.ALT.begin(), var_v.ALT.end());
	//	fprintf(stderr, "%s\n" ,REF.c_str());
	//	fprintf(stderr, "%s\n" ,ALT.c_str());
		int max_search = std::min(ALT.size(), REF.size());
		for(int i = var_v.ref_pos + 1; i < max_search; i++){
			if(REF[i] == ALT[i])	SAM_len++;
			else					break;
		}
	}else{//simple indels
		int var_len = var_v.get_len();
		if(var_len < 0){//del
			int max_search = REF.size() + var_len;
			for(int i = var_v.ref_pos + 1; i < max_search; i++){
				if(REF[i] == REF[i - var_len])	SAM_len++;
				else							break;
			}
		}else if( var_len > 0){ //ins
			bool is_end = false;
			for(int i = 0; i < var_len; i++){
				if(REF[var_v.ref_pos + i + 1] == var_v.ALT[i + 1])		SAM_len++;
				else{													is_end = true;	break;	}
			}
			if(!is_end){
				int max_search = REF.size() - var_len;
				for(int i = var_v.ref_pos + 1; i < max_search; i++){
					if(REF[i] == REF[i + var_len])						SAM_len++;
					else												break;
				}
			}
		}
	}
	return SAM_len;
}

void VAR_FORMAT(HAP_VAR_ITEM & v, std::string & hap_str, int INDEL_offset, int MIN_offset){
	if(v.get_len() == 0) return; //skip SNPs
	if(v.REF.size() !=1 && v.ALT.size() != 1) return;//skip complex vars
	if(v.get_len() > 0){
		while((int)v.ref_pos > MIN_offset && (int)v.ref_pos + INDEL_offset > 0 && v.ALT.back() == hap_str[v.ref_pos + INDEL_offset] && v.REF.back() == hap_str[v.ref_pos + INDEL_offset]){
			v.ref_pos --;
			v.ALT = hap_str[v.ref_pos + INDEL_offset] + v.ALT;
			v.ALT.pop_back();
			v.REF = hap_str[v.ref_pos + INDEL_offset];
		}
	}else{//DEL
		while((int)v.ref_pos > MIN_offset && (int)v.ref_pos + INDEL_offset > 0 && v.ALT.back() == hap_str[v.ref_pos + INDEL_offset] && v.REF.back() == hap_str[v.ref_pos + INDEL_offset]){
			v.ref_pos --;
			v.REF = hap_str[v.ref_pos + INDEL_offset] + v.REF;
			v.REF.pop_back();
			v.ALT = hap_str[v.ref_pos + INDEL_offset];
		}
	}
}

void VAR_FORMAT_HAP(std::vector<HAP_VAR_ITEM> &var_l, std::string & hap_str){
	int INDEL_offset = 0; int MIN_offset = 0;
	for(HAP_VAR_ITEM & v: var_l){
		VAR_FORMAT(v, hap_str, INDEL_offset, MIN_offset);
		INDEL_offset += v.get_len();
		MIN_offset = v.ref_pos;
	}
}

int print_var_list_compect(int argc, char *argv[]){
	std::string INDEX_PATH = argv[optind];
	MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
    idx->load_all_index(INDEX_PATH.c_str(), NULL, true);

    Simple_ref_handler ref1;
	ref1.load_bin_ref(INDEX_PATH.c_str());
	hap_string_loader_single_thread hl_r1;
	hap_string_loader_single_thread hl_DEBUG;
	std::string ref;
	FILE * outf = stdout;
	fprintf(stdout, "Total_wb_size %ld\n", idx->wb_info_size);
	char *var_string = (char *)xcalloc(1000000, 1);
	std::map<std::string, uint32_t> var_map;

	Simple_repeat_checker src;
	src.init(10000);

	for(uint wb_id = 0; wb_id < idx->wb_info_size; wb_id++){
		if(248821 == wb_id)
			fprintf(stderr, " ");
		//else	continue;

		window_block_info &c_wb_info = idx->wb_info[wb_id];
		if(c_wb_info.is_SV()){ break; }

		ref1.load_ref_from_buff(c_wb_info.chrID, c_wb_info.region_st, c_wb_info.region_length, ref);
		String_list_and_var_list & A = hl_r1.get_string_list_and_var_list(wb_id, ref, (idx));
		//if(A.hap_string_l.empty()){
		//	continue;
		//}
		src.simple_repeat_check(ref);

		var_map.clear();
		//basic
		fprintf(outf, "%d %d %d %ld ", wb_id, c_wb_info.chrID, c_wb_info.region_st, A.hap_string_l.size());
		//repeat regions
		fprintf(outf, "%ld ", src.r.size());
		for(auto & r: src.r)
			fprintf(outf, "%d %d ", r.st, r.ed);
		fprintf(outf, "\n");

		//output haplotype
		uint32_t var_ID = 0;
		for(uint hap_ID = 0; hap_ID < A.hap_string_l.size(); hap_ID++){
			if(82 == hap_ID)
				fprintf(stderr, " ");
			fprintf(outf, "%ld ", A.hap_string_l[hap_ID].size());
			//remove duplicated VARs
			remove_duplicated_VARs(A.var_l[hap_ID]);
			//debug: PRE check for all haplotype:
			if(true){
				int var_len = 0;
				for(auto & v: A.var_l[hap_ID]){
					int c_var_len = v.get_len();
					if(ABS(c_var_len) >= 50 ){
						A.var_l[hap_ID].clear();
						break;
					}
					var_len += v.get_len();
				}
				if(!A.var_l[hap_ID].empty()){
					if(A.hap_string_l[hap_ID].size() - var_len != ref.size()){
						//SV check:
						int diff_len = A.hap_string_l[hap_ID].size() - var_len - ref.size();
						if(diff_len <= -50 || diff_len >= 50)
							A.var_l[hap_ID].clear();
						else
							A.var_l[hap_ID].clear();
					}
				}
				//hl_DEBUG.get_string_list_and_var_list(wb_id, ref, (idx));
			}

			VAR_FORMAT_HAP(A.var_l[hap_ID], A.hap_string_l[hap_ID]);
			remove_duplicated_VARs(A.var_l[hap_ID]);
			//store haplotypes
			for(HAP_VAR_ITEM & v: A.var_l[hap_ID]){
				if(v.ALT.back() == 'N')
					fprintf(stderr, " ");
				if(v.get_len() != 0){ //INDELs
					int SAME_length = get_var_ALT_REF_same_length(v, ref);
					if(SAME_length != 0)
						sprintf(var_string,  "%s_%s_%d_%d", v.REF.c_str(), v.ALT.c_str(), v.ref_pos, SAME_length);
					else
						sprintf(var_string,  "%s_%s_%d ", v.REF.c_str(), v.ALT.c_str(), v.ref_pos);
				}else //SNPs
					sprintf(var_string,  "%s_%s_%d ", v.REF.c_str(), v.ALT.c_str(), v.ref_pos);
				std::map<std::string, uint32_t>::iterator it = var_map.find(var_string);
				int32_t c_var_ID = -1;
				if(it == var_map.end()){
					c_var_ID = var_ID;
					var_map[var_string] = var_ID++;
				}else{
					c_var_ID = it->second;
				}
				fprintf(outf, "%d ", c_var_ID);
			}
			fprintf(outf, "\n");
		}
		//output vars
		fprintf(outf, "# %ld\n", var_map.size());
		for(auto & var: var_map){
			fprintf(outf, "@ %s %d\n", var.first.c_str(), var.second);
		}
	}
	return 0;
}
