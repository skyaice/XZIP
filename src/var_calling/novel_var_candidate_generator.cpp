/*
 * novel_var_caller.cpp
 *
 *  Created on: 2023年2月12日
 *      Author: fenghe
 */

#include "novel_var_candidate_generator.hpp"

namespace Aln_online{

void Single_WB_signal_handler_Novel::clear(){
	gap_g.clear();
	gap_um.clear();
	wb_change_detail.clear();
	read_mapping_position.clear();
	is_mapping_position_sorted = false;
	already_assembly = false;
}

void Single_WB_signal_handler_Novel::sort_mapping_position(){
	if(is_mapping_position_sorted) return;
	std::sort(read_mapping_position.begin(), read_mapping_position.end());
	is_mapping_position_sorted = true;
}

//2bit: type
//4bit: length
//10 bit: data
#define TYPE_SNP_change_detail_32 0
#define TYPE_INS_change_detail_32 1
#define TYPE_DEL_change_detail_32 2
#define GET_TYPE_change_detail_32(c32) (c32 >> 30)
#define GET_INDEL_LEN_change_detail_32(c32) ((c32 >> 26) & 0xf)
//#define GET_INDEL_LEN_change_detail_32(c32) (c32 >> 28) & 0x
#define REMOVE_EDGE_READ 3

void Single_WB_signal_handler_Novel::change_detail_parse(std::vector<uint16_t> &ori, std::vector<uint32_t> &rst, int r_read_len){
	rst.clear();
	//INS check:
	int indel_adjust = 0;
	for(uint16_t idx = 0; idx < ori.size(); idx++){
		int change_position = ori[idx] >> 3;
		if(change_position >= r_read_len - REMOVE_EDGE_READ){//remove the change detail near the edge of reads
			continue;
		}
		if((ori[idx] & 0x7) > 4){ //SNP
		//	if(indel_adjust == 0)
		//		rst.emplace_back(ori[idx]);//store directory
		//	else if(indel_adjust > 0){
			if(change_position < REMOVE_EDGE_READ)//remove the (SNP)change detail near(3bp) the edge(begin) of reads
				continue;
			rst.emplace_back(ori[idx] + (indel_adjust << 3));//store directory
//			}else{
//				rst.emplace_back(ori[idx] - ((-indel_adjust) << 3));//store directory
//			}
		}
		else if((ori[idx] & 0x7) == 4){//DEL
			uint16_t del_pos = ori[idx] >> 3;
			uint16_t del_length = 1;
			//search forward
			for(uint j = idx + 1; j < ori.size(); j++){
				if((ori[j] & 0x7) == 4 && (ori[j] >> 3) == del_pos)//DEL in the same pos
					del_length++;
			}
			//store the DEL results
			uint32_t del_rst = (0x2 << 30);//type : 2 bit
			del_rst += del_length << 26; //length; 4bit
			del_rst += (((del_pos+indel_adjust) << 3)); //position: 13 bit
			rst.emplace_back(del_rst);
			//skip parsed data
			idx += (del_length - 1);
			indel_adjust += del_length;
		}else{//INS
			uint16_t ins_pos = ori[idx] >> 3;
			uint16_t ins_length = 1;
			std::vector<uint8_t> ins_string;
			ins_string.emplace_back(ori[idx] & 0x7);
			//search forward
			for(uint j = idx + 1; j < ori.size(); j++){
				if((ori[j] & 0x7) < 4 && (ori[j] >> 3) == (ins_pos + ins_length)){//INS in the nearby
					ins_length++;
					ins_string.emplace_back(ori[j] & 0x7);
				}
			}
			//store the INS results
			uint32_t ins_rst = (0x1 << 30);//type : 2 bit
			ins_rst += ins_length << 26; //length; 4bit
			int ins_string_store_pos = 24;
			for(uint8_t c: ins_string){//INS string: 8bit
				ins_rst += (c & 0x3) << (ins_string_store_pos);
				ins_string_store_pos -= 2;
			}
			ins_rst += ((ins_pos+indel_adjust) << 3); //position: 13 bit
			rst.emplace_back(ins_rst);
			//finally
			idx += (ins_length - 1);
			indel_adjust -= ins_length;
		}
	}
	//remove fake position in the end of DEL reads
	//todo::
}

void Single_WB_signal_handler_Novel::add_signal(Single_read_aln_rst_unpack & r){
	//statistics
	if(r.is_CPX()){
		gap_um.emplace_back();
		gap_um.back().store(r);
	}
	else{
		//P0
		//reset read length, read length is at least the position of change
//		for(uint16_t change: r.change_detail)
//			if(((change >> 3) + 1) > r.read_len){
//				r.read_len = ((change >> 3) + 1);
//				//xassert(r.read_len < 150, "");
//			}
		//r.generate_read_string(r.wb_ID, r.read_len);
		//r.print_detail(stderr);
		//P1: store change detail:
		change_detail_parse(r.change_detail, parse_change_detail_buff, r.read_len);
		for(uint32_t i = 0; i < parse_change_detail_buff.size(); i++){
			//convert the position in read into the position in WB; the lower 3 bits is error type, other bits is position
			uint32_t change = parse_change_detail_buff[i];
			//convert the position in read into the position in WB; the lower 3 bits is error type, other bits is position
			change += (r.read_hamming_mapping_position << 3);
			//store the change detail
			//debug code if(change > 10000){	fprintf(stderr, " "); }
			std::map<uint32_t, uint32_t>::iterator it = wb_change_detail.find(change);
			if(it != wb_change_detail.end())
				it->second ++;
			else
				wb_change_detail[change] = 1;
		}

		//P2: store read start and end information, it is used for calculating depth information
		int read_pos_st = r.read_hamming_mapping_position + REMOVE_EDGE_READ;
		int read_pos_ed = r.read_hamming_mapping_position + r.read_len - REMOVE_EDGE_READ;
		//fprintf(stderr, "read_pos_st %d, read_pos_ed %d \n", read_pos_st, read_pos_ed);
		read_mapping_position.emplace_back(read_pos_st << 1);
		read_mapping_position.emplace_back((read_pos_ed << 1) + 1);
		//P3: copy data to list:
		gap_g.emplace_back(); gap_g.back().store(r);
	}
}

void Nearby_WB_Signal_handler_Novel::init(ALN_ONLINE::MM_idx_loader *wb_idx, ALN_ONLINE::Simple_ref_handler *ref_h, FILE*vcf_out_log, bool show_assembly){
	this->show_assembly = show_assembly;
	this->log_f = vcf_out_log;
	cur_4.resize(4);
	old_4.resize(4);

	result_4.resize(4);
	for(int i = 0; i < 4; i++){
		result_4[i].init(wb_idx, ref_h);//todo::
	}

	//statistics
	total_wb_number = 0;
	total_ASS_number = 0;
}

void Nearby_WB_Signal_handler_Novel::clear(uint32_t wb_basic){
	for(int i = 0; i < 4; i++){
		cur_4[i].clear();
		cur_4[i].wb_id = wb_basic + i;
		cur_4[i].wb_basic = wb_basic;
		result_4[i].clear();
	}
}

//calling variants in for nearby WB
void Nearby_WB_Signal_handler_Novel::generate_novel_vars_candidate_in_4_nearby_WB(){
	for(int i = 0; i < 4; i++){
		//for the final one of OLD and the first one of new
		Single_WB_signal_handler_Novel * pre = (i == 0) ? &(old_4[3]):&(cur_4[i - 1]);
		Single_WB_signal_handler_Novel * next = &(cur_4[i]);
		total_ASS_number += result_4[i].call_novel_vars_single(pre, next, log_f, show_assembly);
	}
	total_wb_number += 4;
}

uint8_t charToDna5n[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,//4 the last but one
    /*   		 A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//'Z'
    /*          	      T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*   		 a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*          	      t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig){
	int contig_seq_len = contig_string.size();
	const char * contig_seq = contig_string.c_str();
	xassert(nullptr != contig_seq, "");
	//store bin contig
	bin_contig.resize(contig_seq_len);
	for (int i = 0; i < contig_seq_len; ++i)
		bin_contig[i] = charToDna5n[(uint8_t)contig_seq[i]];
}

}


