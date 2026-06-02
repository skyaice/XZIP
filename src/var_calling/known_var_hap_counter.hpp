/*
 * hap_count.h
 *
 *  Created on: 2022年6月3日
 *      Author: zyx
 */

#ifndef HAP_COUNT_H_
#define HAP_COUNT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <vector>
extern "C"{
	#include "../clib/utils.h"
}
#include <cstring>

#define ADD_BASE 0x0101010101010101
#define HELP_BASE 0xffffffffffffffff
#define HIGH_DIGIT_CHECK  0x8080808080808080
#define ADD_BASE_TAIL 0x0000000000000001
#define ADD_BASE_HEAD 0x0100000000000000
#define CHECK_TAIL 0x0000000000000080
#define CHECK_HEAD 0x8000000000000000
//#define READ_LENGTH 150

class Hap_counter{

	uint32_t window_id;


	std::vector<uint64_t> save_list;
	uint32_t total_hap_len;

	void add_helper(int head_tail, uint64_t pos, uint64_t *add_pos, uint64_t *save_pos);

public:

	std::vector<uint64_t> hap_count_list;

	std::vector<std::vector<uint16_t>> realign_hap_list;
	std::vector<std::vector<uint16_t>> realign_hap_list_final;

	//get the reads list that support the haplotype
	void get_hap_read_list(uint32_t hap_ID, std::vector<int> &read_pos_list){
		read_pos_list.clear();
		std::vector<uint16_t> &c_hap = realign_hap_list[hap_ID];
		uint32_t hap_len = c_hap.size();
		for(int k = 0; k < c_hap[0]; k++)
			read_pos_list.emplace_back(0);
		for(uint j = 1; j < hap_len; j++){
			int read_at_pos = c_hap[j] - c_hap[j - 1];
			for(int k = 0; k < read_at_pos; k++)
				read_pos_list.emplace_back(j);
		}
	}

	bool is_hap_with_read(uint32_t hap_ID){
		std::vector<uint16_t> &c_hap = realign_hap_list[hap_ID];
		uint32_t hap_len = c_hap.size();
		for(uint j = 0; j < hap_len; j++){
			if(c_hap[j] != 0)
				return true;
		}
		return false;
	}

	void print_realign_hap_list_final(uint32_t hap_num, FILE * log_f){
		fprintf(log_f, "T:\t");
		for(uint i = 0;i< 300; i++){
			fprintf(log_f, "%d\t", i);
		}
		fprintf(log_f, "\n");
		for(uint i = 0;i< hap_num; i++){
			fprintf(log_f, "%d:\t", i);
			for(uint j:realign_hap_list_final[i]){
				fprintf(log_f, "%d\t", j);
			}
			fprintf(log_f, "\n");
		}
	}

	void print_realign_hap_list_ori(uint32_t hap_num, FILE * log_f){
		fprintf(log_f, "T:\t");
		for(uint i = 0;i< 300; i++){
			fprintf(log_f, "%d\t", i);
		}
		fprintf(log_f, "\n");
		for(uint i = 0;i< hap_num; i++){
			fprintf(log_f, "%d:\t", i);
			for(uint j:realign_hap_list[i]){
				fprintf(log_f, "%d\t", j);
			}
			fprintf(log_f, "\n");
		}
		std::vector<int> read_pos_list;

		for(uint i = 0;i< hap_num; i++){
			fprintf(log_f, "%d:\t", i);
			get_hap_read_list(i, read_pos_list);
			for(int j: read_pos_list)
				fprintf(log_f, "%d\t", j);
			fprintf(log_f, "\n");
		}
	}

public:

	inline void cp_from_hap_count_list(uint16_t * realign_hap_list, uint32_t hap_index, int from, int to){
		uint64_t source_lower_7bit = hap_count_list[hap_index];
		uint64_t source_upper_8bit = save_list[hap_index];
		for(int n=from; n < to; n++){
			realign_hap_list[n - from] =    (source_lower_7bit>>((7-n)*8))&0xff;
			realign_hap_list[n - from] += (((source_upper_8bit>>((7-n)*8))&0xff) << 7);
		}
	}

	uint32_t get_window_id(){ return window_id; }

	FILE* log_f;

	void init(uint32_t total_hap_len_, uint32_t window_id_, FILE* log_f_){
		total_hap_len = total_hap_len_;
		window_id = window_id_;
		uint64_t int64_bit_len = ((total_hap_len_ + 8)/8);
		if(int64_bit_len > hap_count_list.size()){
			hap_count_list.resize(int64_bit_len);
			save_list.resize(int64_bit_len);
		}
		memset(&(hap_count_list[0]),0,(int64_bit_len << 3));
		memset(&(save_list[0]),0,(int64_bit_len << 3));
		log_f = log_f_;
	}

	inline uint16_t get_depth(uint32_t hap_ID, uint32_t POS){
		return realign_hap_list_final[hap_ID][POS];
	}

	void add_signal(uint64_t hap_pos, uint32_t read_length){
		//reset for the reference alignment
		if(hap_pos + read_length > total_hap_len){
			fprintf(log_f, "OVER BOUND REF @  window_id %d, total_hap_len %d hap_pos %ld\n", window_id, total_hap_len, hap_pos);
				return;
		}
		uint64_t head_hap_pos = hap_pos % 8;
		uint64_t *add_pos = &hap_count_list[0] + hap_pos/8;
		uint64_t *save_pos = &save_list[0] + hap_pos/8;

		if(head_hap_pos!=0){
			add_helper(1, head_hap_pos, add_pos, save_pos);
			add_pos = add_pos + 1;
			save_pos = save_pos + 1;
		}
		int count = (read_length-(8-head_hap_pos))/8;
		if(head_hap_pos == 0)
			count ++;
		for( int i=0; i<count; i++){
			add_helper(1, 0, add_pos, save_pos);
			add_pos = add_pos + 1;
			save_pos = save_pos + 1;
		}
		uint64_t  last= (read_length-(8-head_hap_pos))%8;
		if(last!=0){
			add_helper(0, last, add_pos,save_pos);
		}
	}
};

#endif /* HAP_COUNT_H_ */
