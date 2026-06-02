/*
 * occ_RST_DEF.cpp
 *
 *  Created on: 2022年12月10日
 *      Author: fenghe
 */

#include "num_def.hpp"
#include "occ_RST_DEF.hpp"
#include<vector>

//bit-level functions
//PART1: basic value: for cluster store type and BWT grow index
int calculate_gap_aln_store_type(uint64_t cluster_number)
{
	if(cluster_number == 0)
		return GAP_ALN_STORE_TYPE_NULL;
	if(cluster_number == 1)
		return GAP_ALN_STORE_TYPE_UNIQ;
	else if(cluster_number <= GAP_ALN_STORE__SMALL_MAX__)//<= 10 right-now
		return GAP_ALN_STORE_TYPE_SMALL;
	else
		return GAP_ALN_STORE_TYPE_BIG;
}

//store in upper 1~18 bit in K
void set_gap_bwt_GROW_idx(
	//input:
	uint32_t gap_bwt_GROW_idx,
	//output:
	bwtint_t &new_k){
	//clear:
	// 1000 0000 0000 0000 0001 1111 1111 1111 1111 1111
	// 0x80001fffff
	new_k &= 0x80001fffff;
	//set the data:
	new_k += (((uint64_t) gap_bwt_GROW_idx) << 21);
}

//store in upper 25~26 bit in K
void set_gap_aln_store_type(
	//input:
	uint32_t gap_aln_store_type,
	//output:
	bwtint_t &new_k){
	//clear:
	// 1111 1111 1111 1111 1111 1111 1001 1111 1111 1111 //40 bit:
	// 0xffffff9fff
	new_k &= 0xffffff9fff;
	//set the data:
	new_k += ((gap_aln_store_type) << 13);
}

//that is "seed length"
//store in upper 19~24 bit in K
void set_gap_bwt_grow_length(
	//input:
	uint8_t gap_bwt_grow_length,
	//output:
	bwtint_t &new_k){
	//clear:
	// 1111 1111 1111 1111 1110 0000 0111 1111 1111 1111
	// 0xffffe07fff
	new_k &= 0xffffe07fff;
	//set the data:
	new_k += (((uint64_t)gap_bwt_grow_length) << 15);
}

void get_gap_aln_metadata(
		//input:
		bwtint_t old_k,
		//output:
		uint32_t &gap_bwt_GROW_idx, int32_t &gap_aln_store_type, uint8_t &gap_bwt_grow_length
		){
	gap_bwt_GROW_idx = (old_k >> 21) & 0x3ffff;
	gap_aln_store_type = (old_k >> 13) & 0x3;
	gap_bwt_grow_length = (old_k >> 15) & 0x3f;
}

int32_t get_gap_aln_store_type_FOR_PAIRING(	bwtint_t old_k, bwtint_t old_l){
	return (old_k >> 13) & 0x3;
}

//PART2: for uniq clusters, store directly in the same place
void set_uniq_cluster_data(
		//input:
		Gap_aln_Cluster_t &c,
		//output:
		bwtint_t &new_k, bwtint_t &new_l){
	//clear
	// 1111 1111 1111 1111 1111 1111 1110 0000 0000 0000
	new_k &= 0xffffffe000;
	// 1000 0000 0000 0000 0000 0000 0000 0000 0000 0000
	new_l &= 0x8000000000;
	//set the data
	//score
	new_k += (((uint64_t)c.score) << 3);
	//WB ID：
	new_l += (((uint64_t)c.WB_ID) << 13);
	//read_to_wb_offset
	new_l += (((uint64_t)c.read_to_wb_offset) << 3);
	//is_rev
	new_l += (((uint64_t)c.is_rev) << 2);
}

//not need to clear data list
void get_uniq_cluster_data(
		//input:
		bwtint_t old_k, bwtint_t old_l,
		//output:
		std::vector<Gap_aln_Cluster_t> & c_v){
	//score
	uint16_t score = (old_k >> 3) & 0x3ff;
	uint32_t WB_ID = (old_l >> 13) & 0x3ffffff;
	uint8_t is_rev = (old_l >> 2) & 0x1;
	uint16_t read_to_wb_offset = (old_l >> 3) & 0x3ff;
	c_v.emplace_back(WB_ID, is_rev, score, read_to_wb_offset);
}

void get_uniq_cluster_data_FOR_PAIRING(
		//input:
		bwtint_t old_k, bwtint_t old_l,
		//output:
		int32_t &wb_id, int8_t &is_rev, int16_t &score, uint16_t &read_to_wb_offset){
	score = (old_k >> 3) & 0x3ff;
	wb_id = (old_l >> 13) & 0x3ffffff;
	is_rev = (old_l >> 2) & 0x1;
	read_to_wb_offset = (old_l >> 3) & 0x3ff;
}

//PART3: for small clusters
void set_small_cluster_meta(
		//input:
		uint32_t cluster_store_bg_idx,
		uint32_t cluster_store_ed_idx,
		//output:
		bwtint_t &new_l){
	//0x0111 1111 1111 1111 1111
	//0x7ffff
	if(cluster_store_ed_idx > 0x7ffff ){
		fprintf(stderr, "%d ", cluster_store_ed_idx);
		xassert(cluster_store_ed_idx <= 0x7ffff, "FATAL ERROR: ERR RATE in reads is too high for handling!\n At most store 512k aln clusters in total for each block\n");
	}
	
	//clear
	// 1000 0000 0000 0000 0000 0000 0000 0000 0000 0000
	new_l &= 0x8000000000;
	//set the data
	new_l += (((uint64_t)cluster_store_bg_idx) << 20);
	new_l += (((uint64_t)cluster_store_ed_idx) << 1);
}

void get_small_cluster_meta(
		//input:
		bwtint_t old_l,
		//output:
		uint32_t &cluster_store_bg_idx,
		uint32_t &cluster_store_ed_idx){
	cluster_store_bg_idx = (old_l >> 20) & 0x7ffff;
	cluster_store_ed_idx = (old_l >> 1) & 0x7ffff;
}

//PART4: for BIG clusters
void set_big_cluster_meta(
		//input:
		uint64_t __SSD__FILE__OFFSET__,
		uint32_t __SSD__FILE__SIZE__,
		//output
		bwtint_t &old_k,
		bwtint_t &old_l
){
	//clear
	// 1111 1111 1111 1111 1111 1111 1110 0000 0000 0000
	old_k &= 0xffffffe000;
	// 1000 0000 0000 0000 0000 0000 0000 0000 0000 0000
	old_l &= 0x8000000000;
	//set:
	__SSD__FILE__SIZE__ /= sizeof(Gap_aln_Cluster_t);
	xassert(__SSD__FILE__SIZE__ <= 0xfff, "At most store 4096 aln clusters");
	old_k += ((uint64_t)__SSD__FILE__SIZE__ << 1);
	old_k += ((uint64_t)__SSD__FILE__OFFSET__ >> 39);
	old_l += ((uint64_t)__SSD__FILE__OFFSET__) & 0x7fffffffff;
}

void get_big_cluster_meta(
		//input:
		bwtint_t new_k, bwtint_t new_l,
		//output:
		uint64_t &__SSD__FILE__OFFSET__,
		uint32_t &__SSD__FILE__SIZE__
		){
	__SSD__FILE__SIZE__ = (new_k >> 1) & 0xfff;
	__SSD__FILE__OFFSET__ = (new_k << 39) & 0xffffffffff;
	__SSD__FILE__OFFSET__ += ((new_l) & 0x7fffffffff);
	__SSD__FILE__SIZE__ *= sizeof(Gap_aln_Cluster_t);
}

void set_big_cluster_meta_FOR_SA(uint32_t big_cluster_store_bg_idx, uint32_t big_cluster_store_ed_idx, uint64_t &k, uint64_t &l){
	//setK
	// #							  ##
	// 1XXX XXXX XXXX XXXX XXXX XXXX X11X XXXX XXXX XXXX
	k = 0x8000006000;
	k += (big_cluster_store_bg_idx & 0x1fff);//13 bit
	k += (((uint64_t)big_cluster_store_bg_idx << 2) & 0xffff8000);//other 17 bit

	//setL:
	// #
	// 0XXX XXXX XXXX XXXX XXXX XXXX XXXX XXXX XXXX XXXX
	l = 0x0;
	l += big_cluster_store_ed_idx;
}

void get_big_cluster_meta_FOR_SA(uint32_t &big_cluster_store_bg_idx, uint32_t &big_cluster_store_ed_idx, uint64_t k, uint64_t l){
	//load BG:
	big_cluster_store_bg_idx = (k & 0x1fff);
	big_cluster_store_bg_idx += ((k>>2) & 0xffffe000);
	//load ED:
	big_cluster_store_ed_idx = (l & 0xffffffff);
}


