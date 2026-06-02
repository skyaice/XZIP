/*
 * aln_RST_DEF.cpp
 *
 *  Created on: 2023年2月12日
 *      Author: fenghe
 */

#include "aln_RST_DEF.hpp"
#include "index_building/haplotype_online.hpp"

namespace Aln_online {

void set_bit_48_little_endian(uint64_t store, uint8_t *r){
	r[5] = (store >> 40) & 0xff;;
	r[4] = (store >> 32) & 0xff;;
	r[3] = (store >> 24) & 0xff;;
	r[2] = (store >> 16) & 0xff;;
	r[1] = (store >> 8) & 0xff;;
	r[0] = (store) & 0xff;;
}

void load_bit_48_little_endian(uint64_t &store, uint8_t *r){
	store = 0;
	store += (r[5]); store <<= 8;
	store += (r[4]); store <<= 8;
	store += (r[3]); store <<= 8;
	store += (r[2]); store <<= 8;
	store += (r[1]); store <<= 8;
	store += (r[0]);
}

void set_bit_40_little_endian(uint64_t store, uint8_t *r){
	r[4] = store >> 32;
	r[3] = (store >> 24) & 0xff;
	r[2] = (store >> 16) & 0xff;
	r[1] = (store >> 8) & 0xff;
	r[0] = store & 0xff;
}

void load_bit_40_little_endian(uint64_t &store, uint8_t *r){
	store = 0;
	store += (r[4]); store <<= 8;
	store += (r[3]); store <<= 8;
	store += (r[2]); store <<= 8;
	store += (r[1]); store <<= 8;
	store += (r[0]);
}

void set_bit_32_little_endian(uint64_t store, uint8_t *r){
	r[3] = (store >> 24) & 0xff;
	r[2] = (store >> 16) & 0xff;
	r[1] = (store >> 8)  & 0xff;
	r[0] = store & 0xff;
}

void load_bit_32_little_endian(uint64_t &store, uint8_t *r){
	store = 0;
	store += (r[3]); store <<= 8;
	store += (r[2]); store <<= 8;
	store += (r[1]); store <<= 8;
	store += (r[0]);
}

void Single_read_aln_rst_t1_48_exact::set(int32_t WB_ID, int32_t WB_POS){
	uint64_t store = ((uint64_t)WB_ID << 22) + (WB_POS & 0x3fffff);
	set_bit_48_little_endian(store, r);
}

void Single_read_aln_rst_t1_48_exact::restore(int32_t &WB_ID, int32_t &WB_POS){
	uint64_t store;
	load_bit_48_little_endian(store, r);
	WB_ID = store >> 22;
	WB_POS = store & 0x3fffff;
}

//	//move the data
void Single_read_aln_rst_unpack::store(Single_read_aln_rst_unpack & from){
	wb_ID = from.wb_ID;
	from.wb_ID = -1;
	read_hamming_mapping_position = from.read_hamming_mapping_position;
	read_len = from.read_len;
	std::swap(change_detail, from.change_detail);
	std::swap(read_string, from.read_string);
}

void Single_read_aln_rst_unpack::clear(){
	wb_ID = 0;
	read_hamming_mapping_position = 0;
	read_len = 0;
	change_detail.clear();
	read_string.clear();
}
//
void Single_read_aln_rst_unpack::print_detail(FILE * o){
	//basic:
	fprintf(o, "wb_ID %d pos %d read_len %d change_detail: N: %ld ", wb_ID,  read_hamming_mapping_position, read_len, change_detail.size());
	for(uint16_t c: change_detail)
	{
		fprintf(o, "%d, %d, %d ;",c >> 3, (c >> 3) + read_hamming_mapping_position, c &0x7);
	}
	for(uint8_t c: read_string)
		fprintf(o, "%c","ACGT"[c]);
	fprintf(o, "\n");
}

void Single_read_aln_rst_unpack::generate_read_string( void *wb_idx_p, int read_length){/////~~~~~~~~~~~~~~~~
	if(is_CPX()){//store type: CPX
		//remove C head
		int c_head_length = 0;
		for(uint8_t c: read_string){
			if(c == 1) c_head_length ++;
			else break;
		}
		if(c_head_length > 1){
			read_hamming_mapping_position += c_head_length;
			read_string.erase(read_string.begin(), read_string.begin() + c_head_length);
		}
		//remove C tail
		int c_tail_length = 0;
		for(int i = read_string.size() - 1; i >= 0; i--){
			if(read_string[i] == 1) c_tail_length ++;
			else break;
		}
		if(c_tail_length > 1){
			read_string.erase(read_string.end() - c_tail_length, read_string.end());
		}
		read_len = read_string.size();
		//read_hamming_mapping_position = -1;
		return;
	}
	//other conditions
	ALN_ONLINE::MM_idx_loader *wb_idx = (ALN_ONLINE::MM_idx_loader *)wb_idx_p;
	//loading the reference:
	int chr_ID; uint wb_bg_offset = 0;
	wb_idx->convert_wb_coordinate_to_ref_coordinate(wb_ID, chr_ID, wb_bg_offset);
	wb_bg_offset += read_hamming_mapping_position;
	wb_idx->ref.load_ref_bin_from_buff(chr_ID, wb_bg_offset, read_length, read_string);

	//debug code: ref string:
	if(false){
		fprintf(stderr, "Ref string ");
		for(uint8_t c: read_string)
			fprintf(stderr, "%c","ACGT"[c]);
		fprintf(stderr, "\n");
	}

	this->read_len = read_length;
	//todo:: change the read string to new read string;
	{
		//for CPX stored result
		//      000:0 INS: A
		//      001:1 INS: C
		//      010:2 INS: G
		//      011:3 INS: T
		//      100:4 DEL
		//      101:5 SNP: + 1
		//      102:6 SNP: + 2
		//      103:7 SNP: + 3
		for(uint16_t change: change_detail){
			uint16_t pos = (change >> 3) & 0xff; //8 bit
			if(pos >= read_string.size())
				continue;
			uint8_t change_detail = (change & 0x7);
			if(change_detail > 4){//SNP
				read_string[pos] += (change_detail - 4);
				read_string[pos] %= 4;
			}else if(change_detail < 4){//INS
				read_string.insert(read_string.begin() + pos, change_detail);
			}else{//DEL
				read_string.erase(read_string.begin() + pos);
			}
		}
		if(read_string.size() > read_length){
			read_string.resize(read_length);
		}
	}

	//debug code: read string:
	if(false){
		fprintf(stderr, "Read string ");
		for(uint8_t c: read_string)
			fprintf(stderr, "%c","ACGT"[c]);
		fprintf(stderr, "\n");
	}
}

#define MAX_STORE_NM_48_bit 1
#define MAX_STORE_NM_80_bit 4
#define MAX_STORE_NM_160_bit 11

void Single_read_aln_rst_t2_48_gap::set(int32_t wb_ID, int read_hamming_mapping_position, int16_t change_detail){
	uint64_t store =
			//WARNING:: the wb_ID must store in the upper 26 bit of the 48 bit,  this is for "SORTING by WB"
			((uint64_t)wb_ID << 22) + //26 bit
			((read_hamming_mapping_position & 0x3ff) << 12) + //10 bit
			((change_detail & 0x7ff) << 1);//11 bit
			//blank: 1 bit
	set_bit_48_little_endian(store, r);
}

void Single_read_aln_rst_t2_48_gap::restore(int32_t &WB_ID, int &read_hamming_mapping_position, uint16_t &change_detail){
	uint64_t store;
	load_bit_48_little_endian(store, r);
	WB_ID = store >> 22;
	read_hamming_mapping_position = (store >> 12) & 0x3ff; // 10 bit
	change_detail = (store >> 1) & 0x7ff;
}

void Single_read_aln_rst_t3_80_gap::set(int32_t wb_ID, int read_hamming_mapping_position, std::vector<uint16_t>& change_detail){
	//36 bit
	uint64_t store_36_p1;
	store_36_p1 = //36 bit:
			((uint64_t)wb_ID << 10) + //26 bit
			((read_hamming_mapping_position & 0x3ff)); //10 bit
	//44( = 4 * 11) bit
	xassert(change_detail.size() <= MAX_STORE_NM_80_bit, "");
//	{
//		for(uint16_t i = 0; i < change_detail.size() - 1; i++)
//			xassert((change_detail[i] >> 3) >= (change_detail[i + 1] >> 3), " ");
//	}

	uint64_t store_44_p2 = 0;
	for(uint16_t i = 0; i < MAX_STORE_NM_80_bit; i++){
		uint16_t c_c = (i < change_detail.size())?(change_detail[i]):(0x7ff);
		store_44_p2 <<= 11;
		(store_44_p2 += c_c);
	}
	//store the upper 4 bit in "store_44_p2" to the tail of store_36_p1
	store_36_p1 <<= 4; store_36_p1 += (store_44_p2 >> 40) & 0xf;
	store_44_p2 &= 0xffffffffff;
	//store the final results
	set_bit_40_little_endian(store_36_p1, r);
	set_bit_40_little_endian(store_44_p2, r + 5);
}

void Single_read_aln_rst_t3_80_gap::restore(int32_t &wb_ID, int &read_hamming_mapping_position,
		std::vector<uint16_t>& change_detail){
	uint64_t store_36_p1; uint64_t store_44_p2;
	load_bit_40_little_endian(store_36_p1, r);
	load_bit_40_little_endian(store_44_p2, r + 5);
	//re-store the upper 4 bit in "store_44_p2" from the tail of store_36_p1
	store_44_p2 += ((store_36_p1 & 0xf) << 40);
	store_36_p1 >>= 4; store_36_p1 &= 0xfffffffff;
	//load from 36 bit data
	wb_ID = (store_36_p1 >> 10); //26 but
	read_hamming_mapping_position = (store_36_p1) & 0x3ff; //10 bit
	//load from 44 bit data
	for(int i = MAX_STORE_NM_80_bit - 1; i >= 0; i--){
		uint16_t cur_change_detail = store_44_p2 & 0x7ff;
		if(cur_change_detail != 0x7ff)
			change_detail.emplace_back(cur_change_detail);
		store_44_p2 >>= 11;
	}
}

void Single_read_aln_rst_t4_160_gap::set(int32_t wb_ID, int read_hamming_mapping_position, std::vector<uint16_t>& change_detail){
	//36 bit
		uint64_t store_36_p1;
		store_36_p1 = //36 bit:
				((uint64_t)wb_ID << 10) + //26 bit
				((read_hamming_mapping_position & 0x3ff)); //10 bit
		//44( = 4 * 11) bit
		xassert(change_detail.size() <= MAX_STORE_NM_160_bit, "");

		//P1: 5 NM
		uint64_t store_44_p2[3] = {0};
		for(uint16_t i = 0; i < MAX_STORE_NM_160_bit; i++){
			uint16_t c_c = (i < change_detail.size())?(change_detail[i]):(0x7ff);
			store_44_p2[i >> 2] <<= 11;
			(store_44_p2[i >> 2] += c_c);
		}

		//store the upper 4 bit in "store_44_p2[0]" to the tail of store_36_p1
		//store the upper 4 bit in "store_44_p2[1]" to the tail of store_44_p2[2]
		uint8_t store_44_p2_0_HEAD = ((store_44_p2[0] >> 40) & 0xf); store_44_p2[0] &= 0xffffffffff;//40bit
		uint8_t store_44_p2_1_HEAD = ((store_44_p2[1] >> 40) & 0xf); store_44_p2[1] &= 0xffffffffff;//40bit
		//uint8_t store_44_p2_2_HEAD = ((store_44_p2[2] >> 40) & 0xf); store_44_p2[2] &= 0xffffffffff;//40bit
		store_36_p1 <<= 4; store_36_p1 += store_44_p2_0_HEAD;
		store_44_p2[2] <<=4; store_44_p2[2] += store_44_p2_1_HEAD;

		//store the final results
		set_bit_40_little_endian(store_36_p1, r);
		set_bit_40_little_endian(store_44_p2[0], r + 5);
		set_bit_40_little_endian(store_44_p2[1], r + 10);
		set_bit_40_little_endian(store_44_p2[2], r + 15);
}

void Single_read_aln_rst_t4_160_gap::restore(int32_t &wb_ID, int &read_hamming_mapping_position, std::vector<uint16_t>& change_detail){

	uint64_t store_36_p1; uint64_t store_44_p2[3];
	load_bit_40_little_endian(store_36_p1, r);
	load_bit_40_little_endian(store_44_p2[0], r + 5);
	load_bit_40_little_endian(store_44_p2[1], r + 10);
	load_bit_40_little_endian(store_44_p2[2], r + 15);

	//re-store the upper 4 bit in "store_44_p2" from the tail of store_36_p1
	uint64_t store_44_p2_0_HEAD = ((store_36_p1 & 0xf) << 40); 	  store_36_p1 >>= 4; 	store_36_p1 &= 0xfffffffff;
	uint64_t store_44_p2_1_HEAD = ((store_44_p2[2] & 0xf) << 40); store_44_p2[2] >>= 4; store_44_p2[2] &= 0xfffffffff;

	store_44_p2[0] += store_44_p2_0_HEAD;
	store_44_p2[1] += store_44_p2_1_HEAD;

	//load from 36 bit data
	wb_ID = (store_36_p1 >> 10); //26 but
	read_hamming_mapping_position = (store_36_p1) & 0x3ff; //10 bit
	//load from 44 bit data; store in reverse
	for(int i = MAX_STORE_NM_160_bit - 1; i >= 0; i--){
		uint16_t cur_change_detail = store_44_p2[i >> 2] & 0x7ff;
		if(cur_change_detail != 0x7ff)
			change_detail.emplace_back(cur_change_detail);
		store_44_p2[i >> 2] >>= 11;
	}
}

//
void Single_read_aln_rst_t5_cpx::set(int32_t wb_ID, int aln_position_in_WB, uint8_t* read_string, int read_len){
	xassert(read_len <= 256, "read length max support: 256bp!");
	xassert(aln_position_in_WB <= 1024, "read max support: 1024bp!");
	uint64_t store_48 = ((uint64_t)wb_ID << 22) + ((read_len & 0x3ff) << 12) + (aln_position_in_WB << 2);
	set_bit_48_little_endian(store_48, r);
	int byte_c = 0; int bit_c = 0; int base_c = 0;
	for(;base_c < read_len;){
		read_string_bin[byte_c] += ((read_string[base_c++]&0x3) << (6 - bit_c));//todo:: N is set to C
		bit_c += 2;
		if(bit_c == 8){
			bit_c = 0;
			byte_c++;
		}
	}
}

void Single_read_aln_rst_t5_cpx::restore(int32_t &wb_ID, int &aln_position_in_WB, std::vector<uint8_t> &read_string, int &read_len){
	uint64_t store_48;
	load_bit_48_little_endian(store_48, r);
	wb_ID = store_48 >> 22;
	read_len = (store_48 >> 12) & 0x3ff;
	aln_position_in_WB = (store_48 >> 2) & 0x3ff;
	int byte_c = 0; int bit_c = 0; int base_c = 0;
	for(;base_c < read_len; base_c++){
		read_string.emplace_back(((read_string_bin[byte_c] >> (6 - bit_c)) & 0x3));
		bit_c += 2;
		if(bit_c == 8){
			bit_c = 0;
			byte_c++;
		}
	}
}


//sort by the first 24 bit of the MAP_rst 48 - 3*8 = 24
#define sort_key_MAP_rst_exact_48_bit(a) ((*((uint64_t *)((a).r)) >> 24)&0xffffff)
KRADIX_SORT_INIT(map_r_x_exact_48, Single_read_aln_rst_t1_48_exact, sort_key_MAP_rst_exact_48_bit, 3)
void sortint_exact_48bit_results_for_sample(Single_read_aln_rst_t1_48_exact * v, uint64_t result_size){
	//SORT FOR exact mapping results
	radix_sort_map_r_x_exact_48(v, v + result_size);
}

//sort by the first 24 bit of the MAP_rst 48 - 3*8 = 24
#define sort_key_MAP_rst_gap_48_bit(a) ((*((uint64_t *)((a).r)) >> 24)&0xffffff)
KRADIX_SORT_INIT(map_r_x_gap_48, Single_read_aln_rst_t2_48_gap, sort_key_MAP_rst_gap_48_bit, 3)
void sortint_gap_48bit_results_for_sample(Single_read_aln_rst_t2_48_gap * v, uint64_t result_size){
	//SORT FOR exact mapping results
	radix_sort_map_r_x_gap_48(v, v + result_size);
}

//sort by the first 24 bit of the MAP_rst 40 - 3*8 = 16
#define sort_key_MAP_rst_gap_80_bit(a) ((*((uint64_t *)((a).r)) >> 16)&0xffffff)
KRADIX_SORT_INIT(map_r_x_gap_80, Single_read_aln_rst_t3_80_gap, sort_key_MAP_rst_gap_80_bit, 3)
void sortint_gap_80bit_results_for_sample(Single_read_aln_rst_t3_80_gap * v, uint64_t result_size){
	//SORT FOR exact mapping results
	radix_sort_map_r_x_gap_80(v, v + result_size);
}

//sort by the first 24 bit of the MAP_rst: 40 - 3*8 = 16
#define sort_key_MAP_rst_gap_160_bit(a) ((*((uint64_t *)((a).r)) >> 16)&0xffffff)
KRADIX_SORT_INIT(map_r_x_gap_160, Single_read_aln_rst_t4_160_gap, sort_key_MAP_rst_gap_160_bit, 3)
void sortint_gap_160bit_results_for_sample(Single_read_aln_rst_t4_160_gap * v, uint64_t result_size){
	//SORT FOR exact mapping results
	radix_sort_map_r_x_gap_160(v, v + result_size);
}

//sort by the first 24 bit of the MAP_rst 48 - 3*8 = 24
#define sort_key_MAP_rst_gap_CPX_bit(a) ((*((uint64_t *)((a).r)) >> 24)&0xffffff)
KRADIX_SORT_INIT(map_r_x_gap_cpx, Single_read_aln_rst_t5_cpx, sort_key_MAP_rst_gap_CPX_bit, 3)
void sortint_gap_cpx_results_for_sample(Single_read_aln_rst_t5_cpx * v, uint64_t result_size){
	//SORT FOR exact mapping results
	radix_sort_map_r_x_gap_cpx(v, v + result_size);
}

}
