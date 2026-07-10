/*
 * aln_result_struct.hpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_ALN_RST_DEF_HPP_
#define SRC_BWT_ONLINE_ALN_RST_DEF_HPP_

#include<vector>
#include<stdlib.h>
extern "C"
{
#include "clib/utils.h"
}

namespace Aln_online {

void set_bit_48_little_endian(uint64_t store, uint8_t *r);
void load_bit_48_little_endian(uint64_t &store, uint8_t *r);
void set_bit_40_little_endian(uint64_t store, uint8_t *r);
void load_bit_40_little_endian(uint64_t &store, uint8_t *r);
void set_bit_32_little_endian(uint64_t store, uint8_t *r);
void load_bit_32_little_endian(uint64_t &store, uint8_t *r);

//store EXACT aln results
struct Single_read_aln_rst_t1_48_exact{
    uint8_t r[6];
    //26 + 22
    void set(int32_t WB_ID, int32_t WB_POS);
    void restore(int32_t &WB_ID, int32_t &WB_POS);
};

struct Single_read_aln_rst_unpack{
    int32_t wb_ID;
    int read_hamming_mapping_position;
    int read_len;
    //for GAP stored results
    std::vector<uint16_t> change_detail;
    //for CPX stored result
    std::vector<uint8_t> read_string;

    Single_read_aln_rst_unpack(){
        wb_ID = -1;
        read_hamming_mapping_position = -1;
        read_len = -1;
    }

//	//move the data
    void store(Single_read_aln_rst_unpack & from);
    void clear();
    bool is_CPX(){	return (!read_string.empty());}
    void generate_read_string(void *wb_idx, int read_len);
    void print_detail(FILE * o);
};

//store GAP aln (NM ==1) results
// total: 48 bit
// wb_ID： 26 bit
// read_hamming_mapping_position: 10 bit
// check bits: 2bit
// change_detail: 10 bit

struct Single_read_aln_rst_t2_48_gap{
    uint8_t r[6];
    void set(int32_t wb_ID, int read_hamming_mapping_position, int16_t change_detail);
    void restore(int32_t &WB_ID, int &read_hamming_mapping_position, uint16_t &change_detail);
    void unpack(Single_read_aln_rst_unpack &r, int default_read_len){
        r.clear();  r.change_detail.emplace_back();
        restore(r.wb_ID, r.read_hamming_mapping_position, r.change_detail[0]);
        r.read_len = default_read_len;
    }
};

//store GAP aln (NM == 2/3/4) results
// total: 80 bit
// wb_ID： 26 bit
// read_hamming_mapping_position: 10 bit
// check bits: 2 bit
// blank: 2bit:
// change_detail1: 11 bit
// change_detail2: 11 bit
// change_detail3: 11 bit
// change_detail4: 11 bit

struct Single_read_aln_rst_t3_80_gap{
    uint8_t r[10];
//	48 bit； 32 + 16
//	80 bit； 32 + 48
//32 ==:
//	WB： 22 bit；
//	POS： 10 bit；（MAX 1024 bp）
    void set(int32_t wb_ID, int read_hamming_mapping_position, std::vector<uint16_t>& change_detail);
    void restore(int32_t &wb_ID, int &read_hamming_mapping_position, std::vector<uint16_t>& change_detail);
    void unpack(Single_read_aln_rst_unpack &r, int default_read_len){
        r.clear();
        restore(r.wb_ID, r.read_hamming_mapping_position, r.change_detail);
        r.read_len = default_read_len;
    }
};

//store GAP aln (NM == 5~11) results
// total: 160 bit
//meta: 36 bit
// wb_ID： 26 bit
// read_hamming_mapping_position: 10 bit

// change_detail1: 11 bit
// change_detail2: 11 bit
// change_detail3: 11 bit
// change_detail4: 11 bit
// change_detail5: 11 bit
// change_detail6: 11 bit
// change_detail7: 11 bit
// change_detail8: 11 bit
// change_detail9: 11 bit
// change_detail10: 11 bit
// change_detail11: 11 bit

//total: 36 + 121 = 157
struct Single_read_aln_rst_t4_160_gap{
    uint8_t r[20];
//	48 bit； 32 + 16
//	80 bit； 32 + 48
//32 ==:
//	WB： 22 bit；
//	POS： 10 bit；（MAX 1024 bp）
    void set(int32_t wb_ID, int read_hamming_mapping_position, std::vector<uint16_t>& change_detail);
    void restore(int32_t &wb_ID, int &read_hamming_mapping_position, std::vector<uint16_t>& change_detail);
    void unpack(Single_read_aln_rst_unpack &r, int default_read_len){
        r.clear();
        restore(r.wb_ID, r.read_hamming_mapping_position, r.change_detail);
        r.read_len = default_read_len;
    }
};


//store GAP aln (NM => 5 or UM) results
//meta: (48 byte in total)
//wb_ID : 26 bit
//read length: 10 bit
//suggest read alignment position: 10 bit
//INFO: for many unmapped reads, the suggest alignment position is existed, for others, the value will be set to 0x3ff;
//blank: 2 bit
//[64byte:] at most support: 256 byte data
//
struct Single_read_aln_rst_t5_cpx{
    uint8_t r[6];//AT least 6 byte
    uint8_t read_string_bin[64];
    void set(int32_t wb_ID, int aln_position_in_WB, uint8_t* read_string, int read_len);
    void restore(int32_t &wb_ID, int &aln_position_in_WB, std::vector<uint8_t> &read_string, int &read_len);
    void unpack(Single_read_aln_rst_unpack &r){
        r.clear();
        restore(r.wb_ID, r.read_hamming_mapping_position, r.read_string, r.read_len);
    }
};

struct Single_sample_align_rst_in_1M_block_BUFF{
    std::vector<Single_read_aln_rst_t1_48_exact> exact_rst_list;
    std::vector<Single_read_aln_rst_t2_48_gap> gap_NM_ONE_rst_list;
    std::vector<Single_read_aln_rst_t3_80_gap> gap_NM_SMALL_rst_list;
    std::vector<Single_read_aln_rst_t4_160_gap> gap_NM_MIDDLE_rst_list;
    std::vector<Single_read_aln_rst_t5_cpx> gap_NM_BIG_rst_list;
    void clear(){
        exact_rst_list.clear();
        gap_NM_ONE_rst_list.clear();
        gap_NM_SMALL_rst_list.clear();
        gap_NM_MIDDLE_rst_list.clear();
        gap_NM_BIG_rst_list.clear();
    }
};

struct OCC_rst_80{
    uint8_t bwt_bg[5];
    uint8_t bwt_ed[5];
    void set(uint64_t k, uint64_t l){
        set_bit_40_little_endian(k, bwt_bg);
        set_bit_40_little_endian(l, bwt_ed);
    }
    void restore(uint64_t &k, uint64_t &l){
        load_bit_40_little_endian(k, bwt_bg);
        load_bit_40_little_endian(l, bwt_ed);
    }
    inline bool is_exact_aln(){ return (bwt_bg[4] <= bwt_ed[4]); }
};

//sort by the first 24 bit of the MAP_rst
void sortint_exact_48bit_results_for_sample(Single_read_aln_rst_t1_48_exact * v, uint64_t result_size);
//sort by the first 24 bit of the MAP_rst
void sortint_gap_48bit_results_for_sample(Single_read_aln_rst_t2_48_gap * v, uint64_t result_size);
//sort by the first 24 bit of the MAP_rst
void sortint_gap_80bit_results_for_sample(Single_read_aln_rst_t3_80_gap * v, uint64_t result_size);
//sort by the first 24 bit of the MAP_rst
void sortint_gap_160bit_results_for_sample(Single_read_aln_rst_t4_160_gap * v, uint64_t result_size);
//sort by the first 24 bit of the MAP_rst
void sortint_gap_cpx_results_for_sample(Single_read_aln_rst_t5_cpx * v, uint64_t result_size);

}  // namespace Aln_online


#endif /* SRC_BWT_ONLINE_ALN_RST_DEF_HPP_ */
