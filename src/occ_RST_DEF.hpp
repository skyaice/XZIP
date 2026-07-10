/*
 * occ_RST_DEF.hpp
 *
 *  Created on: 2022年12月10日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_OCC_RST_DEF_HPP_
#define SRC_BWT_ONLINE_OCC_RST_DEF_HPP_

extern "C"
{
#include "clib/utils.h"
}

#include "BWT_idx/bwt.hpp"

//total 64 bit:
struct Gap_aln_Cluster_t{
    uint32_t WB_ID;//26 bit
    uint32_t is_rev:1;//1 bit
    uint32_t score:10;//10 bit
    uint32_t read_to_wb_offset:10;//10 bit
    Gap_aln_Cluster_t(uint32_t WB_ID, uint32_t is_rev, uint32_t score, uint32_t read_to_wb_offset){
        this->WB_ID = WB_ID;
        this->is_rev = is_rev;
        this->score = score;
        this->read_to_wb_offset = read_to_wb_offset;
    }

    Gap_aln_Cluster_t(){
        WB_ID = 0;
        is_rev = 0;
        score = 0;
        read_to_wb_offset = 0;
    }

    Gap_aln_Cluster_t(Gap_aln_Cluster_t * c){
        this->WB_ID = c->WB_ID;
        this->is_rev = c->is_rev;
        this->score = c->score;
        this->read_to_wb_offset = c->read_to_wb_offset;
    }

    void copy(Gap_aln_Cluster_t &c){
        this->WB_ID = c.WB_ID;
        this->is_rev = c.is_rev;
        this->score = c.score;
        this->read_to_wb_offset = c.read_to_wb_offset;
    }
};

struct SSD_FILE_BLOCK_META{
    uint64_t __SSD__FILE__OFFSET__;
    uint32_t __SSD__FILE__SIZE__;

    SSD_FILE_BLOCK_META(uint64_t __SSD__FILE__OFFSET__, uint32_t __SSD__FILE__SIZE__){
        this->__SSD__FILE__OFFSET__ = __SSD__FILE__OFFSET__;
        this->__SSD__FILE__SIZE__ = __SSD__FILE__SIZE__;
    }
};

//bit-level functions
//PART1: basic value: for cluster store type and BWT grow index
int calculate_gap_aln_store_type(uint64_t cluster_number);

void set_gap_bwt_GROW_idx(
    //input:
    uint32_t gap_bwt_GROW_idx,
    //output:
    bwtint_t &new_k);

void set_gap_aln_store_type(
    //input:
    uint32_t gap_aln_store_type,
    //output:
    bwtint_t &new_k);

void set_gap_bwt_grow_length(
    //input:
    uint8_t gap_bwt_grow_length,
    //output:
    bwtint_t &new_k);

void get_gap_aln_metadata(
        //input:
        bwtint_t old_k,
        //output:
        uint32_t &gap_bwt_GROW_idx, int32_t &gap_aln_store_type, uint8_t &gap_bwt_grow_length
        );

int32_t get_gap_aln_store_type_FOR_PAIRING(
        //input:
        bwtint_t old_k, bwtint_t old_l
        );


//PART2: for uniq clusters, store directly in the same place
void set_uniq_cluster_data(
        //input:
        Gap_aln_Cluster_t &c,
        //output:
        bwtint_t &new_k, bwtint_t &new_l);

void get_uniq_cluster_data(
        //input:
        bwtint_t old_k, bwtint_t old_l,
        //output:
        std::vector<Gap_aln_Cluster_t> & c_v);

void get_uniq_cluster_data_FOR_PAIRING(
        //input:
        bwtint_t old_k, bwtint_t old_l,
        //output:
        int32_t &wb_id, int8_t &is_rev, int16_t &score, uint16_t &read_to_wb_offset);

//PART3: for small clusters
void set_small_cluster_meta(
        //input:
        uint32_t cluster_store_bg_idx,
        uint32_t cluster_store_ed_idx,
        //output:
        bwtint_t &new_l);

void get_small_cluster_meta(
        //input:
        bwtint_t old_l,
        //output:
        uint32_t &cluster_store_bg_idx,
        uint32_t &cluster_store_ed_idx);

//PART4: for BIG clusters
void set_big_cluster_meta(
        //input:
        uint64_t __SSD__FILE__OFFSET__,
        uint32_t __SSD__FILE__SIZE__,
        //output
        bwtint_t &old_k,
        bwtint_t &old_l
);

void get_big_cluster_meta(
        //input:
        bwtint_t new_k, bwtint_t new_l,
        //output:
        uint64_t &__SSD__FILE__OFFSET__,
        uint32_t &__SSD__FILE__SIZE__
        );


void set_big_cluster_meta_FOR_SA(
        //input:
        uint32_t big_cluster_store_bg_idx,
        uint32_t big_cluster_store_ed_idx,
        //output:
        uint64_t &k,
        uint64_t &l);

void get_big_cluster_meta_FOR_SA(
        //output:
        uint32_t &big_cluster_store_bg_idx,
        uint32_t &big_cluster_store_ed_idx,
        //input:
        uint64_t k,
        uint64_t l);


#endif /* SRC_BWT_ONLINE_OCC_RST_DEF_HPP_ */
