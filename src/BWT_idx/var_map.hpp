/*
 * var_map.hpp
 *
 *  Created on: 2022年8月7日
 *      Author: fenghe
 */

#ifndef VAR_MAP_HPP_
#define VAR_MAP_HPP_

#include <cstring>
#include <vector>
#include <map>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include "../CPPLIB//tools.hpp"
extern "C" {
#include "../clib/utils.h"
//#include "../clib/uthash.h"
}

struct VCF_item{
    uint32_t chr_ID;
    uint32_t var_pos;
    const char *ref_p;
    const char *alt_p;
    int QUAL;
    int PASS_TYPE; // 0: NOT_PASS;1: PASS; 2:LOW_COV
    int GT_TYPE;// 0: 0/0; 1: 0/1; 2: 1/1
    int var_len;
    int ref_supp;
    int var_supp;

    uint32_t var_ID;
    uint32_t var_wb;

    int is_novel_var;// 0: known; 1: novel+CHANGE 2: novel+ASS;

public:

    VCF_item(){}

    void set(
            uint32_t chr_ID_,uint32_t var_pos_,
            const char *ref_p_,const char *alt_p_,
            int QUAL_, int PASS_TYPE_, int GT_TYPE_,
            int var_len_,
            int ref_supp_,int var_supp_,
            int var_ID_,int var_wb_,
            int is_novel_var
            ){
        this->chr_ID = chr_ID_; this->var_pos = var_pos_;
        this->ref_p = ref_p_; this->alt_p = alt_p_;
        this->QUAL = QUAL_;
        this->PASS_TYPE = PASS_TYPE_; this->GT_TYPE = GT_TYPE_;
        this->var_len = var_len_;
        this->ref_supp = ref_supp_; this->var_supp = var_supp_;

        this->var_ID = var_ID_;
        this->var_wb = var_wb_;

        this->is_novel_var = is_novel_var;
    }

    bool var_is_same(VCF_item &c){
        return (this->var_pos == c.var_pos && strcmp(this->ref_p, c.ref_p) == 0 && strcmp(this->alt_p, c.alt_p) == 0);
    }

    void VCF_dump(FILE* out){
        VCF_dump_core(out, "debug");
    }

    void VCF_dump_core(FILE* out, const char * var_name_buff){

        if(chr_ID < 22) 		fprintf(out, "chr%d\t", chr_ID + 1); 			//CHROM
        else if(chr_ID == 22) 	fprintf(out, "chrX\t"); 			//CHROM
        else if(chr_ID == 23)	fprintf(out, "chrY\t"); 			//CHROM
        else if(chr_ID == 24)	fprintf(out, "chrM\t"); 			//CHROM
        else{					xassert(0, "UNKNOWN CHR");}

        if(82007683 == var_pos)
            fprintf(stderr, " ");

        fprintf(out, "%d\t", var_pos);		//POS
        fprintf(out, "%s\t", var_name_buff);//ID
        fprintf(out, "%s\t", ref_p);		//REF
        fprintf(out, "%s\t", alt_p);		//ALT
        fprintf(out, "30\t");				//QUAL::30
        switch(PASS_TYPE){					//FILTER
        case 0: fprintf(out, "NOT_PASS\t");	break;
        case 1: fprintf(out, "PASS\t");	break;
        case 2: fprintf(out, "LOW_COV\t");	break;
        default:fprintf(out, "NOT_PASS\t");	break;
        }
        //INFO
        if(var_len == 0)	{		fprintf(out, "TYPE=SNP\t");		}
        else if(var_len > 0){		fprintf(out, "TYPE=INS\t");		}
        else 				{		fprintf(out, "TYPE=DEL\t");		}
        fprintf(out, "GT:DR:DA:NK\t");//									//FORMAT
        switch(GT_TYPE){											//SAMPLE::GT
        case 0: fprintf(out, "0/0:");	break;
        case 1: fprintf(out, "0/1:");	break;
        case 2: fprintf(out, "1/1:");	break;
        default:fprintf(out, "0/0:");	break;
        }
        fprintf(out, "%d:%d:", ref_supp, var_supp);		//SAMPLE::DP
        switch(is_novel_var){
        case 0: fprintf(out, "KNOWN"); break;
        case 1: fprintf(out, "NOVEL_ALN"); break;
        case 2: fprintf(out, "NOVEL_ASS"); break;
        }

        fprintf(out, "\n");		//SAMPLE::DP
    }
};

struct Variant{
private:
    //std::string ref_str;
    //std::string alt_str;
    union{
        char * ref_alt_buff;//
        char ref_alt_buff_1[8];
    };

    int8_t ref_len;
    int8_t alt_len;
    uint8_t REF_ALT_SAME_len;

public:
    uint32_t ref_pos;
    void set(uint32_t ref_pos_, std::string &ref_str_, std::string &alt_str_, uint32_t REF_ALT_SAME_len_){
        ref_pos = ref_pos_;
        if(ref_str_.size() > 50) {ref_str_.clear();  ref_str_.push_back('N');}
        if(alt_str_.size() > 50) {alt_str_.clear();  alt_str_.push_back('N');}

        ref_len = ref_str_.size();
        alt_len = alt_str_.size();
        char *ref_alt = NULL;
        if(alt_len + ref_len <= 6){
            ref_alt = ref_alt_buff_1;
        }else{
            ref_alt_buff = (char *)xmalloc(ref_len + alt_len + 2);
            ref_alt = (ref_alt_buff);
        }
        memcpy(ref_alt,ref_str_.c_str(),ref_len);
        ref_alt[ref_len] = 0;
        memcpy(ref_alt + ref_len + 1, alt_str_.c_str(), alt_len);
        ref_alt[ref_len + 1 + alt_len] = 0;
        if(REF_ALT_SAME_len_ > 255)
            REF_ALT_SAME_len = 255;
        else
            REF_ALT_SAME_len = REF_ALT_SAME_len_;

    }
    int get_length(){ return alt_len - ref_len; }
    int get_REF_ALT_SAME_len(){ return REF_ALT_SAME_len; }
    bool isSAME_ref_alt(Variant & B){ return ((strcmp(get_ref(),B.get_ref()) == 0) && (strcmp(get_alt(),B.get_alt()) == 0)); }
    bool isNoVar(){ return (strcmp(get_ref(),get_alt()) == 0); }
    inline const char *get_ref(){
        if(alt_len + ref_len <= 6){
            return ref_alt_buff_1;
        }else{
            return ref_alt_buff;
        }
    }
    inline const char *get_alt(){
        return get_ref() + ref_len + 1;
    }

    int8_t get_ref_len(){return ref_len;}
    int8_t get_alt_len(){return alt_len;}

    void print(FILE *log_f){
        fprintf(log_f, " [pos %u ", ref_pos);
        fprintf(log_f, " ref %s  ", get_ref());
        fprintf(log_f, " alt %s ] \t", get_alt());
    }
} ;

typedef struct {
    std::vector<uint16_t> variant_list;
    uint32_t hap_length;
    void print(){
        fprintf(stderr, "hap_length: %u ", hap_length);
        for(uint32_t var : variant_list)
            fprintf(stderr, "var ID %u\t", var );
        fprintf(stderr, "\n");
    }

} hap_t;

struct Window_Repeat_region{ int st; int ed; Window_Repeat_region(int st_, int ed_):st(st_),ed(ed_){} };

struct Window_t{
    Window_t(){
        chr_ID = 0;
        st_pos = 0;
        total_hap_length = 300;
        hap_list.emplace_back();
        hap_list.back().hap_length = 300;
        repeat_region_len = 0;
    }
    std::vector<hap_t> hap_list;
    std::vector<Variant> var_list;

    //repeat regions
    //Window_Repeat_region *repeat_v;
    uint32_t repeat_region_len;


    // typedef struct {
    // 	char key[50];  // 假设键的最大长度为 10
    // 	int value;
    // 	UT_hash_handle hh;  // 使该结构体可哈希
    // } HashTableEntry;
    std::unordered_map<std::string, int> variant_hash;
    void build_variant_hash(){
        std::string alt;
        variant_hash.clear();
        variant_hash.reserve(var_list.size());

        for(int i =0; i<var_list.size();i++){
            Variant &var_item = var_list[i];
            uint32_t pos = var_item.ref_pos;
            int8_t var_length = var_item.get_length();
            if(var_length>0){//ins
                alt = var_item.get_alt()[var_item.get_alt_len()-var_length, var_item.get_alt_len()];

            }else if(var_length<0){//del

                alt = var_item.get_ref()[ var_item.get_ref_len()+var_length, var_item.get_ref_len()];

            }else{//mismatch
                alt = var_item.get_alt();
            }
            char key_str[50];
            memset(key_str,0,50);
            snprintf(key_str, sizeof(key_str), "%d_%s", var_item.ref_pos, alt.c_str());
            //std::string s = std::string(key_str);
            variant_hash.insert(std::pair<std::string, int>(key_str, i));


        }
    }
    int get_variant(std::string alt_str){
        auto it = variant_hash.find(alt_str);
        if(it!=variant_hash.end()){
            return it->second;
        }
        return -1;

    }

    void get_hap_string(int hap_ID, std::string &ref, std::string &hap_string) {
        bool trace_this_call = (hap_ID == 2 && ref.size() == 300 && var_list.size() == 25);

        if (trace_this_call) {
            std::cerr << "\n[GHS] ===== ENTER get_hap_string =====" << std::endl;
            std::cerr << "[GHS] hap_ID=" << hap_ID
                      << ", hap_list.size()=" << hap_list.size()
                      << ", var_list.size()=" << var_list.size()
                      << ", ref.size()=" << ref.size() << std::endl;
        }

        hap_string = ref;
        int ALT_store_offset = 0;

        if (hap_ID < 0 || hap_ID >= (int)hap_list.size()) {
            std::cerr << "[GHS][FATAL] hap_ID out of range: hap_ID=" << hap_ID
                      << ", hap_list.size()=" << hap_list.size() << std::endl;
            return;
        }

        std::vector<uint16_t> &var_l = hap_list[hap_ID].variant_list;

        if (trace_this_call) {
            std::cerr << "[GHS] variant_list.size()=" << var_l.size() << std::endl;
        }

        for (size_t var_id = 0; var_id < var_l.size(); var_id++) {
            uint16_t real_var_idx = var_l[var_id];

            if (trace_this_call) {
                std::cerr << "[GHS] --- var_id=" << var_id
                          << ", real_var_idx=" << real_var_idx
                          << ", var_list.size()=" << var_list.size()
                          << ", ALT_store_offset=" << ALT_store_offset
                          << std::endl;
            }

            if (real_var_idx >= var_list.size()) {
                std::cerr << "[GHS][FATAL] var_list index out of range: real_var_idx="
                          << real_var_idx << ", var_list.size()=" << var_list.size()
                          << ", hap_ID=" << hap_ID << ", var_id=" << var_id << std::endl;
                return;
            }

            Variant &c_var = var_list[real_var_idx];

            int var_len = c_var.get_length();
            int ref_pos = c_var.ref_pos;
            int alt_store_pos = ref_pos + ALT_store_offset;

            const char *ref_p = c_var.get_ref();
            const char *alt_p = c_var.get_alt();

            if (ref_p == nullptr) {
                std::cerr << "[GHS][FATAL] c_var.get_ref() returned nullptr"
                          << " (hap_ID=" << hap_ID
                          << ", var_id=" << var_id
                          << ", real_var_idx=" << real_var_idx << ")" << std::endl;
                return;
            }

            if (alt_p == nullptr) {
                std::cerr << "[GHS][FATAL] c_var.get_alt() returned nullptr"
                          << " (hap_ID=" << hap_ID
                          << ", var_id=" << var_id
                          << ", real_var_idx=" << real_var_idx << ")" << std::endl;
                return;
            }

            size_t ref_cstr_len = std::strlen(ref_p);
            size_t alt_cstr_len = std::strlen(alt_p);
            int ref_len_api = c_var.get_ref_len();

            if (trace_this_call) {
                std::cerr << "[GHS] ref_pos=" << ref_pos
                          << ", alt_store_pos=" << alt_store_pos
                          << ", var_len=" << var_len
                          << ", ref_cstr_len=" << ref_cstr_len
                          << ", alt_cstr_len=" << alt_cstr_len
                          << ", ref_len_api=" << ref_len_api
                          << ", hap_string.size()=" << hap_string.size()
                          << std::endl;
                std::cerr << "[GHS] ref_p=\"" << ref_p << "\"" << std::endl;
                std::cerr << "[GHS] alt_p=\"" << alt_p << "\"" << std::endl;
            }

            if (var_len == 0) {
                if (alt_store_pos < 0 || alt_store_pos >= (int)hap_string.size()) {
                    if (trace_this_call) {
                        std::cerr << "[GHS][WARN] SNP/MNP alt_store_pos out of range, continue"
                                  << std::endl;
                    }
                    continue;
                }

                size_t replace_len = std::min(ref_cstr_len, alt_cstr_len);

                if (trace_this_call && ref_cstr_len != alt_cstr_len) {
                    std::cerr << "[GHS][WARN] ref/alt cstring length mismatch in SNP/MNP: "
                              << "ref_cstr_len=" << ref_cstr_len
                              << ", alt_cstr_len=" << alt_cstr_len << std::endl;
                }

                size_t max_writable = hap_string.size() - (size_t)alt_store_pos;
                replace_len = std::min(replace_len, max_writable);

                for (size_t i = 0; i < replace_len; i++) {
                    hap_string[alt_store_pos + i] = alt_p[i];
                }
            } else {
                if (alt_store_pos < 0 || alt_store_pos > (int)hap_string.size()) {
                    if (trace_this_call) {
                        std::cerr << "[GHS][WARN] INDEL alt_store_pos out of range, continue"
                                  << std::endl;
                    }
                    continue;
                }

                if (ref_len_api < 0) {
                    std::cerr << "[GHS][FATAL] c_var.get_ref_len() < 0 : " << ref_len_api
                              << " (hap_ID=" << hap_ID
                              << ", var_id=" << var_id
                              << ", real_var_idx=" << real_var_idx << ")" << std::endl;
                    return;
                }

                ALT_store_offset += var_len;

                if (trace_this_call) {
                    std::cerr << "[GHS] before erase: pos=" << alt_store_pos
                              << ", erase_len=" << ref_len_api
                              << ", hap_string.size()=" << hap_string.size() << std::endl;
                }

                hap_string.erase((size_t)alt_store_pos, (size_t)ref_len_api);

                if (trace_this_call) {
                    std::cerr << "[GHS] after erase: hap_string.size()="
                              << hap_string.size() << std::endl;
                    std::cerr << "[GHS] before insert: pos=" << alt_store_pos
                              << ", alt=\"" << alt_p << "\"" << std::endl;
                }

                hap_string.insert((size_t)alt_store_pos, alt_p);

                if (trace_this_call) {
                    std::cerr << "[GHS] after insert: hap_string.size()="
                              << hap_string.size() << std::endl;
                }
            }
        }

        if (trace_this_call) {
            std::cerr << "[GHS] SUCCESS, final hap_string.size()="
                      << hap_string.size() << std::endl;
            std::cerr << "[GHS] ===== LEAVE get_hap_string =====\n" << std::endl;
        }
    }

    void store_original_repeat_region(std::vector<std::string> & header_line, std::vector<Window_Repeat_region> &original_repeat_region){
        original_repeat_region.clear();
        uint32_t repeat_ori_BUFF_size = atoi(header_line[4].c_str());
        if(repeat_ori_BUFF_size == 0)
            return;
        //store the original regions
        for(uint rr_idx = 0; rr_idx < repeat_ori_BUFF_size; rr_idx++)
            original_repeat_region.emplace_back(atoi(header_line[5 + rr_idx*2].c_str()), atoi(header_line[6 + rr_idx*2].c_str()));
    }

    void store_global_repeat_region(std::vector<Window_Repeat_region> &original_repeat_region){
        if(original_repeat_region.empty())
            return;
        for(uint rep_ID = 0; rep_ID < original_repeat_region.size(); rep_ID ++)
            repeat_region_len += original_repeat_region[rep_ID].ed - original_repeat_region[rep_ID].st;
    }

    uint32_t chr_ID;
    uint32_t st_pos;
    uint32_t total_hap_length;

    void print(FILE * log_f){
        uint32_t hap_id = 0;
        uint32_t bg_pos = 0;
        for(auto &hap: hap_list){
            uint32_t ed_pos = bg_pos + hap.hap_length;
            fprintf(log_f, "hap_id %d, bg_pos %d, ed_pos %d  ",hap_id, bg_pos, ed_pos);
            fprintf(log_f, "hap_length: %u ", hap.hap_length);
            for(uint32_t var : hap.variant_list){
                var_list[var].print(log_f);
            }
            fprintf(log_f, "\n");
            hap_id++;
            bg_pos += hap.hap_length;
        }
        // var list:
        for(auto &var: var_list){
            var.print(log_f);
            fprintf(log_f, "\n");
        }
    }

    static void load_variant_map(char *fn, uint32_t DEBUG_MIN_WB_INX, uint32_t DEBUG_MAX_WB_INX, std::vector<Window_t> &window_info){
            DEBUG_MAX_WB_INX += 10;
            char *temp = (char *)xmalloc(1000000); //1M
            std::vector<std::string> variant_tmp;
            std::vector<std::string> var_split_tmp;
            //std::vector<std::string> var_map_data;
            //try to load map file
            {
                //get file names
                FILE * try_f = xopen(fn, "rb");
                fclose(try_f);
            }
            char *map_load_temp = new char[MAX_LINE_LENGTH];//10M
            std::ifstream map_f(fn);
            //load the title line:
            map_f.getline(map_load_temp, MAX_LINE_LENGTH);
            split_string(variant_tmp, temp, map_load_temp , " ");
            uint32_t window_block_number  = atoi(variant_tmp[1].c_str());
            if(true)
                window_block_number = MIN(DEBUG_MAX_WB_INX, window_block_number); //debug:
            window_info.resize(window_block_number);
            std::vector<Window_Repeat_region> original_repeat_region_BUFF;
            //load_strings_from_file(fn, var_map_data, MAX_MAP_FILE_LINE_NUM);
            while(true){
                map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                if(map_f.eof() && (map_load_temp[0]) == 0) break;

                split_string(variant_tmp, temp, map_load_temp , " ");
                uint32_t window_id  = atoi(variant_tmp[0].c_str());
                if(window_id >= DEBUG_MAX_WB_INX) break;//debug:
                auto & c_wb = window_info[window_id];
                c_wb.chr_ID  = atoi(variant_tmp[1].c_str());
                c_wb.st_pos  = atoi(variant_tmp[2].c_str());
                bool wb_skip = (window_id  < DEBUG_MIN_WB_INX );
                if(wb_skip){
                    uint32_t hap_number  = atoi(variant_tmp[3].c_str());
                    for(uint hap_idx = 0; hap_idx < hap_number; hap_idx++)
                        map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                    map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                    split_string(variant_tmp, temp, map_load_temp , " ");
                    uint32_t var_number  = atoi(variant_tmp[1].c_str());
                    for(uint var_idx = 0; var_idx < var_number; var_idx++)
                        map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                    continue;
                }

                //repeat original repeat regions
                c_wb.store_original_repeat_region(variant_tmp, original_repeat_region_BUFF);
                //haplotypes
                uint32_t hap_number  = atoi(variant_tmp[3].c_str());
                c_wb.hap_list.resize(hap_number + 1);
                auto & c_hap_v = c_wb.hap_list;
                for(uint hap_idx = 0; hap_idx < hap_number; hap_idx++){
                    map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                    split_string(variant_tmp, temp, map_load_temp , " ");
                    c_hap_v[hap_idx].hap_length =  atoi(variant_tmp[0].c_str());
                    c_wb.total_hap_length += c_hap_v[hap_idx].hap_length;
                    //c_hap_v[hap_idx].variant_list.resize(variant_tmp.size() - 1);
                    for(uint i = 1; i < variant_tmp.size(); i++){
                        //c_hap_v[hap_idx].variant_list[i - 1] = atoi(variant_tmp[i].c_str());
                        c_hap_v[hap_idx].variant_list.emplace_back( atoi(variant_tmp[i].c_str()));
                    }
                }
                //adding haplotype for pure reference
                c_hap_v[hap_number].hap_length = 300;
                //store var list
                map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                split_string(variant_tmp, temp, map_load_temp , " ");
                uint32_t var_number  = atoi(variant_tmp[1].c_str());
                c_wb.var_list.resize(var_number);
                auto & c_var_v = c_wb.var_list;
                for(uint var_idx = 0; var_idx < var_number; var_idx++){
                    map_f.getline(map_load_temp, MAX_LINE_LENGTH);
                    split_string(variant_tmp, temp, map_load_temp , " ");
                    uint32_t var_id = atoi(variant_tmp[2].c_str());
                    split_string(var_split_tmp, temp, variant_tmp[1].c_str() , "_");
                    uint32_t REF_ALT_SAME_len = 0;
                    if(var_split_tmp.size() > 3)
                        REF_ALT_SAME_len = atoi(var_split_tmp[3].c_str());
                    c_var_v[var_id].set(atoi(var_split_tmp[2].c_str()), var_split_tmp[0], var_split_tmp[1], REF_ALT_SAME_len);
                }
                //store global region after all info is stored
                c_wb.store_global_repeat_region(original_repeat_region_BUFF);
            }
            map_f.close();
            free(temp);
        }

};




#endif /* VAR_MAP_HPP_ */
