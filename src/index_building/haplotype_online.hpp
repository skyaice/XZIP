/*
 * haplotype.hpp
 *
 *  Created on: 2022年4月11日
 *      Author: fenghe
 */

#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <zlib.h>
#include "../CPPLIB/tools.hpp"
#include <fstream>
extern "C"{
    #include "../clib/utils.h"
    #include "../clib/bam_file.h"
    #include "../clib/desc.h"
}
#include "../CPPLIB/get_option_cpp.hpp"
#include "../minimizer/minimizer.hpp"
//#include "deBGA_index.hpp"

#define vcf_ref_line 3
#define vcf_alt_line 4
#define filter_line_number 6
#define INFO_line_idx 7
#define genotype_begin_index 9

#define WITH_VAR 1
#define WITHOUT_VAR 0

namespace ALN_ONLINE{
struct GT_analysis{

    GT_analysis(std::vector<std::string> &item_value, int chrID_, int pos_){
        POS = pos_;
        chrID = chrID_;
        REF.clear(); REF.append(item_value[vcf_ref_line]);
        ALT.clear(); ALT.append(item_value[vcf_alt_line]);
        int sample_numer = item_value.size() - genotype_begin_index;
        if(item_value.back().size() < 3){ sample_numer -= 1; }
        GT.clear();
        for(int i = 0; i < sample_numer; i++){
            char GT1 = item_value[genotype_begin_index + i][0];
            char GT2 = item_value[genotype_begin_index + i][2];
            if(GT1 != '.' && GT1 != '0') GT.emplace_back(WITH_VAR); else GT.emplace_back(WITHOUT_VAR);
            if(GT2 != '.' && GT2 != '0') GT.emplace_back(WITH_VAR); else GT.emplace_back(WITHOUT_VAR);
        }
    }
    void print(){
        fprintf(stderr, "chrID %d, POS %d, ALT %s , REF %s \n", chrID, POS, ALT.c_str(), REF.c_str());
    }
    int chrID;
    int POS;

    std::string ALT;
    std::string REF;
    std::vector<char> GT;
};

struct VCF_Reader{
    //flags
    std::map<std::string, uint32_t> filter_map;
    bool is_header_line;
    //buffs
    char *temp;//10M
    char *analysis_line;//10M
    gzFile vcf_file;
    std::vector<std::string> item_value;
    //results
    std::vector<GT_analysis> region_GT;
    int region_GT_true_size = 0;
    int var_st_idx = 0;
    int var_ed_idx = 0;

    int max_load_sample = 0;
    int sample_numer = 0;
    int max_loading_line;

    int total_line_numer;

    bool reach_EOF;

    void init(char * vcf_fn, int max_load_sample_, int max_loading_line_){
        is_header_line = true;
        temp = new char[MAX_LINE_LENGTH];//10M
        analysis_line = new char[MAX_LINE_LENGTH];//10M
        vcf_file = gzopen(vcf_fn, "rb");
        max_load_sample = max_load_sample_;
        if(max_load_sample == -1)
            max_load_sample = 200000000;
        max_loading_line = max_loading_line_;
        total_line_numer = 0;
        reach_EOF = false;
    }

    void destroy(){
        free(temp);
        free(analysis_line);
        gzclose(vcf_file);
    }

    int get_sample_number(){
        return sample_numer;
    }

    int get_chr_ID(){
        if(region_GT.empty()){
            return -1;
        }
        return region_GT[0].chrID;
    }

    //return pos
    void load_line(int & chr_ID, int & pos);
    //return true load var number
    void load_vcf_in_region(int chr_ID, int region_st, int region_ed);
    bool is_empty(){ return (var_ed_idx <= var_st_idx); }
};

struct Modify_info{
    Modify_info(int original_pos_){
        original_pos = original_pos_;
        is_modified = false;
    }
    Modify_info(){
        original_pos = -1;
        is_modified = false;
    }
    char is_modified;
    int original_pos;
    void print(){
        fprintf(stderr, " %s @ pos %d\n", is_modified?"TRUE":"False", original_pos);
    }
};
struct Hap_seq_handler{
    char * ref_pointer;
    int true_region_load_len;
    int load_region_st;
    int load_region_ed;

    int ALT_store_offset;
    int kmer_size;
    int chr_ID;

    void set_ref( int chr_ID_, int load_region_st_, int load_region_ed_, int step_length_, int kmer_size_){
        load_region_st = load_region_st_;
        load_region_ed = load_region_ed_;
        kmer_size = kmer_size_;
        with_wrong_check = false;
        chr_ID = chr_ID_;
        get_reference();
    }

    //results
    //haplotype string
    std::string HAP_str;
    //std::string modify_info;
    std::vector<Modify_info> modify_info;
    int var_number;   		//flags
    bool with_wrong_check;	//flags

    faidx_t * fai;//index file for fasta
    void init(char * reference_fn){
        //open and load reference
        fai = fai_load(reference_fn);
        ref_pointer = NULL;
    }

    void destroy(){
        fai_destroy(fai);
    }

    void get_reference();

    void modified_info_init(){
        modify_info.clear();
        int modify_info_size = HAP_str.size();
        for(int i = 0; i < modify_info_size; i++){
            modify_info.emplace_back(load_region_st + i);
        }
    }
    void get_haplotype_seq(std::vector<GT_analysis> &region_GT, int h, int region_GT_st, int region_GT_ed);

    std::string & store_hap_string(){ return HAP_str; }

    void store_alt_string(std::unordered_set<std::string> & to);

private:
};

//128 bit each block, 20M * 16 = 320M
struct window_block_info{
    window_block_info(uint8_t chrID_, uint32_t region_st_, uint32_t region_ed_, uint64_t global_offset_){
        global_offset = global_offset_;
        region_st = region_st_;
        region_length = region_ed_ - region_st_;
        chrID = chrID_;
        flag = 0;
    }

    void set(uint8_t chrID_, uint32_t region_st_, uint32_t region_ed_, uint64_t global_offset_){
        global_offset = global_offset_;
        region_st = region_st_;
        region_length = region_ed_ - region_st_;
        chrID = chrID_;
        flag = 0;
    }
    window_block_info(const window_block_info & from){
        global_offset = from.global_offset;
        region_st = from.region_st;
        region_length = from.region_length;
        chrID = from.chrID;
        flag = from.flag;
    }
    window_block_info(){
        global_offset = 0;
        region_st = 0;
        region_length = 0;
        chrID = 0;
        flag = 0;
    }

    void show_wb(FILE * out){
        fprintf(out, "@ window block, global_offset %ld, region_st %d, region_length %d, chrID %d, flag %d\n", global_offset, region_st,region_length, chrID, flag);
    }

    bool is_SV(){
        return ((flag & 0x1) == 1);
    }

    uint64_t global_offset; //window block DATA in hap_database
    uint32_t region_st; //window block POS in reference
    uint16_t region_length;//window block length
    uint8_t chrID;
    uint8_t flag; //the bit7:( == 0) when it is from SNP/INDEL; ( == 1) when it is from SV;
};

#define WINDOW_SIZE 300
#define WINDOW_OVERLAP 150

struct HAP_string_single_VCF_builder_PARA{
    char* ref_fn;
    char * vcf_fn;
    int window_size;
    int overlap;
    int step_length;
    char * wb_dump_fn;

    int get_option(int argc, char *argv[]){
        options_list l;
        l.show_command(stderr, argc + 1, argv - 1);
        l.add_title_string("\n");
        l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  panG_SNP_INDEL  [Options] [ref.fa] [vcf_file.vcf(.gz)] [wb_dump_fn.bin] \n");
        l.add_title_string("  Basic:   \n");
        l.add_title_string("    [wb_dump_fn.bin]  FILE   the output file\n");
        //thread number
        l.add_option("window_size", 	'w', "WINDOW_SIZE of haplotype window", true, WINDOW_SIZE ); l.set_arg_pointer_back((void *)&window_size);
        l.add_option("overlap",  		'o', "overlap of haplotype window", true, WINDOW_OVERLAP); l.set_arg_pointer_back((void *)&overlap);
        if(l.default_option_handler(argc, argv)) {
            l.show_c_value(stderr);
            return 1;
        }
        l.show_c_value(stderr);

        if (argc - optind < 3)
            return l.output_usage();

        ref_fn = argv[optind];
        vcf_fn = argv[optind + 1];
        wb_dump_fn = argv[optind + 2];

        step_length = window_size - overlap;
        return 0;
    }
};

struct HAP_string_single_VCF_builder{
private:
    /***
     * data or buffs to build a windows block
     */
    std::unordered_set<std::string> unique_alt_string_set;//unused by now
    std::unordered_set<std::string> unique_hap_seq_set;
    std::string ref_hap_seq_string;//the hap seq that without VARs
    //compact of haplotype
    std::vector<std::vector<uint8_t>> hap_seq_list;
    //buff for store each window block
    std::vector<uint8_t> window_block_buff;

    void store_match(std::vector<uint8_t> & hap_seq, int length);
    void store_SNP(std::vector<uint8_t> & hap_seq, char alt_char);
    void store_del(std::vector<uint8_t> & hap_seq, int length);
    void store_ins(std::vector<uint8_t> & hap_seq, int length, std::string & alt);

    /**
     * input: a vcf reads, and all variance in a window
     * output:
     * (1) an minimizer list for alt-string (todo::)
     * (2) a compact data struct for hap-seq
     * * ***/
    std::vector<uint8_t> & build_alt_string_idx(VCF_Reader &vcf_reader, Hap_seq_handler &hh, int load_region_st_, int load_region_ed_, int step_length_, int kmer_size_);
private:
    /***
     * data or buffs to build all windows blocks using a VCF file
     */
    //used to store all hap_string
    std::vector<uint8_t> final_window_block_buff;
    std::vector<window_block_info> final_window_block_info;
    void dump_final_window_block( char * wb_dump_fn);
public:
    /*****
     * Main function of the single vcf file builder, it using (SNP and INDEL) vcf file, reference file as input, and dump compact haplotype data into a file name : wb_dump_fn
     */
    void build_hap( char* reference_fn, char * vcf_fn,int step_length, int overlap, char * wb_dump_fn){
        build_hap_string_run_core(reference_fn, vcf_fn , step_length, overlap, 22);
        //dump the final windows block to files
        if(wb_dump_fn != NULL){ dump_final_window_block(wb_dump_fn); }
    }
    void build_hap_string_run_core( char* reference_fn,  char * vcf_fn,int step_length, int overlap, int kmer_size){
        //used to read vcf files
        VCF_Reader vcf_reader;
        vcf_reader.init(vcf_fn, -1, -1);
        bool load_region_already_set = false;
        //haplotype string generator
        Hap_seq_handler hh;
        hh.init(reference_fn);
        uint64_t global_offset = 0;
        //init set chrID and
        for(int chrID = -1, load_region_st = -1; load_region_st < MAX_int32t; load_region_st += step_length){
            int load_region_ed = load_region_st + step_length + overlap; // the last kmer: P >= 999 && P < 1021
            vcf_reader.load_vcf_in_region(chrID, load_region_st, load_region_ed);
            //auto set chrID and load_region_st
            if(load_region_already_set == false){
                if(vcf_reader.reach_EOF) break;
                load_region_already_set = true;
                chrID = vcf_reader.get_chr_ID();
                load_region_st = ((vcf_reader.region_GT[0].POS/step_length ) * step_length - step_length) + (1);
                continue;
            }else if(chrID != vcf_reader.get_chr_ID()){
                load_region_already_set = false;
                chrID = vcf_reader.get_chr_ID();
                continue;
            }

            std::vector<uint8_t> &window_block = build_alt_string_idx(vcf_reader, hh, load_region_st, load_region_ed, step_length, kmer_size);
            final_window_block_info.emplace_back(chrID, load_region_st, load_region_ed, global_offset);
            //final_window_block_info.back().show_wb(stderr);	fprintf(stderr, "global_offset @ %ld \n", global_offset);
            global_offset += window_block.size();
            final_window_block_buff.insert(final_window_block_buff.end(), window_block.begin(), window_block.end());
        }
    }
};

struct Hap_seq_address{
    uint8_t *p;
    uint64_t len;
    void clear(){
        p = NULL;
        len = 0;
    }
};

struct Hap_seq_modify_info{
    Hap_seq_modify_info(int original_pos_){
        original_pos = original_pos_;
        is_modified = false;
    }
    Hap_seq_modify_info(){
        original_pos = -1;
        is_modified = false;
    }
    char is_modified;
    int original_pos;
    void print(){
        fprintf(stderr, " %s @ pos %d\n", is_modified?"TRUE":"False", original_pos);
    }
};

// X：[30 bit hash value] + [1 bit REF or ALT] + [1 bit Reverse or Forward]
// Y + Z：[25bit window ID @ MAX 32M space] + [19bit: hap ID @ MAX 0.5M space] + [20bit: hap offset @ MAX 1M]
struct mm_96_t{
    mm_96_t(uint32_t x_, uint32_t y_,uint32_t z_){
        x = x_;
        y = y_;
        z = z_;
    }
    uint32_t x, y, z;
};

struct Minimizer_index_builder{

private:
        std::unordered_set<std::string> alt_string_uniq_set;
        std::unordered_set<std::string> alt_string_uniq_set_old;
        int window_size_minimizer;
        int kmer_size_minimizer;
        int kmer_size_alt_string;

        Minimizer_generater mg;
public:
        void new_window_block_clear(int window_size_minimizer_, int kmer_size_minimizer_, int kmer_size_alt_string_){
            std::swap(alt_string_uniq_set, alt_string_uniq_set_old);
            alt_string_uniq_set.clear();
            window_size_minimizer = window_size_minimizer_;
            kmer_size_minimizer = kmer_size_minimizer_;
            kmer_size_alt_string = kmer_size_alt_string_;
            mg.init(window_size_minimizer, kmer_size_minimizer, false);
        }
private:
        std::string alt_string;
        std::string alt_string_with_pos;
public:
        void minimizer_generater_hap_seq(std::vector<Hap_seq_modify_info> &modify_info, std::vector<char> &hap_seq,
                uint64_t window_ID, uint64_t hap_ID, std::vector<mm128_t> &mm_v, bool is_from_SV);

        void minimizer_generater_unitig(std::string &UNITIG_str, uint64_t unitig_ID, std::vector<mm128_t> &mm_v);

        void mm_sort(std::vector<mm128_t> &mm_v);

private:

public:
        void building_mm_idx_and_dump(FILE * mm_v_f, FILE * mm_idx_f, std::vector<mm128_t> &sort_mm_v, int kmer_size_minimizer);
    };

struct HAP_VAR_ITEM{
    HAP_VAR_ITEM(uint32_t ref_pos_, std::string &REF_, std::string &ALT_){
        ref_pos = ref_pos_;
        std::swap(REF,REF_);
        std::swap(ALT,ALT_);
    }
    uint32_t ref_pos = 0;
    std::string REF;
    std::string ALT;

    int get_len(){
        return ALT.size() - REF.size();
    }
    void print(FILE * out){
        fprintf(out, "%s_%s_%d\t", REF.c_str(), ALT.c_str(), ref_pos);
    }

    bool try_to_combine(HAP_VAR_ITEM & B){
        //skip combine SNPs
        //if(get_len() == 0 || B.get_len() == 0)
        //	return false;
        //if(ref_pos != B.ref_pos)return false;
        if(false){
            fprintf(stderr, "Before: ");
            print(stderr);
            B.print(stderr);
            fprintf(stderr, "\n ");
        }
        int pos_diff = B.ref_pos - ref_pos; pos_diff = MAX(0, pos_diff);
        REF.insert(REF.begin() + 1 + pos_diff, B.REF.begin() + 1, B.REF.end());
        ALT.insert(ALT.end(), B.ALT.begin() + 1, B.ALT.end());

        //if(ALT.size() > 1){// && B.ALT.size() > 1){ //both INS
        //	ALT.insert(ALT.end(), B.ALT.begin() + 1, B.ALT.end());
        //}else{
        //	ALT.insert(ALT.begin() + 1, B.ALT.begin() + 1, B.ALT.end());
        //}

        //remove same string in the tail
        while(true){
            if(REF.size() > 1 && ALT.size() > 1 && REF.back() == ALT.back()){
                REF.pop_back();
                ALT.pop_back();
            }else
                break;
        }
        if(false){
            fprintf(stderr, "After: ");
            print(stderr);
            fprintf(stderr, "\n ");
        }
        return true;
    }

    bool is_no_var(){
        return (REF == ALT);
    }
};

struct Window_block_handler{
private:
    //
    std::vector<char> rst_haplotype_BUFF;
    /***
     * unfold a hap seq and generate a haplotype string
     * input:: (1) a pointer to the bin hap seq
     * input:: (2) length of the bin hap seq
     * input:: (3) the string of reference
     *
     * Output: the haplotype sequence
     */

    std::vector<Hap_seq_modify_info> modify_info;
public:
    std::vector<char> & get_string(){ return rst_haplotype_BUFF; }
    std::vector<Hap_seq_modify_info> & get_modify_info(){ return modify_info; }
    //debug code:: print modify inf0
    void printf_modify_info(FILE* out){
        for(auto & h : rst_haplotype_BUFF){ fprintf(out, "%c", h); }
        fprintf(out, "\n");
        for(auto & mi : modify_info){ fprintf(out, "%c", (mi.is_modified == true)?'1':'0'); }
        fprintf(out, "\n");
    }

    std::vector<HAP_VAR_ITEM> var_l;
    std::vector<char> & unfold_hap_seq(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, bool generate_modify_string, bool generate_var_list, int kmer_size_alt_string);

    std::vector<char> & unfold_hap_seq_SV(uint8_t * bin_hap_seq, uint64_t hap_seq_len, std::string & ref, window_block_info & c_wb, bool generate_modify_string, uint32_t kmer_size_alt_string);

    uint8_t * window_block_p_old = NULL;
    std::vector<Hap_seq_address> block_all_hap_list;
    void unfold_window_block_all(uint8_t * window_block_p);
    /***
     * unfold one hap seq in a window block
     * input:: (1) a pointer to the window block
     * input:: (2) ID of the hap seq
     *
     * Output: block_i_rst_address of the haplotype sequence
     */
    Hap_seq_address block_i_rst_address;
    void unfold_window_block_i(uint8_t * window_block_p, uint64_t hap_i);
};

struct Simple_ref_handler1{
    std::vector<std::string> all_ref_data_list;
    std::vector<std::string> all_ref_name_list;
    //window block part loader
    void load_reference_from_file(const char * ref_fn){
        faidx_t * fai = fai_load(ref_fn);
        uint32_t ref_n = faidx_nseq(fai);
        //debug code::
        int debug_load_ref = -1;
//		if(true){ debug_load_ref = 9; }//debug code
        for(uint32_t i = 0; i < ref_n; i++){
            all_ref_name_list.emplace_back(faidx_iseq(fai, i));
            int true_region_load_len = 0;
            //debug code
            if(debug_load_ref != -1 && debug_load_ref != (int)i){
                all_ref_data_list.emplace_back();
                continue;
            }
            char * ref_pointer = fai_fetch(fai, faidx_iseq(fai, i), &true_region_load_len);
            all_ref_data_list.emplace_back(ref_pointer);
            if(ref_pointer != NULL) free(ref_pointer);
        }
        fai_destroy(fai);
        fprintf(stderr, "all_ref_name_list.size() %ld\n", all_ref_name_list.size());
    }

    //window block part loader
    void load_ref_from_buff(int chr_ID, uint32_t st_pos, int length, std::string & rst){
        rst.clear();
        if(chr_ID < 0 || chr_ID >= (int)all_ref_data_list.size())	return;
//		for(int i = 0; i < length; i++){
//			fprintf(stderr, "%c", all_ref_data_list[chr_ID][st_pos - 1 + i]);
//		}
        rst.insert(rst.begin(), all_ref_data_list[chr_ID].begin() + st_pos - 1, all_ref_data_list[chr_ID].begin() + st_pos - 1 + length);
        return ;
    }

    std::string & get_chr_name(int chr_ID){
        return all_ref_name_list[chr_ID];
    }
    uint64_t get_chr_length(int chr_ID){
        return all_ref_data_list[chr_ID].size();
    }
    uint64_t get_chr_N(){
        return all_ref_data_list.size();
    }

};
struct GAPRecord{
    int tid;
    uint32_t start;
    uint32_t end;

    bool operator<(const GAPRecord& other) const {
    if (tid == other.tid) {
        return start < other.start;
    }
    return tid < other.tid;
    }

};

struct Simple_ref_handler{
public:
    std::vector<GAPRecord> gap_records;
    void load_gap_file(const char * path_name) {
        char full_fn[1024] = {0};


        strcpy(full_fn, path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        fprintf(stderr, "gap file %s\n", full_fn);
        strcat(full_fn, "/ref.gap.bed");

        std::ifstream infile(full_fn);
        std::string line;


        while (std::getline(infile, line)) {
            GAPRecord record;
            std::istringstream iss(line);
            iss >> record.tid >> record.start >> record.end;
            gap_records.push_back(record);
        }

        infile.close();

        // 排序
        std::sort(gap_records.begin(), gap_records.end());

    }



private:
    //flag:
    bool ref_store_in_binary;
    //reference stored in ACGT
    std::vector<std::string> all_ref_data_list_ACGT;


    //reference store in binary
    struct ref_info{
        ref_info(){
            ref_length = 0;
            ref_offset = 0;
        }
        ref_info(std::string &s){
            std::vector<std::string> item_value;
            char * temp = (char *)xcalloc(1024,1);
            split_string(item_value, temp, s.c_str(), " ");
            name = item_value[0];
            ref_length = atol(item_value[1].c_str());
            ref_offset = atol(item_value[2].c_str());
            if(temp != NULL){
                free(temp); temp = NULL;
            }
        }
        uint32_t ref_length;
        uint64_t ref_offset;
        std::string name;
        void dump(FILE * out){
            fprintf(out, "%s %d %ld\n", name.c_str(), ref_length, ref_offset);
        }
    };
    std::vector<ref_info> all_ref_data_list;
    uint8_t *bin_ref_data;

    inline uint8_t char2bin(char c){
        switch(c){
        case 'A': case 'a': return 0;  break;
        case 'C': case 'c': return 1;  break;
        case 'G': case 'g': return 2;  break;
        case 'T': case 't': return 3;  break;
        case 'N': case 'n': return 2;  break;
        }
        return 0;
    }
public:
    //dump the reference info into bin_ref_info_f, dump the bin ref into bin_ref_f
    void dump_bin_ref(FILE * bin_ref_f, FILE * bin_ref_info_f){
        xassert(ref_store_in_binary == false, "");
        uint64_t total_uint8_t_size = all_ref_data_list.back().ref_offset + (all_ref_data_list.back().ref_length + 3)/4;
        bin_ref_data = (uint8_t *)xcalloc(total_uint8_t_size, 1);
        for(uint32_t chr_ID = 0;chr_ID < all_ref_data_list_ACGT.size(); chr_ID ++){
            uint64_t ref_global_offset = all_ref_data_list[chr_ID].ref_offset;
            std::string& ref_string = all_ref_data_list_ACGT[chr_ID];
            for(uint32_t offset = 0; offset < ref_string.size(); offset++){
                uint8_t bin_char = char2bin(ref_string[offset]);
                bin_ref_data[((offset) >> 2) + ref_global_offset] |= bin_char << ((3 - ((offset) & 0X3)) << 1);
            }
        }

        //dump files
        for(uint32_t chr_ID = 0;chr_ID < all_ref_data_list_ACGT.size(); chr_ID ++){
            all_ref_data_list[chr_ID].dump(bin_ref_info_f);
        }
        fwrite(bin_ref_data, 1, total_uint8_t_size, bin_ref_f);
        if(bin_ref_data != NULL){
            free(bin_ref_data); bin_ref_data = NULL;
        }
    }

    void load_bin_ref(const char * path_name){
        ref_store_in_binary = true;
        std::vector<std::string> chr_if_str;
        //load info file
        char full_fn[1024] = {0};
        strcpy(full_fn, path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        fprintf(stderr, "refinfo %s\n", full_fn);
        strcat(full_fn, "/ref_info.bin");
        load_strings_from_file(full_fn, chr_if_str, MAX_uint64_t);
        all_ref_data_list.clear();
        for(auto & s: chr_if_str){
            all_ref_data_list.emplace_back(s);
        }
        //load data file:
        strcpy(full_fn, path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        strcat(full_fn, "/ref.bin");
        FILE * ref_data_f = xopen(full_fn, "rb");
        uint64_t total_uint8_t_size = all_ref_data_list.back().ref_offset + (all_ref_data_list.back().ref_length + 3)/4;
        bin_ref_data = (uint8_t *) xmalloc(total_uint8_t_size);
        xread(bin_ref_data, 1, total_uint8_t_size, ref_data_f);
        fclose(ref_data_f);
    }

    //window block part loader, simple load reference from ref.fa
    void load_reference_from_file(const char * ref_fn){
        ref_store_in_binary = false;
        faidx_t * fai = fai_load(ref_fn);
        uint32_t ref_n = faidx_nseq(fai);
        uint64_t total_uint8_t_size = 0;
//		if(true){ debug_load_ref = 9; }//debug code
        for(uint32_t chr_ID = 0; chr_ID < ref_n; chr_ID++){
            int true_region_load_len = 0;
            char * ref_pointer = fai_fetch(fai, faidx_iseq(fai, chr_ID), &true_region_load_len);
            all_ref_data_list_ACGT.emplace_back(ref_pointer);

            std::string& ref_string = all_ref_data_list_ACGT[chr_ID];
            all_ref_data_list.emplace_back();
            all_ref_data_list.back().name.append(faidx_iseq(fai, chr_ID));
            all_ref_data_list.back().ref_length = ref_string.size();
            all_ref_data_list.back().ref_offset = total_uint8_t_size;
            total_uint8_t_size += (ref_string.size() + 3)/4;

            if(ref_pointer != NULL) free(ref_pointer);
        }
        fai_destroy(fai);
    }

//	void load_ref_from_buff_global_offset(uint64_t global_offset, int length, std::string & rst, deBGA_INDEX * dbi){
//		int chr_ID, POS;
//		dbi->get_chr_ID_POS(global_offset, chr_ID, POS);
//		load_ref_from_buff(chr_ID, POS, length, rst);
//	}

    //window block part loader
    void load_ref_from_buff(int chr_ID, uint32_t st_pos, int length, std::string & rst){
        rst.clear();
        if(chr_ID < 0 || chr_ID >= (int)all_ref_data_list.size())	return;
        if(ref_store_in_binary){
            int32_t max_load = (all_ref_data_list[chr_ID].ref_length - st_pos + 1);
            max_load = MIN(max_load, length);
            if(max_load < 0)
                return;
            rst.resize(max_load);
            uint64_t c_ref_global_offset = all_ref_data_list[chr_ID].ref_offset;
            for(int32_t i = 0; i < max_load; i++){
                uint8_t b_c = ((bin_ref_data[((st_pos + i - 1) >> 2) + c_ref_global_offset]) >> (((3 - ((st_pos + i - 1) & 0X3)) << 1))) & 0x3;
                rst[i] = "ACGT"[b_c];
            }
        }else{
            rst.insert(rst.begin(), all_ref_data_list_ACGT[chr_ID].begin() + st_pos - 1, all_ref_data_list_ACGT[chr_ID].begin() + st_pos - 1 + length);
        }
    }

    //window block part loader
    void load_ref_bin_from_buff(int chr_ID, uint32_t st_pos, int length, std::vector<uint8_t> & rst){
        rst.clear();
        if(chr_ID < 0 || chr_ID >= (int)all_ref_data_list.size())	return;
        if(ref_store_in_binary){
            int32_t max_load = (all_ref_data_list[chr_ID].ref_length - st_pos + 1);
            max_load = MIN(max_load, length);
            if(max_load < 0)
                return;
            rst.resize(max_load);
            uint64_t c_ref_global_offset = all_ref_data_list[chr_ID].ref_offset;
            for(int32_t i = 0; i < max_load; i++){
                uint8_t b_c = ((bin_ref_data[((st_pos + i - 1) >> 2) + c_ref_global_offset]) >> (((3 - ((st_pos + i - 1) & 0X3)) << 1))) & 0x3;
                rst[i] = b_c;
            }
        }else{
            xassert(0, "Run only when ref stored in binary!");
        }
    }

    std::string & get_chr_name(int chr_ID){
        return all_ref_data_list[chr_ID].name;
    }
    uint64_t get_chr_length(int chr_ID){
        return all_ref_data_list[chr_ID].ref_length;
    }
    uint64_t get_chr_N(){
        return all_ref_data_list.size();
    }
};

#define MINIMIZER_LEN 22
#define MINIMIZER_WINDOW_SIZE 5

struct Hap_seq_minimizer_idx_global_builder_PARA{

    int minimizer_len;//kmer size when building minimizers
    int minimizer_window_size;

    char * ref_fn;
    char * wb_fn_list_fn;
    char * unitig_fn;
    std::vector<std::string> wb_fn_list;
    char * idx_dir;

    int get_option(int argc, char *argv[]){
        options_list l;
        l.show_command(stderr, argc + 1, argv - 1);
        l.add_title_string("\n");
        l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  index_mm  [Options] [ref.fa] [wb_dump_fn_list.txt] [unitig.txt] [index_dir]\n");
        l.add_title_string("  Basic:   \n");
        l.add_title_string("    [wb_dump_fn_list.txt]  FILE   The file name list of the [wb_dump_fn.bin] files\n");
        l.add_title_string("    [unitig.txt]  FILE   The UNITIG file.\n");
        l.add_title_string("    [index_dir]  The dumped index directory.\n");
        //thread number
        l.add_option("minimizer_len", 	'm', "length of the minimizer", true, MINIMIZER_LEN ); l.set_arg_pointer_back((void *)&minimizer_len);
        l.add_option("minimizer_window_size", 	'w', "window_size of minimizer", true, MINIMIZER_WINDOW_SIZE ); l.set_arg_pointer_back((void *)&minimizer_window_size);
        if(l.default_option_handler(argc, argv)) {
            l.show_c_value(stderr);
            return 1;
        }
        l.show_c_value(stderr);

        if (argc - optind < 3)
            return l.output_usage();

        ref_fn = argv[optind];
        wb_fn_list_fn = argv[optind + 1];
        unitig_fn = argv[optind + 2];
        load_strings_from_file(wb_fn_list_fn, wb_fn_list, MAX_uint64_t);
        idx_dir = argv[optind + 3];
        return 0;
    }
};


struct Hap_files_merger_PARA{
    char * ref_fn;
    char * wb_fn_list_fn;
    std::vector<std::string> wb_fn_list;
    char * idx_dir;

    int get_option(int argc, char *argv[]){
        options_list l;
        l.show_command(stderr, argc + 1, argv - 1);
        l.add_title_string("\n");
        l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  hap_files_merge  [Options] [ref.fa] [wb_dump_fn_list.txt] [index_dir]\n");
        l.add_title_string("  Basic:   \n");
        l.add_title_string("    [wb_dump_fn_list.txt]  FILE   The file name list of the [wb_dump_fn.bin] files\n");
        l.add_title_string("    [index_dir]  The dumped index directory.\n");
        if(l.default_option_handler(argc, argv)) {
            l.show_c_value(stderr);
            return 1;
        }
        l.show_c_value(stderr);

        if (argc - optind < 3)
            return l.output_usage();

        ref_fn = argv[optind];
        wb_fn_list_fn = argv[optind + 1];
        load_strings_from_file(wb_fn_list_fn, wb_fn_list, MAX_uint64_t);
        idx_dir = argv[optind + 2];
        return 0;
    }
};

struct Hap_seq_idx_global_builder{
    struct Single_file_loader{
        uint64_t window_block_size;
        window_block_info * wb_info;
        uint64_t window_block_data_size;
        uint8_t * wb_data;

        void destroy(){
            if(wb_info != NULL) {free(wb_info); wb_info = NULL;}
            if(wb_data != NULL) {free(wb_data); wb_data = NULL;}
        }
        void load_data(const char * wb_dump_fn){
            window_block_size = 0;
            wb_info = NULL;
            window_block_data_size = 0;
            wb_data = NULL;
            if(wb_dump_fn != NULL){
                FILE * dump_f = xopen(wb_dump_fn, "rb");//equal to : fopen(wb_dump_fn, wb_dump_fn);
                window_block_size = CPP_vector_load_bin(dump_f, (void **)&wb_info, sizeof(window_block_info));
                window_block_data_size = CPP_vector_load_bin(dump_f, (void **)&wb_data, 1);
                fclose(dump_f);
            }
        }
    };
    Single_file_loader sfl;

    //reference list, when building index ,just simple loading all the reference into memorys
    Simple_ref_handler ref;
private:
    std::vector<window_block_info> wb_info;

    FILE* bin_data_open(const char * idx_dir, const char * fn, const char *mode){
        std::string fn_string;
        fn_string.append(idx_dir);
        fn_string.append(fn);
        return xopen(fn_string.c_str(), mode);//fopen64 in the core function
    }

    void ALT_string_index_building(std::vector<std::string> & wb_fn_list, std::vector<mm128_t> *mm_v, std::vector<uint64_t> &chr_bg_wb_ID, FILE *wb_data_f,
            Minimizer_index_builder *mib, int kmer_size_alt_string, int kmer_size_minimizer, int window_size_minimizer){
        //init:: window block handler
        //std::
        chr_bg_wb_ID.clear();
        Window_block_handler wbh;
        uint64_t wb_data_global_offset = 0;
        wb_info.clear();
        std::string c_ref;
        int chr_ID = -1;
        for(auto & wb_fn : wb_fn_list){
            fprintf(stderr, "Current processing: %s\n", wb_fn.c_str());

            xassert(chr_ID != -2, "The data file from SV must be at the end of file lists, only one file (from SV) is accepted");
            sfl.load_data(wb_fn.c_str());
            //S1: adding blank window block in the begin of each chromosome
            //check::
            chr_ID++;
            chr_bg_wb_ID.emplace_back(wb_info.size());

            bool file_from_SV = (sfl.wb_info[0].flag & 0x1);

            if(file_from_SV == false){
                xassert((uint8_t)chr_ID == sfl.wb_info[0].chrID, "Data in wb_fn_list must begin with chr1, and followed by chr2, chr3 .....");
            }else{
                chr_ID = -2;//when the file is SV-based wb, skip check and reset the chr_ID
            }

            uint8_t chr_ID = sfl.wb_info[0].chrID;
            uint16_t step_length = sfl.wb_info[0].region_length/2;
            uint16_t wb_length = sfl.wb_info[0].region_length;

            uint64_t new_wb_in_the_begin = 0;
            if(file_from_SV == false){
                new_wb_in_the_begin = (sfl.wb_info[0].region_st - 1) / (step_length);
                for(uint64_t i = 0; i < new_wb_in_the_begin; i++){
                    wb_info.emplace_back(chr_ID, step_length * i + 1, step_length * i + wb_length + 1, wb_data_global_offset);
                    wb_data_global_offset += 2;
                }
                //write blank wb data:
                uint8_t * blank_wb_data = (uint8_t * )xcalloc(new_wb_in_the_begin, 2);
                fwrite(blank_wb_data,1,new_wb_in_the_begin * 2, wb_data_f);
                free(blank_wb_data);
            }

            //S2: adding true windows blocks
            for(uint64_t i = 0; i < sfl.window_block_size; i++){
                wb_info.emplace_back(sfl.wb_info[i]);
                wb_info.back().global_offset += wb_data_global_offset;
            }
            wb_data_global_offset += sfl.window_block_data_size;
            //dump the windows block data into file
            fwrite(sfl.wb_data,1,sfl.window_block_data_size, wb_data_f);

            if(file_from_SV == false){
                //S3:adding blank window block in the end of each chromosome
                uint64_t new_wb_in_the_end = (ref. get_chr_length(chr_ID) - wb_info.back().region_st + 1000) / (step_length);
                int back_wb_st = wb_info.back().region_st + step_length;
                for(uint64_t i = 0; i < new_wb_in_the_end; i++){
                    wb_info.emplace_back(chr_ID, step_length * i + back_wb_st, step_length * i + wb_length + back_wb_st, wb_data_global_offset);
                    wb_data_global_offset += 2;
                }
                //write blank wb data:
                uint8_t * blank_wb_data = (uint8_t * )xcalloc(new_wb_in_the_end, 2);
                fwrite(blank_wb_data,1,new_wb_in_the_end * 2, wb_data_f);
                free(blank_wb_data);
            }

            //building minimizer index
            if(mib != NULL){
                uint64_t window_ID = chr_bg_wb_ID.back() + new_wb_in_the_begin;
                //S5：unfold all haplotypes and building all minimizer index
                for(uint64_t i = 0; i < sfl.window_block_size; i++){
                    window_block_info & c_wb = sfl.wb_info[i];
                    window_ID ++;
                    int hap_ID = 0;
                    if(false)//debug code
                        c_wb.show_wb(stdout);//debug code
                    //ref info:
                    //unfold all window block and build minimizer index for all alt-strings
                    wbh.unfold_window_block_all(sfl.wb_data + c_wb.global_offset);
                    if(wbh.block_all_hap_list.empty()){ continue; }
                    //load reference:
                    ref.load_ref_from_buff(c_wb.chrID, c_wb.region_st, c_wb.region_length, c_ref);
                    if(false)//debug code
                        fprintf(stdout, "\nRef hap: %s\n", c_ref.c_str());

                    std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
                    //clear minimizer builder:
                    mib->new_window_block_clear(window_size_minimizer, kmer_size_minimizer, kmer_size_alt_string);
    //				if(c_wb.region_st == 75523800){ fprintf(stderr, "DEBUG~~~~~"); } // debug code
                    for(auto & hap:hap_l){//load all haplotypes
                        //unfold_hap_seq_SV
                        bool wb_store_SV = c_wb.flag & 0x1;
                        std::vector<char> & hap_string = (wb_store_SV)?wbh.unfold_hap_seq_SV(hap.p, hap.len, c_ref, c_wb, true, kmer_size_alt_string):wbh.unfold_hap_seq(hap.p, hap.len, c_ref, true, false, kmer_size_alt_string);
                        //for(char A: hap_string){ if(A < 'A' || A > 'T') fprintf(stderr, ""); } //debug code
                        //hap_string.emplace_back(0);
                        hap_ID ++;
                        if(false){//debug code
                            //fprintf(stdout, "@1 HAP_ID_%d_OFFSET_0\n", ++hap_idx);
                            fprintf(stdout, "%s\n+\n", &(hap_string[0]));
                            for(uint32_t q_i = 0; q_i < hap_string.size() - 1; q_i++){ fprintf(stdout, "B"); }
                            fprintf(stdout, "\n");
                            if(true){//debug code
                                wbh.printf_modify_info(stdout);
                            }
                        }
                        //minimizer generator: core code
                        {
                            std::vector<Hap_seq_modify_info> & modify_info = wbh.get_modify_info();
                            //generate minimizer:
                            mib->minimizer_generater_hap_seq(modify_info, hap_string, window_ID, hap_ID, *mm_v, file_from_SV);
                        }
                    }
                }
            }

            sfl.destroy();//free memory
        }
        //store the total wb_info size
        chr_bg_wb_ID.emplace_back(wb_info.size());
    }

public:
    void merge_window_block_and_build_mm_index(char * ref_fn, std::vector<std::string> & wb_fn_list, const char * unitig_fn, const char * idx_dir,
            int kmer_size_minimizer, int window_size_minimizer){

        FILE * wb_data_f = bin_data_open(idx_dir, "/wb_data.bin", "wb");
        FILE * wb_info_f = bin_data_open(idx_dir, "/wb_info.bin", "wb");

        FILE * alt_mm_f = bin_data_open(idx_dir, "/alt_mm.bin", "wb");
        FILE * alt_mm_idx_f = bin_data_open(idx_dir, "/alt_mm_idx.bin", "wb");

        FILE * ref_mm_f = bin_data_open(idx_dir, "/ref_mm.bin", "wb");
        FILE * ref_mm_idx_f = bin_data_open(idx_dir, "/ref_mm_idx.bin", "wb");

        //dump bin reference
        FILE * bin_ref_f = bin_data_open(idx_dir, "/ref.bin", "wb");
        FILE * bin_ref_info_f = bin_data_open(idx_dir, "/ref_info.bin", "wb");

        //dump bin reference
        FILE * chr_bg_wb_ID_f = bin_data_open(idx_dir, "/chr_bg_wb_ID.bin", "w");

        //para::
        int kmer_size_alt_string = kmer_size_minimizer;//kmer size when generate alt strings
        //init:: reference file
        fprintf(stderr, "load_reference_from_file %s\n", ref_fn);
        ref.load_reference_from_file(ref_fn);
        //
        std::vector<uint64_t> chr_bg_wb_ID;
        //
        //minimizer list:
        std::vector<mm128_t> mm_v;

        Minimizer_index_builder mib;

        //merge and dump files
        if(true){//debug set false
            ALT_string_index_building(wb_fn_list, &mm_v, chr_bg_wb_ID, wb_data_f, &mib, kmer_size_alt_string, kmer_size_minimizer, window_size_minimizer);
            //dump window block info vector into file
            CPP_vector_dump_bin(wb_info_f, &(wb_info[0]), wb_info.size() * sizeof(window_block_info));
            for(uint i = 0; i < chr_bg_wb_ID.size(); i++){
                fprintf(chr_bg_wb_ID_f, "%ld\n", chr_bg_wb_ID[i]);
            }
            fclose(chr_bg_wb_ID_f);
            fclose(wb_data_f);
            fclose(wb_info_f);
        }
        //generating mm index for hap seq
        if(true){//debug set false
            //sort minimizer
            mib.mm_sort(mm_v);
            //building index for sorted minimizer:
            mib.building_mm_idx_and_dump(alt_mm_f, alt_mm_idx_f, mm_v, kmer_size_minimizer);
            mm_v.clear();
            fclose(alt_mm_f);
            fclose(alt_mm_idx_f);
        }

        //generating mm index for UNITIG
        if(true){
            char *temp = new char[MAX_LINE_LENGTH];//10M
            std::ifstream UNITIG_f(unitig_fn);
            mm_v.clear();
            uint32_t unitig_ID = 0;
            std::string UNITIG_str;
            mib.new_window_block_clear(window_size_minimizer, kmer_size_minimizer, 0);
            while(true){
                UNITIG_f.getline(temp, MAX_LINE_LENGTH);
                if(UNITIG_f.eof())	break;
                UNITIG_str.clear();
                UNITIG_str.append(temp);
                mib.minimizer_generater_unitig(UNITIG_str, unitig_ID++, mm_v);
                //debug code
                if(unitig_ID%100 == 0)
                    fprintf(stderr, "unitig_ID %d, mm_v.size()  %ld, UNITIG_str.size()  %ld\n",unitig_ID, mm_v.size(), UNITIG_str.size());
            }
            //sort minimizer
            mib.mm_sort(mm_v);
            mib.building_mm_idx_and_dump(ref_mm_f, ref_mm_idx_f, mm_v, kmer_size_minimizer);
            mm_v.clear();
            fclose(ref_mm_f);
            fclose(ref_mm_idx_f);
            UNITIG_f.close();
        }

        if(true){//dump binary reference ??? todo::
            ref.dump_bin_ref(bin_ref_f, bin_ref_info_f);
            fclose(bin_ref_f);
            fclose(bin_ref_info_f);
        }
    }

    void merge_window_block(char * ref_fn, std::vector<std::string> & wb_fn_list, const char * idx_dir){
            //dump WB database
            FILE * wb_data_f = bin_data_open(idx_dir, "/wb_data.bin", "wb");
            FILE * wb_info_f = bin_data_open(idx_dir, "/wb_info.bin", "wb");
            //dump bin reference
            FILE * bin_ref_f = bin_data_open(idx_dir, "/ref.bin", "wb");
            FILE * bin_ref_info_f = bin_data_open(idx_dir, "/ref_info.bin", "wb");
            //dump bin reference
            FILE * chr_bg_wb_ID_f = bin_data_open(idx_dir, "/chr_bg_wb_ID.bin", "w");

            //init:: reference file
            fprintf(stderr, "load_reference_from_file %s\n", ref_fn);
            ref.load_reference_from_file(ref_fn);
            std::vector<uint64_t> chr_bg_wb_ID;
            //merge and dump files
            if(true){//debug set false
                ALT_string_index_building(wb_fn_list, NULL, chr_bg_wb_ID, wb_data_f, NULL, 0, 0, 0);
                //dump window block info vector into file
                CPP_vector_dump_bin(wb_info_f, &(wb_info[0]), wb_info.size() * sizeof(window_block_info));
                for(uint i = 0; i < chr_bg_wb_ID.size(); i++){
                    fprintf(chr_bg_wb_ID_f, "%ld\n", chr_bg_wb_ID[i]);
                }
                fclose(chr_bg_wb_ID_f);
                fclose(wb_data_f);
                fclose(wb_info_f);
            }

            if(true){
                ref.dump_bin_ref(bin_ref_f, bin_ref_info_f);
                fclose(bin_ref_f);
                fclose(bin_ref_info_f);
            }
        }

};

//struct get_haplotype_string_run{
//	struct run_option{
//		char * vcf_fn = NULL;
//		int64_t max_loading_line = 0;
//		char * reference_fn = NULL;
//		int max_sample_number = 0;
//		bool show_st_ed_info = 0;
//
//		int window_size;//300
//		int overlap;//150
//		int step_length = 0;
//		int kmer_size = 0;
//
//		void get_alt_string_option(int argc, char *argv[]) {
//			fprintf(stderr, "VERSION 1.0.6\n\n");
//			fprintf(stderr, "USAGE: P_NAME anno_analysis [AVR.vcf] [max_load_line] [reference.fa] [window_size] [overlap] \n\n");
//			fprintf(stderr, "[max_load_line]: set to -1 to load all vcf records\n");
//			//parameters
//			//register all analysis items
//			vcf_fn = argv[1];
//			max_loading_line = atol(argv[2]);
//			reference_fn = argv[3];
//
//			window_size = atol(argv[4]);
//			overlap = atol(argv[5]);
//
//			max_sample_number = -1;
//			//the begin read position in a annoFile
//			if(max_sample_number == -1){
//				max_sample_number = 200000000;
//			}
//
//			kmer_size = 22;
//			step_length = window_size - overlap;
//		}
//	};
//
//	int run(int argc, char *argv[]) {
//		run_option o;
//		o.get_alt_string_option(argc, argv);
//
//		//used to read vcf files
//		VCF_Reader vcf_reader;
//		vcf_reader.init(o.vcf_fn, o.max_sample_number, o.max_loading_line);
//
//		//haplotype string generator
//		Hap_seq_handler hap_handler;
//		hap_handler.init(o.reference_fn);
//
//		//count the number of haplotype string
//		std::vector<int> hap_number_count;
//
//		std::vector<int> alt_string_length_count;
//		for(int i = 0; i < o.step_length + 100; i++)
//			alt_string_length_count.emplace_back(0);
//
//		std::unordered_set<std::string> unique_alt_kmer_set;
//		std::unordered_set<std::string> unique_hap_seq_set;
//		std::string ref_hap_seq_string;//the hap seq that without VARs
//
//		//used to get split string
//		std::vector<std::string> alt_string_item_value;
//		char *temp = new char[MAX_LINE_LENGTH];//10M
//
//		int total_hap_seq_number = 0;
//		int total_alt_string_number = 0;
//
//		bool load_region_already_set = false;
//		bool show_strings = false;
//		bool show_stats = false;
//
//		//init set chrID and
//		for(int chrID = -1, load_region_st = -1; load_region_st < MAX_int32t; load_region_st += o.step_length){
//			int load_region_ed = load_region_st + o.step_length + o.overlap; // the last kmer: P >= 999 && P < 1021
//			vcf_reader.load_vcf_in_region(chrID, load_region_st, load_region_ed);
//			//auto set chrID and load_region_st
//			if(load_region_already_set == false){
//				if(vcf_reader.reach_EOF) break;
//				load_region_already_set = true;
//				chrID = vcf_reader.get_chr_ID();
//				load_region_st = (vcf_reader.region_GT[0].POS/o.step_length ) * o.step_length - o.step_length;
//				continue;
//			}else if(chrID != vcf_reader.get_chr_ID()){
//				load_region_already_set = false;
//				chrID = vcf_reader.get_chr_ID();
//				continue;
//			}
//
//			if(hap_number_count.empty()){
//				for(int i = 0; i < vcf_reader.sample_numer*2; i++)
//					hap_number_count.emplace_back(0);
//			}
//
//			if(vcf_reader.var_ed_idx != vcf_reader.var_st_idx){
//				//clear result list
//				unique_hap_seq_set.clear();
//				unique_alt_kmer_set.clear();
//				ref_hap_seq_string.clear();
//				//S1: loading reference
//				hap_handler.set_ref(vcf_reader.get_chr_ID(), load_region_st, load_region_ed, o.step_length, o.kmer_size);
//				//S2.0 store ref string
//				ref_hap_seq_string.append(hap_handler.ref_pointer);
//				unique_hap_seq_set.insert(ref_hap_seq_string);
//
//				//S2: generate alt string for each haplotype
//				int GT_size = vcf_reader.get_sample_number() * 2;
//				for(int h = 0; h < GT_size; h++){
//					//generate haplotype seq for this hap
//					hap_handler.get_haplotype_seq(vcf_reader.region_GT, h, vcf_reader.var_st_idx, vcf_reader.var_ed_idx);
//					//store hap_string
//					unique_hap_seq_set.insert(hap_handler.HAP_str);
//					//generate alt string
//					hap_handler.store_alt_string(unique_alt_kmer_set);
//				}
//				//S2.1 wrong check
//				if(false && hap_handler.with_wrong_check)
//					for(GT_analysis & var : vcf_reader.region_GT){	var.print(); }
//				//output hap_string into file
//				if(show_strings){
//					for(auto & hap_seq: unique_hap_seq_set){
//						total_hap_seq_number++;
//						if(total_hap_seq_number%10 == 1){fprintf(stdout, ">1\n");}
//						fprintf(stdout, "%sN", hap_seq.c_str());
//						if(total_hap_seq_number%10 == 0){fprintf(stdout, "\n");}
//					}
//				}
//				//output unique kmers
//				if(show_strings){
//					for(auto & kmer: unique_alt_kmer_set){
//						total_alt_string_number++;
//						if(total_alt_string_number%10 == 1){fprintf(stdout, ">1\n");}
//						if(o.show_st_ed_info){
//							fprintf(stdout, "%s\t", kmer.c_str());
//						}else{
//							split_string(alt_string_item_value, temp, kmer.c_str(), "_");
//							fprintf(stdout, "%sN", alt_string_item_value[3].c_str());
//							alt_string_length_count[alt_string_item_value[3].size()]++;
//						}
//						if(total_alt_string_number%10 == 0){fprintf(stdout, "\n");}
//					}
//				}
//				// run log:
//				fprintf(stderr, "load_region_st %d load_region_ed %d %s, VAR count %d ,hap_count %ld \n", load_region_st, load_region_ed, hap_handler.ref_pointer, vcf_reader.var_ed_idx - vcf_reader.var_st_idx, unique_hap_seq_set.size() - 1);
//				if(show_stats){
//					//count total number of hap-seq
//					hap_number_count[unique_hap_seq_set.size() - 1] ++;
//					//count size distribution for hap string
//					for(auto & kmer: unique_alt_kmer_set){
//						split_string(alt_string_item_value, temp, kmer.c_str(), "_");
//						alt_string_length_count[alt_string_item_value[3].size()]++;
//					}
//				}
//			}
//			if(show_stats && vcf_reader.var_ed_idx == vcf_reader.var_st_idx) hap_number_count[0]++;
//			if(vcf_reader.reach_EOF) break;
//		}
//
//		fprintf(stdout, "\n");
//		if(show_stats){
//			//size counting
//			for(uint i = 0; i < hap_number_count.size();i++) fprintf(stderr, "count %d @ %d \n", hap_number_count[i], i);
//			for(uint i = 0; i < hap_number_count.size();i++) fprintf(stderr, "%d ", hap_number_count[i]);
//			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOGs END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			for(uint i = 0; i < alt_string_length_count.size();i++) fprintf(stderr, "count %d @ %d \n", alt_string_length_count[i], i);
//			for(uint i = 0; i < alt_string_length_count.size();i++) fprintf(stderr, "%d ", alt_string_length_count[i]);
//			fprintf(stderr, "\n");
//		}
//
//		vcf_reader.destroy();
//		hap_handler.destroy();
//		//output:
//		return 0;
//	}
//};

#define SV_IDX_bucket_size 0x1fff //8K, 3G / 8k = 0.375M(blocks)
struct MM_idx_loader{
    //PART1： deBGA index
    //(2)UNITIG IDs
    //deBGA_INDEX unitig_idx;
    //PART2:UNITIG MM index
    mm_96_t* ref_mm;		uint64_t ref_mm_size;
    uint32_t* ref_mm_idx;	uint64_t ref_mm_idx_size;
    //PART3: hap strings
    window_block_info* wb_info; uint64_t wb_info_size;//the total number of window block
    uint8_t * wb_data; uint64_t wb_data_size;
    //PART3.1: hap search
    std::vector<int> chr_bg_wb_ID;
    int total_chr_number = 0;//total number of chr, normally, this number is 24 (22 + X + Y)
    uint32_t SV_wb_bg = 0; //the begin window block id of SV
    uint32_t SV_wb_ed = 0;//the end window block id of SV
    int step_length = 0;
    //uint32_t* sv_wb_idx;	uint64_t sv_wb_idx_size;
    bool idx_with_SV = false;
    std::vector<std::vector<uint32_t>> sv_search_idx;

    //part4: ALT SEQ MM IDX
    mm_96_t* alt_mm; uint64_t alt_mm_size;
    uint32_t* alt_mm_idx; uint64_t alt_mm_idx_size;
    //part5: reference
    Simple_ref_handler ref;

private:
    //load from file, return the true load data size(in byte)
    uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, int skip_byte, bool must_load){
        fprintf(stderr, "END loading %s\n", fn);
        char full_fn[1024] = {0};
        strcpy(full_fn, path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        strcat(full_fn, fn);
        FILE *fp_us_b = NULL;
        if(must_load){
            fp_us_b = xopen (full_fn, "rb" );
        }else{
            fp_us_b = fopen (full_fn, "rb" );
            if(fp_us_b == NULL){
                (*data) = NULL;
                return 0;
            }
        }
        fseek(fp_us_b, 0, SEEK_END);// non-portable
        uint64_t file_size = ftell(fp_us_b);
        rewind(fp_us_b);
        fseek(fp_us_b, skip_byte, SEEK_SET);// non-portable
        (*data) = (uint64_t* ) xmalloc (file_size - skip_byte);
        xread ((*data), 1, file_size - skip_byte, fp_us_b);
        fclose(fp_us_b);
        return file_size;
    }

    uint64_t load_data_from_file_part(const char * path_name, const char * fn, void ** data, uint64_t bg, uint64_t ed){
        char full_fn[1024] = {0};
        strcpy(full_fn, path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        strcat(full_fn, fn);
        FILE *fp_us_b = xopen (full_fn, "rb" );
        fseek(fp_us_b, bg, SEEK_SET);// non-portable
        (*data) = (uint64_t* ) xmalloc (ed - bg);
        xread ((*data), 1, ed - bg, fp_us_b);
        fclose(fp_us_b);
        return ed - bg;
    }

    void* load_chr_bg_wb_ID(const char * idx_path_name){
        chr_bg_wb_ID.clear();
        char full_fn[1024] = {0};
        strcpy(full_fn, idx_path_name);
        if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
        strcat(full_fn, "chr_bg_wb_ID.bin");
        fprintf(stderr, "full fn : %s\n", full_fn);
        load_int_from_file(full_fn, chr_bg_wb_ID, MAX_uint64_t);
        return (void *)this;
    }
public:
    //bool skip_unitig_mm_idx, when set to be true, only load wb index, when set false, load all index
    //set unitig_path_name to NULL when only load wb data
    void * load_all_index(const char * idx_path_name, const char * unitig_path_name, bool skip_unitig_mm_idx){
        //PART3: hap strings
        ref.load_gap_file(idx_path_name);
        wb_info_size = load_data_from_file(idx_path_name, "wb_info.bin", (void **)&wb_info, 8, true) / sizeof(window_block_info);
        wb_data_size = load_data_from_file(idx_path_name, "wb_data.bin", (void **)&wb_data, 0, true);
        //PART3.1: wb data
        step_length = wb_info[0].region_length/2;

        load_chr_bg_wb_ID(idx_path_name);

        if(wb_info[wb_info_size - 1].is_SV()){
            idx_with_SV = true;
            total_chr_number = chr_bg_wb_ID.size() - 2;
        }
        else
            total_chr_number = chr_bg_wb_ID.size() - 1;

        ref.load_bin_ref(idx_path_name);


        if(idx_with_SV){//with SV, and build index for it
            SV_wb_bg = chr_bg_wb_ID[total_chr_number];
            SV_wb_ed = chr_bg_wb_ID[total_chr_number + 1];
            sv_search_idx.resize(total_chr_number);
            for(uint chr_ID = 0; chr_ID < sv_search_idx.size(); chr_ID++){
                uint32_t block_number = ref.get_chr_length(chr_ID)/(SV_IDX_bucket_size) + 2;
                sv_search_idx[chr_ID].resize(block_number);
            }
            for(uint wb_id = SV_wb_bg; wb_id < SV_wb_ed; wb_id++){
                if(wb_info[wb_id].chrID >= total_chr_number) break;
                uint32_t blockID = wb_info[wb_id].region_st/(SV_IDX_bucket_size);
                sv_search_idx[wb_info[wb_id].chrID][blockID + 1] = wb_id + 1;
            }
            uint32_t old_idx = SV_wb_bg;
            for(std::vector<uint32_t> & chr :sv_search_idx){
                for(uint32_t & idx : chr){
                    if(idx == 0){	idx = old_idx;}
                    else{			old_idx = idx;}
                }
            }
        }else{
            SV_wb_bg = SV_wb_ed = wb_info_size;
        }

        if(!skip_unitig_mm_idx){
            //PART1： deBGA index
            //unitig_idx.load_index_file(unitig_path_name);
            //part5: reference
            //PART2:UNITIG MM index
            ref_mm_size = load_data_from_file(idx_path_name, "ref_mm.bin", (void **)&ref_mm, 0, true) / sizeof(mm_96_t);
            ref_mm_idx_size = load_data_from_file(idx_path_name, "ref_mm_idx.bin", (void **)&ref_mm_idx, 0, true) / sizeof(uint32_t);
            //part4: ALT SEQ MM IDX
            alt_mm_size = load_data_from_file(idx_path_name, "alt_mm.bin", (void **)&alt_mm, 0, true) / sizeof(mm_96_t);
            alt_mm_idx_size = load_data_from_file(idx_path_name, "alt_mm_idx.bin", (void **)&alt_mm_idx, 0, true) / sizeof(uint32_t);
            //part 3.1: sv idx
            //sv_wb_idx_size = load_data_from_file(idx_path_name, "sv_wb_idx.bin", (void **)&sv_wb_idx, 0, false);
        }
        return (void *)this;
    }

    void load_window_ID_idx(const char * idx_path_name){
        //try load wb_info
        wb_info_size = load_data_from_file_part(idx_path_name, "wb_info.bin", (void **)&wb_info, sizeof(window_block_info)*0 + 8, sizeof(window_block_info)*1 + 8) / sizeof(window_block_info);
        step_length = wb_info[0].region_length/2;
        if(wb_info != NULL) {free (wb_info); wb_info = NULL;}
        load_chr_bg_wb_ID(idx_path_name);
    }

    void load_part_index(const char * idx_path_name, uint32_t wb_bg, uint32_t wb_ed){
        wb_info_size = load_data_from_file_part(idx_path_name, "wb_info.bin", (void **)&wb_info, sizeof(window_block_info)*wb_bg + 8, sizeof(window_block_info)*(wb_ed+1) + 8) / sizeof(window_block_info);
        wb_data_size = load_data_from_file_part(idx_path_name, "wb_data.bin", (void **)&wb_data, wb_info[0].global_offset, wb_info[wb_ed - wb_bg].global_offset);
        wb_info_size -= 1;
        uint64_t wb_data_bg_offset = wb_info[0].global_offset;
        for(int i = 0; i < wb_info_size; i++){
            wb_info[i].global_offset -= wb_data_bg_offset;
        }
    }
    //return WB ID:
    int convert_ref_coordinate_to_wb_coordinate(int chr_ID, uint &local_offset, uint REF_BASE){
        //the begin of one CHR
        if(local_offset == 0){
            return chr_bg_wb_ID[chr_ID];
        }
        //IS (local_offset-1); NOT (local_offset):
        //when wb_offset is 0, always reset to 149 and wb_ID-=1;
        uint wb_offset = local_offset;
        int wb_ID = chr_bg_wb_ID[chr_ID] + (local_offset-1) / step_length;
        wb_offset -= (wb_info[wb_ID].region_st - 1);
        if(wb_offset <= 0 || wb_offset >= 1023){ fprintf(stderr, "FATAL ERROR 1: wb_ID %d, POS %d\n", wb_ID, wb_offset);
            if(wb_ID != 0)
                xassert(0,"");
        }
        wb_offset += REF_BASE;
        local_offset = wb_offset;
        return wb_ID;
    }

    //return chr ID and POS
    void convert_wb_coordinate_to_ref_coordinate(int wb_ID, int &chr_ID, uint &wb_bg_offset){
        for(chr_ID = 0; chr_ID < total_chr_number; chr_ID++){
            if(chr_bg_wb_ID[chr_ID + 1] >= wb_ID){
                wb_bg_offset = (wb_ID - chr_bg_wb_ID[chr_ID])*step_length + 1;
                break;
            }
        }
    }

};

struct String_list_and_var_list{
    std::vector<std::string> hap_string_l;
    std::vector<std::vector<HAP_VAR_ITEM>> var_l;
    void print(FILE * out){
        for(int i = 0; i < hap_string_l.size(); i++){
            for(auto & v: var_l[i]){
                v.print(out);
            }
            fprintf(out, "\n: %s\n", hap_string_l[i].c_str());
        }
    }

};

struct hap_string_loader_single_thread{
private:
    //PART3: hap strings
    Window_block_handler wbh;
    //load from file, return the true load data size(in byte)
    std::vector<std::string> hap_string_buff;
    String_list_and_var_list slvl;
    std::vector<std::vector<HAP_VAR_ITEM>> var_l;
    std::vector<std::string> & get_string_list_core(uint32_t window_ID, std::string & window_ref, window_block_info* wb_info, uint8_t* wb_data, bool generate_var_list){
        hap_string_buff.clear();
        window_block_info & c_wb = wb_info[window_ID];
        wbh.unfold_window_block_all(wb_data + c_wb.global_offset);
        //load reference:
        std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
        //clear minimizer builder:
        if(generate_var_list){
            var_l.clear();
        }
        if(c_wb.is_SV()){
            std::vector<char> & hap_string = wbh.unfold_hap_seq_SV(hap_l[0].p, hap_l[0].len, window_ref, c_wb, false, 0);
            hap_string_buff.emplace_back(&(hap_string[0]));
        }else{
            for(auto & hap:hap_l){//load all haplotypes
                std::vector<char> & hap_string = wbh.unfold_hap_seq(hap.p, hap.len, window_ref, false, generate_var_list, 0);
                if(generate_var_list){
                    var_l.emplace_back(wbh.var_l);
                }
                hap_string.emplace_back(0);
                hap_string_buff.emplace_back(&(hap_string[0]));
            }
        }
        if(generate_var_list){
            std::swap(slvl.hap_string_l, hap_string_buff);
            std::swap(slvl.var_l, var_l);
        }
        return hap_string_buff;
    }

    std::vector<char> & get_string_core(uint32_t window_ID, uint32_t hap_ID, std::string & window_ref, window_block_info* wb_info, uint8_t* wb_data){
        window_block_info & c_wb = wb_info[window_ID];
        wbh.unfold_window_block_all(wb_data + c_wb.global_offset);
        //load reference:
        std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
        if(c_wb.is_SV()){
            xassert(hap_ID == 0, "SV has only one haplotype in one window block");
            return wbh.unfold_hap_seq_SV(hap_l[0].p, hap_l[0].len, window_ref, c_wb, false, 0);
        }else{
            return wbh.unfold_hap_seq(hap_l[hap_ID].p, hap_l[hap_ID].len, window_ref, false, false, 0);
        }
    }

    std::vector<std::vector<uint8_t>> hap_bin_buff;
    std::vector<std::vector<uint8_t>> & get_bin_data_core(uint32_t window_ID,  window_block_info* wb_info, uint8_t* wb_data){
        hap_bin_buff.clear();
        window_block_info & c_wb = wb_info[window_ID];
        wbh.unfold_window_block_all(wb_data + c_wb.global_offset);
        //load reference:
        std::vector<Hap_seq_address> & hap_l = wbh.block_all_hap_list;
        //clear minimizer builder:
        for(auto & hap:hap_l){//load all haplotypes
            hap_bin_buff.emplace_back();
            std::vector<uint8_t> c_bin_hap = hap_bin_buff.back();
            c_bin_hap.insert(c_bin_hap.begin(), hap.p, hap.p + hap.len);
            hap_bin_buff.emplace_back(c_bin_hap);
        }
        return hap_bin_buff;
    }

    //step_length is 150 normally
    //return the corresponding wb ID
    //
    int get_windows_ID_core(uint8_t chr_ID, int POS, std::vector<int> &chr_bg_wb_ID, int step_length){
        return chr_bg_wb_ID[chr_ID] + (POS - 1)/step_length;
    }

    std::vector<uint32_t> SV_wb_list_buff;
    //search SVs in a region:
    //search all SVs in region chr_ID:POS-END;
    //the result (WB ID) stored in wb_bg(include) and wb_ed (not include)
    void search_SV_IN_region_core(uint8_t chr_ID, uint32_t POS, uint32_t END,  window_block_info* wb_info, std::vector<std::vector<uint32_t>> &sv_search_idx, uint32_t &wb_bg, uint32_t &wb_ed){
        uint32_t search_bg = sv_search_idx[chr_ID][POS/SV_IDX_bucket_size];
        uint32_t search_ed = sv_search_idx[chr_ID][END/SV_IDX_bucket_size + 1];
        for(uint32_t i = search_bg; i < search_ed;i++){
            if(wb_info[i].region_st >= POS){
                wb_bg = i; break;
            }
        }
        for(uint32_t i = search_ed - 1; i >= search_bg; i--){
            if(wb_info[i].region_st < END){
                wb_ed = i; break;
            }
        }
        wb_ed += 1;
    }
public:
    //function used in python : 1
    int get_windows_ID(uint8_t chr_ID, int POS, void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        return get_windows_ID_core(chr_ID, POS, i->chr_bg_wb_ID, i->step_length);
    }
    //search SVs in a region:
    //search all SVs in region chr_ID:POS-END;
    //the result (WB ID) stored in wb_bg(include) and wb_ed (not include)
    //return the SV number
    std::vector<uint32_t> &search_SV_IN_region(uint8_t chr_ID, uint32_t POS, uint32_t END, void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        uint32_t wb_bg = UINT32_MAX; uint32_t wb_ed = 0;
        search_SV_IN_region_core(chr_ID, POS, END, i->wb_info, i->sv_search_idx, wb_bg, wb_ed);
        SV_wb_list_buff.clear();
        for(uint32_t i = wb_bg; i < wb_ed; i++){
            SV_wb_list_buff.emplace_back(i);
        }
        return SV_wb_list_buff;
    }
    std::vector<char> & get_string(uint32_t window_ID, uint32_t hap_ID, std::string & window_ref, void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        return get_string_core(window_ID, hap_ID, window_ref, i->wb_info, i->wb_data);
    }
    //function used in python : 1
    String_list_and_var_list & get_string_list_and_var_list(uint32_t window_ID, std::string & window_ref, void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        get_string_list_core(window_ID, window_ref, i->wb_info, i->wb_data, true);
        return slvl;
    }

    //function used in python : 1
    std::vector<std::string> & get_string_list(uint32_t window_ID, std::string & window_ref, void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        return get_string_list_core(window_ID, window_ref, i->wb_info, i->wb_data, false);
    }

    //function used in python : 2
    std::vector<std::vector<uint8_t>> & get_bin_data_list(uint32_t window_ID,  void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        return get_bin_data_core(window_ID,  i->wb_info, i->wb_data);
    }

    uint64_t get_hap_N(uint32_t window_ID,  void* i_){
        MM_idx_loader * i = (MM_idx_loader *)i_;
        hap_string_buff.clear();
        window_block_info & c_wb = i->wb_info[window_ID];
        wbh.unfold_window_block_all(i->wb_data + c_wb.global_offset);
        return wbh.block_all_hap_list.size();
    }

public:
    //debug code, not use
    void printf_modify_info_after_get_string(){
        wbh.printf_modify_info(stderr);
    }

};

struct Mm_seed_core{
    void set(bool is_same_direction_, int read_pos_, uint32_t seed_idx_){
        is_same_direction = is_same_direction_;
        read_pos = read_pos_;
        seed_idx = seed_idx_;
    }
    bool is_same_direction;
    int read_pos;
    uint32_t seed_idx;
    void print(FILE * out){
        fprintf(out, "read_pos: %d, is_same_direction %d ", read_pos, is_same_direction);
    }
};

struct Mm_seed_ALT{
    Mm_seed_ALT(bool is_same_direction_, int read_pos_, uint32_t window_ID_, uint32_t hap_ID_, uint32_t alt_offset_, uint32_t seed_idx_){
        mm.set(is_same_direction_, read_pos_, seed_idx_);
        window_ID = window_ID_;
        hap_ID = hap_ID_;
        alt_offset = alt_offset_;
    }
    Mm_seed_core mm;
    uint32_t window_ID;
    uint32_t hap_ID;
    uint32_t alt_offset;
    void print(FILE * out){
        mm.print(out);
        fprintf(out, "window_ID: %d, hap_ID: %d, alt_offset: %d\n", window_ID, hap_ID, alt_offset);
    }
};

struct Mm_seed_REF{
    Mm_seed_REF(bool is_same_direction_, int read_pos_, uint32_t unitig_ID_, uint32_t alt_offset_, uint32_t seed_idx_){
        mm.set(is_same_direction_, read_pos_, seed_idx_);
        unitig_ID = unitig_ID_;
        unitig_offset = alt_offset_;
    }
    Mm_seed_core mm;
    uint32_t unitig_ID;
    uint32_t unitig_offset;
    void print(FILE * out){
        mm.print(out);
        fprintf(out, "unitig_ID: %d, unitig_offset: %d\n", unitig_ID, unitig_offset);
    }
};

struct Seed_handler{
public:
    std::vector<uint8_t> seed_with_rst_bool;
    std::vector<uint32_t> seed_hit_number;//the hit number of each minimizer, used only in the alt seq search, each windows used as one results
    std::vector<uint8_t> seed_hit_right;//whether the seed is right, used for evaluation

    //std::set<uint32_t> uniq_winddow_set;//the hit number of each minimizer, used only in the alt seq search, each windows used as one results

    std::vector<Mm_seed_ALT> same_dir_seed_alt;
    std::vector<Mm_seed_ALT> not_same_dir_seed_alt;

    std::vector<Mm_seed_REF> same_dir_seed_ref;
    std::vector<Mm_seed_REF> not_same_dir_seed_ref;

    uint32_t max_ref_hit = 0;
    uint32_t max_alt_hit = 0;

private:
    //read minimizer buff;
    Minimizer_generater mg;
    std::vector<mm128_t> read_mm_v;

    int index_length = 0;//26 bit
    int minimizer_length = 0; //at most 56 bit
    int mini_left_length = 0; // at most 30 bp
    int mini_left_mask = 0;//0x3fffffff when left == 30; and 0x3ffff when left = 18;
public:
    //debug code: eval
    int get_mm_kmer_len(){  return minimizer_length/2;  }
    void show_read_mm(){
        for(auto & m: read_mm_v){
            m.print(stderr);
        }
    }

public:
    void init(int window_size_, int mini_kmer_len_, int is_hpc_, uint32_t max_ref_hit_,	uint32_t max_alt_hit_){
        mg.init(window_size_, mini_kmer_len_, is_hpc_);
        index_length = 26;//26 bit
        minimizer_length = mini_kmer_len_ * 2; //at most 56 bit
        mini_left_length = minimizer_length - index_length; // at most 30 bp
        mini_left_mask = (0x1 << mini_left_length) - 1;//0x3fffffff when left == 30; and 0x3ffff when left = 18;

        max_ref_hit = max_ref_hit_;
        max_alt_hit = max_alt_hit_;
    }

    void read_mm(const char * query, int q_len){
        read_mm_v.clear();
        mg.mm_sketch(query, q_len, read_mm_v);
        seed_with_rst_bool.resize(read_mm_v.size());
        memset(&(seed_with_rst_bool[0]), 0, seed_with_rst_bool.size());
        seed_hit_number.resize(read_mm_v.size());

        seed_hit_right.resize(read_mm_v.size());
        memset(&(seed_hit_right[0]), 0, seed_hit_right.size());
    }
private:
    int binsearch_range(uint32_t key, mm_96_t* v, int64_t n,  int64_t *range, int8_t k_off)
    {
        int64_t l=0, r=n-1, m;
        uint32_t tmp = 0;
        range[0] = range[1] = -1;

        while (l <= r)
        {
            m = (l+r)/2;
            tmp = (v[m].x >> k_off);
            if (tmp == key)
            {
                range[0] = range[1] = m;
                //run low bound
                int64_t sl=l, sr=m-1, sm;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = (v[sm].x >> k_off);
                    if (tmp == key){ range[0] = sm; sr = sm-1;}
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }
                //run upper bound
                sl = m+1; sr = r;
                while (sl <= sr){
                    sm = (sl+sr)/2;
                    tmp = (v[sm].x >> k_off);
                    if (tmp == key){ range[1] = sm; sl = sm+1; }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }
                return 1;
            }
            else if (tmp > key) r = m - 1;
            else l = m + 1;
        }
        return -1;
    }

private:
    void bin_search_idx(uint32_t*mm_idx, mm_96_t* mm, mm128_t *read_mm, mm_96_t* mm_bg, mm_96_t* mm_ed, uint64_t &read_is_forward, uint32_t &read_offset){
        mm_bg = NULL; mm_ed = NULL;

        uint64_t read_hash_key = (read_mm->x >> 8);
        uint64_t read_idx_bucket = (read_hash_key >> mini_left_length);
        uint64_t read_left_key = (read_hash_key & mini_left_mask);

        //search mm idx
        //binary search the mm seed
        uint32_t idx_bucket_bg = mm_idx[read_idx_bucket];
        uint32_t idx_bucket_ed = mm_idx[read_idx_bucket + 1];
        //search rst:
        int64_t search_range_rst[2];

        if(-1 == binsearch_range(read_left_key, mm + idx_bucket_bg, idx_bucket_ed - idx_bucket_bg, search_range_rst, 2))
            return;

        mm_bg = mm + idx_bucket_bg + search_range_rst[0];
        mm_ed = mm + idx_bucket_bg + search_range_rst[1];

        read_is_forward = (read_mm->y & 0x1);
        read_offset = ((read_mm->y & 0xffffffff) >> 1);
    }

public:
    uint64_t read_mm_size(){ return read_mm_v.size(); }
    uint32_t get_read_mm_offset(int read_mm_idx){
        return ((read_mm_v[read_mm_idx].y & 0xffffffff) >> 1);
    }

    void search_seed_single_minimizer_ref(int read_mm_idx, MM_idx_loader * idx, bool is_append,
            std::vector<Mm_seed_REF> &same_dir_seed_ref,
            std::vector<Mm_seed_REF> &not_same_dir_seed_ref){

        if(!is_append){
            same_dir_seed_ref.clear();
            not_same_dir_seed_ref.clear();
        }

        mm_96_t* mm_bg, * mm_ed;
        uint64_t read_is_forward = 0; uint32_t read_offset = 0;
        bin_search_idx(idx->ref_mm_idx, idx->ref_mm, &(read_mm_v[read_mm_idx]) , mm_bg, mm_ed, read_is_forward, read_offset);

        if(mm_bg == NULL){	return; }

        for(mm_96_t* c_mm = mm_bg;c_mm <= mm_ed; c_mm++){
            uint64_t mm_is_forward = (c_mm->x & 0x1);
            uint64_t UNITIG_ID = (c_mm->y);
            uint64_t UNITIG_offset = (c_mm->z);
            if((read_is_forward == mm_is_forward)){
                same_dir_seed_ref.emplace_back(true, read_offset, UNITIG_ID, UNITIG_offset, read_mm_idx);
            }else{
                not_same_dir_seed_ref.emplace_back(false, read_offset, UNITIG_ID, UNITIG_offset, read_mm_idx);
            }
        }
    }

public:

    void search_seed_single_minimizer_alt(int read_mm_idx, MM_idx_loader * idx, bool is_append, int max_hit,
            std::vector<Mm_seed_ALT> &same_dir_seed_alt, std::vector<Mm_seed_ALT> &not_same_dir_seed_alt){
        if(!is_append){
            same_dir_seed_alt.clear();
            not_same_dir_seed_alt.clear();
        }
        uint64_t read_is_forward = 0; uint32_t read_offset = 0;
        mm_96_t* mm_bg, * mm_ed;
        bin_search_idx(idx->alt_mm_idx, idx->alt_mm, &(read_mm_v[read_mm_idx]) , mm_bg, mm_ed, read_is_forward, read_offset);
        if(mm_bg == NULL){	return; }
        if(mm_ed + 1 > mm_bg + max_hit){ return;}

        for(mm_96_t* c_mm = mm_bg;c_mm <= mm_ed; c_mm++){
            uint64_t mm_is_forward = (c_mm->x & 0x1);
            uint64_t POS_bit_64 = (((uint64_t)c_mm->y) << 32) + ((uint64_t)c_mm->z);
            uint32_t window_ID = (POS_bit_64 >> 39) - 1;
            uint32_t hap_ID = ((POS_bit_64 >> 20) & 0x7ffff) - 1;
            uint32_t alt_offset = ((POS_bit_64) & 0xfffff);
            if((read_is_forward == mm_is_forward)){
                same_dir_seed_alt.emplace_back(true, read_offset, window_ID, hap_ID, alt_offset, read_mm_idx);
            }else{
                not_same_dir_seed_alt.emplace_back(false, read_offset, window_ID, hap_ID, alt_offset, read_mm_idx);
            }
        }
    }

    //same as search seed
    void search_seed_all(MM_idx_loader * idx, bool is_ref_not_alt, bool skip_already_found_seeds){
            //get minimizer list from seed
            int seed_step_length = 0;
            uint32_t*mm_idx = NULL;
            mm_96_t* mm = NULL;

            if(is_ref_not_alt){
                mm_idx = idx->ref_mm_idx;
                mm = idx->ref_mm;
                seed_step_length = 1;
                same_dir_seed_ref.clear();
                not_same_dir_seed_ref.clear();
            }else{
                mm_idx = idx->alt_mm_idx;
                mm = idx->alt_mm;
                seed_step_length = 1;
                same_dir_seed_alt.clear();
                not_same_dir_seed_alt.clear();
            }

            for(uint32_t seed_idx = 0; seed_idx < read_mm_v.size(); seed_idx += seed_step_length){
                if(skip_already_found_seeds && seed_with_rst_bool[seed_idx] > 0)
                    continue;
                if(is_ref_not_alt){
                    search_seed_single_minimizer_ref(seed_idx, idx, true, same_dir_seed_ref, not_same_dir_seed_ref);
                }else{
                    uint64_t same_dir_seed_alt_size = same_dir_seed_alt.size();
                    uint64_t not_same_dir_seed_alt_size = not_same_dir_seed_alt.size();
                    search_seed_single_minimizer_alt(seed_idx, idx, true, 500, same_dir_seed_alt, not_same_dir_seed_alt);
                    seed_hit_number[seed_idx] = same_dir_seed_alt.size() - same_dir_seed_alt_size + not_same_dir_seed_alt.size() - not_same_dir_seed_alt_size;
                }
            }
        }

    void search_seed(MM_idx_loader * idx, bool is_ref_not_alt, bool skip_already_found_seeds){
        //get minimizer list from seed
        int seed_step_length = 0;
        uint32_t*mm_idx = NULL;
        mm_96_t* mm = NULL;

        if(is_ref_not_alt){
            mm_idx = idx->ref_mm_idx;
            mm = idx->ref_mm;
            seed_step_length = 1;
            same_dir_seed_ref.clear();
            not_same_dir_seed_ref.clear();
        }else{
            mm_idx = idx->alt_mm_idx;
            mm = idx->alt_mm;
            seed_step_length = 1;
            same_dir_seed_alt.clear();
            not_same_dir_seed_alt.clear();
        }

        for(uint32_t seed_idx = 0; seed_idx < read_mm_v.size(); seed_idx += seed_step_length){
            if(skip_already_found_seeds && seed_with_rst_bool[seed_idx] > 0)
                continue;

            mm128_t *c_mm = &(read_mm_v[seed_idx]);
            uint64_t read_hash_key = (c_mm->x >> 8);
            uint64_t read_idx_bucket = (read_hash_key >> mini_left_length);
            uint64_t read_left_key = (read_hash_key & mini_left_mask);
            uint64_t read_is_forward = (c_mm->y & 0x1);
            uint32_t read_offset = ((c_mm->y & 0xffffffff) >> 1);

            //search mm idx
            //binary search the mm seed
            uint32_t idx_bucket_bg = mm_idx[read_idx_bucket];
            uint32_t idx_bucket_ed = mm_idx[read_idx_bucket + 1];
            //search rst:
            int64_t search_range_rst[2];

            if(-1 == binsearch_range(read_left_key, mm + idx_bucket_bg, idx_bucket_ed - idx_bucket_bg, search_range_rst, 2))
                continue;

            mm_96_t* mm_bg = mm + idx_bucket_bg + search_range_rst[0];
            mm_96_t* mm_ed = mm + idx_bucket_bg + search_range_rst[1];

            //set the seed already has a result
            if(mm_bg <= mm_ed)
                seed_with_rst_bool[seed_idx] |= ((is_ref_not_alt == true)?1:2);
            //store in the list
            if(is_ref_not_alt){
                for(mm_96_t* c_mm = mm_bg;c_mm <= mm_ed; c_mm++){
                    uint64_t mm_is_forward = (c_mm->x & 0x1);

                    uint64_t UNITIG_ID = (c_mm->y);
                    uint64_t UNITIG_offset = (c_mm->z);

                    if((read_is_forward == mm_is_forward)){
                        same_dir_seed_ref.emplace_back(true, read_offset, UNITIG_ID, UNITIG_offset, seed_idx);
                    }else{
                        not_same_dir_seed_ref.emplace_back(false, read_offset, UNITIG_ID, UNITIG_offset, seed_idx);
                    }
                }
            }else{
                if(mm_ed + 1 > mm_bg + max_alt_hit){ break; }
                for(mm_96_t* c_mm = mm_bg;c_mm <= mm_ed; c_mm++){

                    uint64_t mm_is_forward = (c_mm->x & 0x1);
                    uint64_t POS_bit_64 = (((uint64_t)c_mm->y) << 32) + ((uint64_t)c_mm->z);
                    uint32_t window_ID = (POS_bit_64 >> 39) - 1;
                    uint32_t hap_ID = ((POS_bit_64 >> 20) & 0x7ffff) - 1;
                    uint32_t alt_offset = ((POS_bit_64) & 0xfffff);
                    if((read_is_forward == mm_is_forward)){
                        same_dir_seed_alt.emplace_back(true, read_offset, window_ID, hap_ID, alt_offset, seed_idx);
                    }else{
                        not_same_dir_seed_alt.emplace_back(false, read_offset, window_ID, hap_ID, alt_offset, seed_idx);
                    }
                }
                seed_hit_number[seed_idx] = mm_ed - mm_bg + 1;
            }
        }
    }
};

int hapseq2fa(int argc, char *argv[]);
int print_var_list(int argc, char *argv[]);
int print_var_list_compect(int argc, char *argv[]);
}
#endif /* HAPLOTYPE_HPP_ */
