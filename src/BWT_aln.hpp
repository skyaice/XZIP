/*
 * BWT_aln.hpp
 *
 *  Created on: 2022年5月26日
 *      Author: fenghe
 */

 #ifndef BWT_ALN_HPP_
 #define BWT_ALN_HPP_
 
 #include <stdint.h>
 #include <stdio.h>
 #include <getopt.h>
 #include <stdlib.h>
 #include <cstring>
 #include <ctype.h>
 
 #include <vector>
 #include <map>
 #include "CPPLIB/get_option_cpp.hpp"
 #include "BWT_idx/bwt.hpp"
 #include "process2_sa.hpp"
 extern "C"
 {
 #include "clib/bam_file.h"
 #include "clib/desc.h"
 }
 #include "./htslib/htslib/sam.h"
 #include "process0_idx_loader.hpp"
 #include "occ_RST_DEF.hpp"
 //#include "ReadHandler.hpp"
 
 
 namespace BWT_aln{
 
 #define SEED_OFFSET (K_T - LEN_KMER)
 #define LEN_KMER_LEFT (32 - LEN_KMER) //LEN_KMER_LEFT = 32 - LEN_KMER
 #define SEED_STEP 5
 #define UNI_POS_N_MAX 32
 
 #define GAP_OPEN 16
 #define GAP_EXT 1
 #define GAP_OPEN2 32
 #define GAP_EXT2 0
 #define MATCH_SCORE 2
 #define MISMATCH_SCORE 12
 
 #define ZDROP_SCORE 400
 #define BANDWIDTH 500
 
 #define OUTPUT "./aln.sam"
 
 #define MINIMIZER_LEN_MAP 22
 #define MINIMIZER_WINDOW_SIZE_MAP 5
 
 //SEEDING_LEVEL
 #define SEEDING_LEVEL_NO_SEEDING 0
 #define SEEDING_LEVEL_ONLY_UNITIG 1
 #define SEEDING_LEVEL_ONLY_ALT 2
 #define SEEDING_LEVEL_AUTO 3
 
 #define SEEDING_LEVEL_MIN_OUTPUT_LOG 4
 #define SEEDING_LEVEL_AUTO_AND_OUTPUT_LOG 4
 #define SEEDING_LEVEL_ALL_SEED_AND_OUTPUT_LOG 5
 #define SEEDING_LEVEL_UNITIG_ONLY_AND_OUTPUT_LOG 6
 
 #define MAX_ref_hit 500
 #define MAX_alt_hit 500
 
 
 
 
 struct MAP_PARA{
	 //basic option
	 int thread_n;
 
	 //score option
	 int match_D;
	 int mismatch_D;
	 int gap_open_D;
	 int gap_ex_D;
	 int gap_open2_D;
	 int gap_ex2_D;
	 int bw;
	 int zdrop_D; //for DNA zdrop = 400, 200 for RNA
 
	 int max_use_read;
	 //output option
	 char * reference_fa;
	 char * hap_idex_folder;
	 //char * deBGA_indexDir;
	 char * read_fastq1;
	 char * read_fastq2;
 
	 //debug code
 //	int window_size;
 //	int mini_kmer_len;
	 int seeding_level;
 
	 int max_ref_hit;
	 int max_alt_hit;
 
 
 
	 int get_option(int argc, char *argv[]){
		 options_list l;
		 l.show_command(stderr, argc + 1, argv - 1);
		 l.add_title_string("\n");
		 l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  align  [Options] <reference_fa> <haplotype_idx> [ReadFile_1.fq] [ReadFile_2.fq] \n");
		 l.add_title_string("  Basic:   \n");
		 l.add_title_string("    <reference_fa>  FILEs  the BWA index\n");
		 l.add_title_string("    <haplotype_idx> Folder the haplotype index\n");
 
		 l.add_title_string("    [ReadFiles.fq]  FILES  reads files, FASTQ(A)(or fa.gz/fq.gz) \n");
 
		 //thread number
		 l.add_option("thread", 			't', "Number of threads", true, 4); l.set_arg_pointer_back((void *)&thread_n);
		 //score
		 l.add_option("gap-open1", 		'O', "Gap open penalty 1", true, GAP_OPEN); l.set_arg_pointer_back((void *)&gap_open_D);
			 l.l.back().add_help_msg("a k-long gap costs min{O+k*E,P+k*F}..");
		 l.add_option("gap-open2", 		'P', "Gap open penalty 2.", true, GAP_OPEN2); l.set_arg_pointer_back((void *)&gap_open2_D);
		 l.add_option("gap-extension1", 	'E', "Gap extension penalty 1.", true, GAP_EXT);l.set_arg_pointer_back((void *)&gap_ex_D);
		 l.add_option("gap-extension2", 	'F', "Gap extension penalty 2.", true, GAP_EXT2);l.set_arg_pointer_back((void *)&gap_ex2_D);
		 l.add_option("match-score", 	'M', "Match score for SW-alignment.", true, MATCH_SCORE);l.set_arg_pointer_back((void *)&match_D);
		 l.add_option("mis-score", 		'X', "Mismatch score for SW-alignment.", true, MISMATCH_SCORE);l.set_arg_pointer_back((void *)&mismatch_D);
		 l.add_option("zdrop", 			'Z', "Z-drop score for splice/non-splice alignment.", true, ZDROP_SCORE);l.set_arg_pointer_back((void *)&zdrop_D);
		 l.add_option("band-width", 		'B', "Bandwidth used in chaining and DP-based alignment.", true, BANDWIDTH);l.set_arg_pointer_back((void *)&bw);
		 //output
		 l.add_option("max_use_read", 	'R', "Max number of read to alignment, used for debug or test PG.", true, MAX_int32t); l.set_arg_pointer_back((void *)&max_use_read);
 
		 l.add_option("max_ref_hit", 	'J', "Used in seeding. Max hit number of read a minimizer in the single reference.", true, MAX_ref_hit); l.set_arg_pointer_back((void *)&max_ref_hit);
		 l.l.back().add_help_msg("When hit number is over max hit, no result will return");
		 l.add_option("max_alt_hit", 	'K', "Used in seeding. Max hit number of read a minimizer in the single alternative sequence.", true, MAX_alt_hit); l.set_arg_pointer_back((void *)&max_alt_hit);
		 l.l.back().add_help_msg("When hit number is over max hit, no result will return");
 
		 //debug options
		 //l.add_option("minimizer_len", 	'm', "length of the minimizer", true, MINIMIZER_LEN_MAP ); l.set_arg_pointer_back((void *)&mini_kmer_len);
		 //l.add_option("minimizer_window_size", 	'w', "window_size of minimizer", true, MINIMIZER_WINDOW_SIZE_MAP ); l.set_arg_pointer_back((void *)&window_size);
		 l.add_option("seed_level", 		'S', "Seeding level", true, SEEDING_LEVEL_AUTO ); l.set_arg_pointer_back((void *)&seeding_level);
		 l.l.back().add_help_msg("Detail of seeding level:");
		 l.l.back().add_help_msg("SEEDING_LEVEL: NO_SEEDING 0");
		 l.l.back().add_help_msg("SEEDING_LEVEL: ONLY_UNITIG 1");
		 l.l.back().add_help_msg("SEEDING_LEVEL: ONLY_ALT 2");
		 l.l.back().add_help_msg("SEEDING_LEVEL: AUTO 3");
		 l.l.back().add_help_msg("SEEDING_LEVEL: AUTO and OUTPUT LOG 4");
 
		 if(l.default_option_handler(argc, argv)) {
			 l.show_c_value(stderr);
			 return 1;
		 }
		 l.show_c_value(stderr);
 
		 if (argc - optind < 3)
			 return l.output_usage();
		 xassert((thread_n >= 1) && (thread_n <= 256), "Input error: thread_n cannot be less than 1 or more than 256\n");
 
		 reference_fa = strdup(argv[optind]);
		 hap_idex_folder = strdup(argv[optind + 1]);
		 if (hap_idex_folder[strlen(hap_idex_folder) - 1] != '/') strcat(hap_idex_folder, "/");
		 //deBGA_indexDir = strdup(argv[optind + 1]);
		 read_fastq1 = strdup(argv[optind + 2]);
		 read_fastq2 = strdup(argv[optind + 3]);
 
		 return 0;
	 }
 };
 
 struct BWT_MAP_rst{
	 int32_t wb_ID;
	 int32_t POS;
	 int8_t is_rev;
 
	 void set(std::vector<BWT_MAP_rst> & l, int idx){
		 if(idx < 0){
			 wb_ID = -1;
			 POS = -1;
			 is_rev = -1;
		 }else{
			 wb_ID = l[idx].wb_ID;
			 POS = l[idx].POS;
			 is_rev = l[idx].is_rev;
		 }
	 }
	 BWT_MAP_rst(){
		 wb_ID = -1;
		 POS = -1;
		 is_rev = -1;
	 }
 
	 BWT_MAP_rst(int32_t wb_ID_, int32_t POS_, int8_t ori_){
		 wb_ID = wb_ID_;
		 POS = POS_;
		 is_rev = ori_;
		 if(POS_ < 0 || POS_ > 0x3fffff){ fprintf(stderr, "FATAL ERROR 1: wb_ID %d, POS %d\n", wb_ID, POS); xassert(0,"");}
	 }
 
	 int pair_score(BWT_MAP_rst &r2, int SV_wb_bg){
		 int wb_diff = wb_ID - r2.wb_ID;
		 //SV check
		 int score = 2;
		 if(wb_ID >= SV_wb_bg){
			 score -= 1;
		 }
		 if(r2.wb_ID >= SV_wb_bg){
			 score -= 1;
		 }
		 //POSITION check
		 if(wb_diff < 6 && wb_diff > -6){
			 score += 10;
		 }
		 //direction check
		 if(is_rev != r2.is_rev){
			 if((!is_rev && wb_diff <= 0) || (is_rev && wb_diff >= 0)){
				 score += 5;
			 }
		 }
		 return score;
	 }
 
	 int pair_score(int SV_wb_bg){
		 if(wb_ID >= SV_wb_bg){
			 return 0;
		 }
		 return 1;
	 }
 
	 void print(){
		 fprintf(stderr, "chr_ID %d, POS %d, ori %d\n", wb_ID, POS, is_rev);
	 }
 
 
 };
 
 //data used for each ALIGN-PROCESS thread // 24 threads
 struct Match_block{
	 int start_read_pos;
	 uint start_ref_pos;
	 uint8_t match_length;
 
	 Match_block(){
		 start_read_pos = 0;
		 start_ref_pos = 0;
		 match_length = 0;
	 }
	 Match_block(int &start_read_pos_, uint &start_ref_pos_, uint8_t &match_length_){
		 start_read_pos = start_read_pos_;
		 start_ref_pos = start_ref_pos_;
		 match_length = match_length_;
	 }
 };
 
 struct Classify_buff_pool
 {
	 //seed handler
	 std::vector<BWT_MAP_rst> map_v_r1;
	 std::vector<BWT_MAP_rst> map_v_r2;
 
	 std::set<int32_t> uniq_wb;
 
 //	hap_string_loader_single_thread hl_r1;
 //	hap_string_loader_single_thread hl_r2;
	 //seed index pointer
	 Aln_online::IDX_loader *idx;
	 OL_PAR * p;
	 bam_hdr_t *header ;
 
	 //other buffs :: todo
	 //single end handler
	 //ALIGN_BUFF_per_thread t;
	 Aln_online::SA_result  * sa_r;
 
	 // exact gap statics
	 const char* statics_file_names;
	 FILE* statics_file_handler;
 
	 //std::vector<uint16_t> cigar_out[2];
	 std::vector<uint64_t> SA_pos_result_BUFFs;
	 std::vector<BASIC_MAP_rst> read_map_rst[2];
	 std::vector<uint64_t> pairing_rst_list;		
	 std::vector<Aln_online::Gap_Aln_SEED_t> seed_list_BUFF;
	 std::vector<Aln_online::Gap_aln_hit_t> hit_list_BUFF;
	 std::vector<uint64_t> SA_pos_result;
	 std::vector<Gap_aln_Cluster_t> cluster_list_out_BUFF;
	 std::vector<uint8_t> ref_buff_for_hamming;
	 std::vector<std::tuple<uint8_t, uint8_t, uint8_t, uint8_t>> cigar_rst;
	 std::vector<std::tuple<uint8_t, uint, std::string, uint8_t>> md_rst;
	 std::string md_out[2];
	 std::string cigar_out_str[2];
	 char *seq;
	 //BAM_handler *bh;
	 char hs[150];
	 char ts[150];
	 std::string convert_md_str;
	 std::string ref_str;
	 Aln_online::SIMPLE_EDIT_DISTANCE_handler edit_h[2];
	 uint64_t qual_count[4][151];
	 uint64_t read_qual_count[151];
	 uint64_t acgt_count[4];
	 uint64_t acgt_count2[5];
	 std::vector<Match_block> match_block_list;
	 path_t cigar;
	 uint64_t stats;
 
	 void init(){
		 stats = 0;
		 kv_init(cigar);
		 seq = (char*)xcalloc(1, 200);
		 for(int i = 0; i< 4; i++){
			 memset(qual_count[i], 0, 151*sizeof(uint64_t));
			 
		 }
		 memset(acgt_count, 0, 4*sizeof(uint64_t));
		 memset(acgt_count2, 0, 5*sizeof(uint64_t));
		 memset(read_qual_count, 0, 151*sizeof(uint64_t));
		 match_block_list.resize(100);
		 SA_pos_result_BUFFs.resize(100);
		 convert_md_str.resize(100);
		 cigar_rst.resize(128);
		 cluster_list_out_BUFF.resize(128);
		 ref_buff_for_hamming.resize(256);
		 edit_h[0].init(150);
		 edit_h[1].init(150);
	 }
 
	 void clear(){
		 match_block_list.clear();
		 convert_md_str.clear();
		 //free(bh);
		 free(seq);
		 SA_pos_result_BUFFs.clear();
		 // cigar_out[0].clear();
		 // cigar_out[1].clear();
		 md_rst.clear();
		 cigar_rst.clear();
		 read_map_rst[0].clear();
		 read_map_rst[1].clear();
		 pairing_rst_list.clear();
		 //seed_list_BUFF.resize(128);
		 hit_list_BUFF.clear();
		 SA_pos_result.clear();
		 cluster_list_out_BUFF.clear();
		 ref_buff_for_hamming.clear();
		 edit_h[0].destroy();
		 edit_h[1].destroy();
	 }
 };
 
 //debug code:
 
 #define ALT_HIT_POS_MAX 150
 
 struct DEBUG_minimizer_evaluate_per_read{
 //	uint64_t minimizer_N;//total minimizer of read//
 //	uint64_t UNITIG_seed_N;//total seed of read on the UNITIG
 //	//float UNITIG_seed_rate;//the probability of a minimizer find a seed on the UNITIG
 //	uint64_t UNITIG_seed_N_ALT;//total number minimizer of read that search in the alt string
 //	//float ALT_seed_rate;//the probability of a minimizer find a seed on the ALT-SEQ (when fail finding seeds on the UNITIG)
 //	uint64_t ALT_with_seed_minimizer_N;//total number of minimizes that found seeds on the ALT-SEQ (when fail finding seeds on the UNITIG) in a read
 //	uint64_t ALT_seed_N;//total number of seeds on the ALT-SEQ (when fail finding seeds on the UNITIG) in a read
 //	uint16_t minimizer_hit_pos_alt[ALT_HIT_POS_MAX];//(when hit) the minimizer hit position number on the alt string index, when hit is over 149, set to 149
 //	uint64_t ALT_seed_UNIQUE_N;//(when hit) total number of minimizers on the ALT-SEQ that have unique seed(when fail finding seeds on the UNITIG) in a read
 //
 //	//for simulation reads
 //	uint64_t UNITUG_right_hit_seed_N;
 //	uint64_t UNITUG_wrong_hit_seed_N;
 //
 //	uint64_t ALT_right_hit_seed_N;
 //	uint64_t ALT_wrong_hit_seed_N;
 //
 //	//for SIMU reads
 //	//all count:
 //	uint64_t ALL_right_hit_seed_N;
 //	uint64_t ALL_wrong_hit_seed_N;
 //	uint64_t ALL_No_hit_seed_N;
 
	 //
	 uint64_t BWA_hit_N;
	 int wrong_base_idx;
	 uint8_t wrong_base_qual;
	 uint8_t stop_reason;
 
	 void clear(){
		 memset(this, 0, sizeof(DEBUG_minimizer_evaluate_per_read));
	 }
 
	 void show_data(){
 //		fprintf(stderr, "minimizer_N %ld\n", minimizer_N);
 //		fprintf(stderr, "UNITIG_seed_N %ld\n", UNITIG_seed_N);
 //		fprintf(stderr, "UNITIG_seed_N_ALT %ld\n", UNITIG_seed_N_ALT);
 //		fprintf(stderr, "ALT_with_seed_minimizer_N %ld\n", ALT_with_seed_minimizer_N);
 //		fprintf(stderr, "ALT_seed_N %ld\n", ALT_seed_N);
 //		fprintf(stderr, "ALT_seed_UNIQUE_N %ld\n", ALT_seed_UNIQUE_N);
 //		fprintf(stderr, "UNITUG_right_hit_seed_N %ld\n", UNITUG_right_hit_seed_N);
 //		fprintf(stderr, "UNITUG_wrong_hit_seed_N %ld\n", UNITUG_wrong_hit_seed_N);
 //		fprintf(stderr, "ALT_right_hit_seed_N %ld\n", ALT_right_hit_seed_N);
 //		fprintf(stderr, "ALT_wrong_hit_seed_N %ld\n", ALT_wrong_hit_seed_N);
	 }
 };
 
 #define BWA_hit_N_MAX 100
 struct DEBUG_minimizer_evaluate_overall{
 //	uint64_t minimizer_N[150];//total minimizer number overall
 //	uint64_t UNITIG_seed_N[150];//total seed of on the UNITIG overall
 //	//float UNITIG_seed_rate;//the probability of a minimizer find a seed on the UNITIG
 //	uint64_t UNITIG_seed_N_ALT[150];//total number minimizer of read that search in the alt string overall
 //	//float ALT_seed_rate;//the probability of a minimizer find a seed on the ALT-SEQ (when fail finding seeds on the UNITIG)
 //	uint64_t ALT_with_seed_minimizer_N[150];//total number of minimizes that found seeds on the ALT-SEQ (when fail finding seeds on the UNITIG) in a read
 //	uint64_t ALT_seed_N[ALT_seed_N_MAX];//total number of seeds on the ALT-SEQ (when fail finding seeds on the UNITIG) in a read
 //
 //	uint64_t minimizer_hit_pos_alt[ALT_HIT_POS_MAX];//the minimizer hit position number on the alt string index, when hit is over 19, set to 19
 //	uint64_t ALT_seed_UNIQUE_N[150];//(when hit) total number of minimizers on the ALT-SEQ that have unique seed(when fail finding seeds on the UNITIG) in a read
 //	//float ALT_seed_UNIQUE_rate;//(when hit) the probability of a minimizer find a seed on the ALT-SEQ that is uniq (when fail finding seeds on the UNITIG)
 //
 //	//for simulation reads
 //	uint64_t UNITUG_right_hit_seed_N[150];
 //	uint64_t UNITUG_wrong_hit_seed_N[150];
 //
 //	uint64_t ALT_right_hit_seed_N[150];
 //	uint64_t ALT_wrong_hit_seed_N[150];
 //
 //	//all count:
 //	uint64_t ALL_right_hit_seed_N[150];
 //	uint64_t ALL_wrong_hit_seed_N[150];
 //	uint64_t ALL_No_hit_seed_N[150];
 
	 uint64_t BWA_hit_N[BWA_hit_N_MAX];
	 uint64_t wrong_base_idx[200];
	 uint64_t wrong_base_qual[200];
 
	 uint64_t stop_reason[10];
 
	 void clear(){
		 memset(this, 0, sizeof(DEBUG_minimizer_evaluate_overall));
	 }
	 void add_signal(DEBUG_minimizer_evaluate_per_read &s){
 
		 s.BWA_hit_N = MIN(s.BWA_hit_N, BWA_hit_N_MAX -1 );
		 BWA_hit_N[s.BWA_hit_N]++;
 //
 //		if(s.wrong_base_idx == 0){
 //			fprintf(stderr, "XXXXreadNum 1 %ld %ld\n", s.wrong_base_idx, wrong_base_idx[0]);
 //		}
 
		 s.wrong_base_idx = MIN(s.wrong_base_idx, 199);
		 wrong_base_idx[s.wrong_base_idx]++;
 
 
		 s.wrong_base_qual = MIN(s.wrong_base_qual, 199);
		 wrong_base_qual[s.wrong_base_qual]++;
 
		 s.stop_reason = MIN(s.stop_reason, 9);
		 stop_reason[s.stop_reason]++;
 
	 }
	 void output_uint64_t_vector(uint64_t *v, int total_size, const char *title){
		 for(int i = 0; i < total_size; i++){
			 fprintf(stderr, "%s %d : [ %ld ]\n",title, i, v[i]);
		 }
	 }
	 void output_rst(int seeding_level){
 
		 output_uint64_t_vector(BWA_hit_N, BWA_hit_N_MAX, "BWA_hit_N");
 //		output_uint64_t_vector(wrong_base_idx, 200, "wrong_base_idx");
 //		output_uint64_t_vector(wrong_base_qual, 200, "wrong_base_qual");
 //		output_uint64_t_vector(stop_reason, 9, "stop_reason");
 
		 fprintf(stderr, "\n");
	 }
 };
 
 
 struct OCC_rst{
	 uint64_t bwt_bg;
	 uint64_t bwt_ed;
	 //debug code
	 BWT_MAP_rst map_r;
 };
 
 
 struct CHANGE_detail{
	 uint8_t pos;
	 uint8_t length;
	 uint8_t type;
	 std::vector<uint8_t> change_base;
	 CHANGE_detail(){
		 pos = 0;
		 length = 0;
		 type = 0;
		 change_base.clear();
 
	 }
 
 
	 void set(uint8_t p, uint8_t l, uint8_t t, std::vector<uint8_t> &change_b){
		 pos = p;
		 length = l;
		 type = t;
		 change_base = change_b;
 
	 }
 };
 
 
 struct Read_BWT{
	 uint64_t bwt_k;
	 char seq[151];
	 int l;
	 Read_BWT(){
		 bwt_k = 0;
		 l = 0;
		 memset(seq, 0 , 151);
	 }
 
	 Read_BWT(uint64_t k_, const char *s, int l_){
		 bwt_k = k_;
		 l = l_;
		 int i = 0;
		 int k = l+1;
		 if(l_ == 0){
			 k = 0;
		 }
		 
		 memset(seq, 0 , 151);
		 // if(k<150&&k>30)
		 // {
		 // 	fprintf(stderr, "de");
		 // }
		 for(;k <150; k++){
			 seq[i]=s[k];
			 i++;
		 }
		 seq[i] = '\0';
	 }
 
		 // 从文件读取单个 SAM_OUT 对象
	 friend std::istream& operator>>(std::istream& is, Read_BWT& obj) {
		 is.read(reinterpret_cast<char*>(&obj.bwt_k), sizeof(obj.bwt_k));
		 // 读取其他成员...
		 return is;
	 }
 
	 // 写入单个 SAM_OUT 对象到文件
	 friend std::ostream& operator<<(std::ostream& os, const Read_BWT& obj) {
		 os.write(reinterpret_cast<const char*>(&obj.bwt_k), sizeof(obj.bwt_k));
		 
		 // 写入其他成员...
		 return os;
	 }
 };
 
 
 struct SAM_OUT{
	 uint64_t bwt_k;
	 uint8_t length;
	 kseq_t * read;
	 uint8_t left_length;
 
 
	 //OCC_rst* occ_rst;
	 std::string sam_line;
	  kseq_t *read_p;
	 int sam_flag = 0;
	 int mapped =0;
	 int exc = 0;
 
 
	 SAM_OUT(){
		 //occ_rst = (OCC_rst *)xcalloc(1,sizeof(OCC_rst));
		 bwt_k = 0;
		 length = 0;
		 read = NULL;
		 
		 sam_flag = 0;
		 mapped = 0;
 
	 }
	 void destory(){
		 //bam_destroy1(b);
		 //free(occ_rst);
	 }
 
	 void set_compact_seq(uint64_t bwt_k_, uint8_t l, kseq_t * read_){
		 bwt_k = bwt_k_;
		 length = l;
		 read = read_;
	 }
 
 
	 void set_mate(SAM_OUT *s){
 
		 if(s->sam_flag &BAM_FUNMAP){
			 sam_flag |= BAM_FMUNMAP;
		 }
		 if(s->sam_flag &BAM_FREVERSE){
			 sam_flag |= BAM_FMREVERSE;
		 }
 
 
	 }
 
	 void set_read(kseq_t * read){
		 read_p = read;
	 }
 
	 // 从文件读取单个 SAM_OUT 对象
	 friend std::istream& operator>>(std::istream& is, SAM_OUT& obj) {
		 is.read(reinterpret_cast<char*>(&obj.bwt_k), sizeof(obj.bwt_k));
		 // 读取其他成员...
		 return is;
	 }
 
	 // 写入单个 SAM_OUT 对象到文件
	 friend std::ostream& operator<<(std::ostream& os, const SAM_OUT& obj) {
		 os.write(reinterpret_cast<const char*>(&obj.bwt_k), sizeof(obj.bwt_k));
		 
		 // 写入其他成员...
		 return os;
	 }
 
 
 
	 
	  void write_line(int s, uint flag, std::string  &cigar_p, std::string &md_p, kseq_t * read,
			   std::string &chr_id_p, uint32_t ref_p, int nm , int as,  std::string &mate_chr_id_p, uint32_t mate_ref_p, int rev, int read_length, int exc_);
 
 
 
	 void seq_bin(char  seq, ubyte_t seq_bin){
 
			 switch(seq){
			 case 'A': seq_bin = 0; break;
			 case 'C': seq_bin = 1; break;
			 case 'G': seq_bin = 2; break;
			 case 'T':  seq_bin = 3; break;
			 default :
					 seq_bin = 4; break;
 
		 }
	 }
 
	 char bin_seq(ubyte_t seq_bin){
		 char seq;
		 switch(seq_bin){
			 case '0': seq = 'A'; break;
			 case '1': seq = 'C'; break;
			 case '2': seq = 'G'; break;
			 case '3': seq = 'T'; break;
 
			 }
		 return seq;
	 }
 
 
 };
 
 
 #include <iostream>
 #include <fstream>
 #include <cstring>
 #include <cstdint>
 
 struct Compress_block {
	 uint32_t window_id;
	 uint16_t hap_id;
	 uint16_t hap_offset;
	 uint8_t head_n;
	 uint8_t tail_n;
	 char qname[50];
	 uint8_t lqname;
	 int max_match;
	 std::string quality_score1;
	 uint8_t *qual;
	 uint32_t real_pos;
	 std::string real_seq1;
	 uint16_t flags_1;
	 std::string diff_seq1;
	 std::string diff_base1;
	 int best_begin1;
 
	 uint32_t window_id2;
	 uint16_t hap_id2;
	 uint16_t hap_offset2;
	 uint8_t head_n2;
	 uint8_t tail_n2;
	 char qname2[50];
	 uint8_t lqname2;
	 int max_match2;
	 std::string quality_score2;
	 uint8_t *qual2;
	 uint32_t real_pos2;
	 std::string real_seq2;
	 uint16_t flags_2;
	 std::string diff_seq2;
	 std::string diff_base2;
	 int best_begin2;
	 uint8_t is_paired;
 
	 int hap_offset_flag1,hap_offset_flag2;
 
	 int error_flag = 0;
 
	 Compress_block() {
		 window_id = 0;
		 hap_id = 0;
		 hap_offset = 0;
		 head_n = 0;
		 tail_n = 0;
		 lqname = 0;
		 max_match = 0;
		 real_pos = 0;
		 // head_seq= (char *)xcalloc(1, 150);
		 // tail_seq= (char *)xcalloc(1, 150);
		 // qname = (char *)xcalloc(1, 50);
		 memset(qname, 0, 50);
		 flags_1 = 0;
		 window_id2 = 0;
		 hap_id2 = 0;
		 hap_offset2 = 0;
		 head_n2 = 0;
		 tail_n2 = 0;
		 lqname2 = 0;
		 max_match2 = 0;
		 real_pos2 = 0;
		 // head_seq2= (char *)xcalloc(1, 150);
		 // tail_seq2= (char *)xcalloc(1, 150);
		 // qname2 = (char *)xcalloc(1, 50);
		 memset(qname2, 0, 50);
		 flags_2 = 0;
		 is_paired = 1;
		 error_flag = 0;
	 }
	 void destory(){
		 hap_offset = 0;
		 window_id = 0;
		 hap_offset = 0;
		 head_n = 0;
		 tail_n = 0;
	 }
 
	 void set(uint32_t win_id, uint16_t hap, uint16_t hap_off,  char* hs,  char* ts, int hn, int tn, int ml, char* name, int ql,
		 uint32_t win_id2, uint16_t hap2, uint16_t hap_off2,  char* hs2,  char* ts2, int hn2, int tn2, int ml2, char* name2, int ql2) {
		 window_id = win_id;
		 hap_id = hap;
		 hap_offset = hap_off;
		 head_n = hn;
		 tail_n =tn;
		 memcpy(qname, name, ql);
		 max_match = ml;
		 window_id2 = win_id2;
		 hap_id2 = hap2;
		 hap_offset2 = hap_off2;
		 head_n2 = hn2;
		 tail_n2 =tn2;
		 memcpy(qname2, name2, ql2);
		 max_match2 = ml2;
 
	 }
 
	 void setqname( char* name, int ql, char* name2, int ql2) {
		 memcpy(qname, name, ql);
		 memcpy(qname2, name2, ql2);
	 }
 
	 /*
	 void set_read1(uint64_t win_id, uint64_t hap, uint64_t hap_off,  char* hs,  char* ts, int hn, int tn, int ml, char* name, int ql, char* quality_score, uint8_t* qual_1, uint64_t pos){
		 memset(head_seq, 0, 150);
		 memset(tail_seq, 0, 150);
		 window_id = win_id;
		 hap_id = hap;
		 hap_offset = hap_off;
		 head_n = hn;
		 tail_n =tn;
		 memcpy(head_seq, hs, hn);
		 memcpy(tail_seq, ts, tn);
		 memcpy(qname, name, ql);
		 max_match = ml;
		 memcpy(quality_score1, quality_score, 150);
		 qual = qual_1;
		 real_pos = pos;
	 }
 
	 void set_read2(uint64_t win_id2, uint64_t hap2, uint64_t hap_off2,  char* hs2,  char* ts2, int hn2, int tn2, int ml2, char* name2, int ql2, char* quality_score, uint8_t* qual_2, uint64_t pos2) {
		 memset(head_seq2, 0, 150);
		 memset(tail_seq2, 0, 150);
		 window_id2 = win_id2;
		 hap_id2 = hap2;
		 hap_offset2 = hap_off2;
		 head_n2 = hn2;
		 tail_n2 =tn2;
		 memcpy(head_seq2, hs2, hn2);
		 memcpy(tail_seq2, ts2, tn2);
		 memcpy(qname2, name2, ql2);
		 max_match2 = ml2;
 
		 memcpy(quality_score2, quality_score, 150);
		 qual2 = qual_2;
		 real_pos2 = pos2;
	 }
	 */
 };
 
 
 struct Compress_block_store{
	 uint32_t hap_dev;
	 uint16_t hap_offset;
	 std::string hapdev_num_ACGT;
	 std::string hapoffset_num_ACGT;
	 char head_seq[405]; 
	 char tail_seq[405];  
	 uint8_t head_n;
	 uint8_t tail_n;
	 char qname[50];
	 int lqname;
 
	 uint32_t read_pair_hapid_dev;
	 uint16_t hap_offset2;
	 std::string hapdev_num_ACGT2;
	 std::string hapoffset_num_ACGT2;
	 char head_seq2[405]; 
	 char tail_seq2[405];  
	 uint8_t head_n2;
	 uint8_t tail_n2;
	 char qname2[50];
	 int lqname2;
 
	 bool real_add_flag;
 
	 bool add_flag;
	 int hap_offset_flag1;
	 int hap_offset_flag2;
 
	 std::string real_seq1;
	 std::string real_seq2;
 
 
	 Compress_block_store() {
		 hap_dev = 0;
		 hap_offset = 0;
		 memset(head_seq, 0, sizeof(head_seq));
		 memset(tail_seq, 0, sizeof(tail_seq));
 
		 read_pair_hapid_dev = 0;
		 hap_offset2 = 0;
		 memset(head_seq2, 0, sizeof(head_seq));
		 memset(tail_seq2, 0, sizeof(tail_seq));
		 hap_offset_flag1=0,hap_offset_flag2=0;
	 }
 
	 friend std::ostream& operator<<(std::ostream& os, const Compress_block_store& cb) {
		 os << cb.hapdev_num_ACGT << " "
			<< cb.hap_offset << " "
			<< cb.hapoffset_num_ACGT << " "
			<< cb.head_seq << " "
			<< cb.tail_seq;
		 return os;
	 }
 
	 // 将数字转换为ACGT字符串
	 static void numberToACGT(uint32_t num, std::string & result) {
		 result.clear();
		 if(num==0)
			 result+="A";
		 size_t bit_length = 0;
		 uint32_t numc = num;
		 while (numc > 0) {
			 bit_length += 2;
			 numc >>= 2; // 每次右移两位
		 }
 
		 const char nucleotideMap[4] = {'A', 'C', 'G', 'T'};
		 result.reserve(bit_length / 2); // 提前分配所需空间以避免多次分配
 
		 for (size_t i = 0; i < bit_length; i += 2) {
			 // 取出最低两位
			 uint64_t bits = (num >> (bit_length - 2 - i)) & 0b11;
			 // 根据最低两位的值，映射到相应的核苷酸字符
			 result += nucleotideMap[bits];
		 }
 
	 }
 
	 void store(uint32_t hap_dev_, uint32_t read_pair_hapid_dev_, const Compress_block &cb){
		 
		 hapdev_num_ACGT.clear();
		 hapoffset_num_ACGT.clear();
		 hap_dev =  hap_dev_;
		 hap_offset = cb.hap_offset;
		 numberToACGT(hap_dev_, hapdev_num_ACGT);
		 numberToACGT(cb.hap_offset, hapoffset_num_ACGT);
 
 
		 hapdev_num_ACGT2.clear();
		 hapoffset_num_ACGT2.clear();
		 read_pair_hapid_dev =  read_pair_hapid_dev_;
		 hap_offset2 = cb.hap_offset2;
		 numberToACGT(read_pair_hapid_dev, hapdev_num_ACGT2);
		 numberToACGT(cb.hap_offset2, hapoffset_num_ACGT2);
 
 
	 }
 };
 
 struct Compress_bio_string_block{
	 std::string bio_string;
	 char qname[50];
	 int lqname;
 
	 uint32_t window_id;
	 uint16_t hap_id;
	 uint16_t hap_offset;
	 uint8_t hn;
 
	 uint32_t window_id2;
	 uint16_t hap_id2;
	 uint16_t hap_offset2;
	 uint8_t hn2;
 
	 uint64_t hapid1;
	 uint64_t hapid2;
 
	 std::string real_seq1;
	 std::string real_seq2;
 
 };
 
 struct Compress_final_block{
	 std::string huffman_code;
	 std::string real_bio_string;
	 char qname[50];
	 int lqname;
 
	 uint32_t window_id;
	 uint16_t hap_id;
	 uint16_t hap_offset;
	 uint8_t hn;
 
	 uint32_t window_id2;
	 uint16_t hap_id2;
	 uint16_t hap_offset2;
	 uint8_t hn2;
 
	 uint64_t hapid1;
	 uint64_t hapid2;
 
	 std::string real_seq1;
	 std::string real_seq2;
 
 
 };
 
 struct Compress_final_block_with_huffman_table{
	 std::vector<Compress_final_block> Final_blocks;
	 std::unordered_map<char, std::string> huffmanCode;
 };
 struct HuffmanNode {
	 char data;
	 uint64_t freq;
	 std::shared_ptr<HuffmanNode> left, right;
	 
 
	 HuffmanNode() 
			 : data('\0'),    // 用空字符作为默认数据
			   freq(0),       // 频率默认为0
			   left(nullptr), // 左子节点默认为空
			   right(nullptr) // 右子节点默认为空
		 {}
	 HuffmanNode(char data, unsigned freq) 
		 : data(data), freq(freq), left(nullptr), right(nullptr) {}
 };
 
 // 比较器（用于优先队列）
 struct Compare {
	 bool operator()(const std::shared_ptr<HuffmanNode>& a,
					 const std::shared_ptr<HuffmanNode>& b) {
		 return a->freq > b->freq;  // 小顶堆
	 }
 };
 
 class HuffmanCoder {
 public:
	 // 构建哈夫曼树
	 std::shared_ptr<HuffmanNode> buildHuffmanTree(const std::unordered_map<char, uint64_t>& freqMap) {
		 std::priority_queue<std::shared_ptr<HuffmanNode>, 
					   std::vector<std::shared_ptr<HuffmanNode>>, 
					   Compare> minHeap;
 
		 // 创建叶节点
		 for (const auto& pair : freqMap) {
			 minHeap.push(std::make_shared<HuffmanNode>(pair.first, pair.second));
		 }
 
		 // 处理只有一个字符的情况
		 if (minHeap.size() == 1) {
			 auto left = minHeap.top();
			 minHeap.pop();
			 auto top = std::make_shared<HuffmanNode>('\0', left->freq);
			 top->left = left;
			 minHeap.push(top);
		 }
 
		 // 合并节点直到只剩一个
		 while (minHeap.size() > 1) {
			 auto left = minHeap.top(); minHeap.pop();
			 auto right = minHeap.top(); minHeap.pop();
 
			 auto top = std::make_shared<HuffmanNode>(
				 '\0', left->freq + right->freq);
			 top->left = left;
			 top->right = right;
			 minHeap.push(top);
		 }
 
		 return minHeap.empty() ? nullptr : minHeap.top();
	 }
 
	 // 生成编码表
	 void buildCodes(const std::shared_ptr<HuffmanNode>& root, 
					std::string code,
					std::unordered_map<char, std::string>& huffmanCode) {
		 if (!root) return;
 
		 // 叶节点包含字符
		 if (root->data != '\0') {
			 huffmanCode[root->data] = code.empty() ? "0" : code;
		 }
 
		 buildCodes(root->left, code + "0", huffmanCode);
		 buildCodes(root->right, code + "1", huffmanCode);
	 }
 
	 // 编码字符串
	 std::string encode(const std::string& text, 
				  std::unordered_map<char, std::string>& huffmanCode) {
		 std::string encoded;
		 for (char c : text) {
			 encoded += huffmanCode[c];
		 }
		 return encoded;
	 }
 
	 // 解码字符串
	 std::string decode(const std::string& encoded, 
				  const std::shared_ptr<HuffmanNode>& root) {
		 std::string decoded;
		 auto current = root;
		 
		 for (char bit : encoded) {
			 current = (bit == '0') ? current->left : current->right;
 
			 if (current->data != '\0') {
				 decoded += current->data;
				 current = root;  // 重置到根节点
			 }
		 }
		 return decoded;
	 }
	 std::shared_ptr<HuffmanNode> rebuildTreeFromCodeTable(const std::unordered_map<char, std::string>& code_table) {
		 if (code_table.empty()) {
			 std::cerr << "HuffmanCoder::rebuildTreeFromCodeTable: 编码表为空！" << std::endl;
			 return nullptr;
		 }
 
		 auto root = std::make_shared<HuffmanNode>(); // 根节点
		 root->data = '\0'; // 内部节点数据为空
 
		 for (const auto& pair : code_table) {
			 char character = pair.first;
			 const std::string& code = pair.second;
 
			 auto current_node = root;
			 for (char bit : code) {
				 if (bit == '0') {
					 if (!current_node->left) {
						 current_node->left = std::make_shared<HuffmanNode>();
						 current_node->left->data = '\0';
					 }
					 current_node = current_node->left;
				 } else { // bit == '1'
					 if (!current_node->right) {
						 current_node->right = std::make_shared<HuffmanNode>();
						 current_node->right->data = '\0';
					 }
					 current_node = current_node->right;
				 }
			 }
			 // 在路径的末端节点设置字符
			 current_node->data = character;
		 }
 
		 return root;
	 }
 };
 
 struct CompareSAM_OUT {
	 bool operator()(const std::pair<Read_BWT, int>& a, const std::pair<Read_BWT, int>& b) {
		 return a.first.bwt_k > b.first.bwt_k; // 最小堆
	 }
 };
 struct OCC_rst1{
	 uint8_t bwt_bg[5];
	 uint8_t bwt_ed[5];
	 void set(OCC_rst &r){
		 uint64_t bg = r.bwt_bg;
		 bwt_bg[0] = (bg>>32);
		 bwt_bg[1] = (bg>>24);
		 bwt_bg[2] = (bg>>16);
		 bwt_bg[3] = (bg>>8);
		 bwt_bg[4] = (bg);
 
		 uint64_t ed = r.bwt_ed;
		 bwt_ed[0] = (ed>>32);
		 bwt_ed[1] = (ed>>24);
		 bwt_ed[2] = (ed>>16);
		 bwt_ed[3] = (ed>>8);
		 bwt_ed[4] = (ed);
	 }
	 void restore(OCC_rst &r){
		 uint64_t bg = 0;
		 bg += (bwt_bg[0]); bg <<= 8;
		 bg += (bwt_bg[1]); bg <<= 8;
		 bg += (bwt_bg[2]); bg <<= 8;
		 bg += (bwt_bg[3]); bg <<= 8;
		 bg += (bwt_bg[4]);
		 r.bwt_bg = bg;
 
		 uint64_t ed = 0;
		 ed += (bwt_ed[0]); ed <<= 8;
		 ed += (bwt_ed[1]); ed <<= 8;
		 ed += (bwt_ed[2]); ed <<= 8;
		 ed += (bwt_ed[3]); ed <<= 8;
		 ed += (bwt_ed[4]);
		 r.bwt_ed = ed;
	 }
 };
 
 struct MAP_rst{
	 std::vector<uint16_t> change_detail_list;
	 int score;
	 uint8_t r[6];
	 void set(int32_t chr_ID, int32_t POS){
		 uint64_t store = ((uint64_t)chr_ID << 22) + (POS & 0x3fffff);
		 r[5] = store >> 40;
		 r[4] = store >> 32;
		 r[3] = store >> 24;
		 r[2] = store >> 16;
		 r[1] = store >> 8;
		 r[0] = store;
	 }
	 void restore(int32_t &chr_ID, int32_t &POS){
		 uint64_t store = 0;
 
		 store += (r[5]); store <<= 8;
		 store += (r[4]); store <<= 8;
		 store += (r[3]); store <<= 8;
		 store += (r[2]); store <<= 8;
		 store += (r[1]); store <<= 8;
		 store += (r[0]);
		 chr_ID = store >> 22;
		 POS = store & 0x3fffff;
	 }
 };
 
 //data used for PIPELINE each threads: 3 threads
 struct CLASSIFY_THREAD_DATA
 {
	 bam1_t* bams1;
	 bam1_t* bams2;
	 // kseq_t	*	seqs1; // first  read in a pair
	 // kseq_t	*	seqs2; // second read in a pair
	 //OUTPUT_RST * rst; 	//mapping result for read1 and read2
	 // std::vector<OCC_rst>  sam_rst_r1; //sam result for read1
	 // std::vector<OCC_rst>  sam_rst_r2; //sam result for read2
	 std::vector<Compress_block> compress_block_r;
 
	 std::vector<SAM_OUT>  sam_out_r1; //sam result for read1
	 std::vector<SAM_OUT>  sam_out_r2; //sam result for read2
 
	 char * fq_buff1 = NULL;
	 char * fq_buff2 = NULL;// unmap reads buff
 
	 char * sam_buff = NULL;
 
 
	 std::vector<DEBUG_minimizer_evaluate_per_read>  debug_eval_r1; //
	 std::vector<DEBUG_minimizer_evaluate_per_read>  debug_eval_r2; //
 
	 int 		readNum;
	 void 	*	share_data_pointer;// register for the shared data
 
 };
 
 struct CLASSIFY_SHARE_DATA
 {
	 //idx used for alignment
	 //BWT_IDX_loader *idx;
	 Aln_online::IDX_loader *idx;
	 //shared
	 htsFile *bam_in_fp;//output sam
	 
	 // kstream_t 	*_fp1 = NULL; //fastq files 1
	 // kstream_t 	*_fp2 = NULL; //fastq files 2
	 FILE *sam_out_fp;//output sam
	 FILE *unmapped_fq1_file;
	 FILE *unmapped_fq2_file;
 
	 bam_hdr_t *header;
 
	 OL_PAR 	*o; 		//mapping option
	 Classify_buff_pool 	*buff = NULL;//data used for each classify thread
	 CLASSIFY_THREAD_DATA *data = NULL;//for each pipeline thread
 
	 DEBUG_minimizer_evaluate_overall debug_eval_all;
 
	 std::vector<OCC_rst1> all_occ_rst; //list to store final OCC
	 std::vector<MAP_rst> all_map_r; //list to store final mapping rst
	 bam1_t* input_bam_buff;
	 int input_bam_buff_num;
	 std::vector<FILE*> blockFiles;
 
 };
 
 struct BWT_CLASSIFY_MAIN{
	 CLASSIFY_SHARE_DATA *share = NULL;
	 void init_run(int argc, char *argv[]);
 private:
	 void map_r_sort(std::vector<MAP_rst> & r);
	 static void *classify_pipeline(void *shared, int step, int tid, void *_data);									//pipeline
	 static int load_reads(htsFile *_fp1, bam1_t *_bams1, bam1_t *_bams2, int n_needed, void *shared);	//function in step 1
	 static void inline worker_for(void *_data, long data_index, int thread_index); 									//function in step 2
	 static void output_results(FILE *fp,int readNum, std::vector<Compress_block > &cb, double ctime, std::vector<FILE*> &blockFiles, int read_length); //function in step 3
 };
 
 struct DECOMPRESS_MAIN
 {
	 CLASSIFY_SHARE_DATA *share = NULL;
	 void run(int argc, char *argv[]);
 
 };
 
 struct reverse_inf
 {
	 char qname[50];
 
	 uint64_t window_id;
	 uint64_t hap_id;
	 uint64_t hap_offset;
	 char head_seq[150];
	 char tail_seq[150];
 
	 uint64_t window_id2;
	 uint64_t hap_id2;
	 uint64_t hap_offset2;
	 char head_seq2[150];
	 char tail_seq2[150];
 };
 
 struct reverse_inf_blocks
 {
	 std::vector<reverse_inf> blocks;
	 int block_num;
 };
 
 }
 
 #endif
 