/*
 * SA_process.hpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_PROCESS2_SA_HPP_
#define SRC_BWT_ONLINE_PROCESS2_SA_HPP_

#include "BWT_idx/bwt_cpt_sa.hpp"
#include "aln_para.hpp"
#include "aln_RST_DEF.hpp"
#include "num_def.hpp"
namespace Aln_online{

//the result for each 1M block of SA execute threads
class SA_result{

public:

	std::vector<Single_read_aln_rst_t1_48_exact> exact_rst_list;
	std::vector<Single_read_aln_rst_t2_48_gap> gap_NM_ONE_rst_list;
	std::vector<Single_read_aln_rst_t3_80_gap> gap_NM_SMALL_rst_list;
	std::vector<Single_read_aln_rst_t4_160_gap> gap_NM_MIDDLE_rst_list;
	std::vector<Single_read_aln_rst_t5_cpx> gap_NM_BIG_rst_list;

	void init(){
		exact_rst_list.clear();
		gap_NM_ONE_rst_list.clear();
		gap_NM_SMALL_rst_list.clear();
		gap_NM_MIDDLE_rst_list.clear();
		gap_NM_BIG_rst_list.clear();
	}

	void store_result_exact_alignmet(int32_t wb_ID, int aln_position_in_WB){
		exact_rst_list.emplace_back();
		exact_rst_list.back().set(wb_ID, aln_position_in_WB);
	}

	void store_result_gap_alignmet_NM_UNIQ(int32_t wb_ID, int read_hamming_mapping_position, int16_t change_detail){//, MM_idx_loader * wb_idx){

			gap_NM_ONE_rst_list.emplace_back();
			gap_NM_ONE_rst_list.back().set(wb_ID, read_hamming_mapping_position, change_detail);
//
//			//debug:
//			Single_read_aln_rst_unpack p;
//			p.clear();
//			result_buff_single_sample_in_1M_block[buff_ID].gap_NM_ONE_rst_list.back().unpack(p);
//			p.generate_read_string(wb_idx, 150);
//			p.print_detail();
		}

	void store_result_gap_alignmet_NM_SMALL(int32_t wb_ID, int read_hamming_mapping_position,
			std::vector<uint16_t> &change_detail){//, MM_idx_loader * wb_idx){
		gap_NM_SMALL_rst_list.emplace_back();
		gap_NM_SMALL_rst_list.back().set( wb_ID, read_hamming_mapping_position, change_detail);

//		//debug:
//		Single_read_aln_rst_unpack p;
//		p.clear();
//		result_buff_single_sample_in_1M_block[buff_ID].gap_NM_SMALL_rst_list.back().unpack(p);
//		p.generate_read_string(wb_idx, 150);
//		p.print_detail();
	}

	void store_result_gap_alignmet_NM_MIDDLE(int32_t wb_ID, int read_hamming_mapping_position,
			std::vector<uint16_t> &change_detail){//, MM_idx_loader * wb_idx){
		gap_NM_MIDDLE_rst_list.emplace_back();
		gap_NM_MIDDLE_rst_list.back().set( wb_ID, read_hamming_mapping_position, change_detail);

//		//debug:
//		Single_read_aln_rst_unpack p;
//		p.clear();
//		result_buff_single_sample_in_1M_block[buff_ID].gap_NM_SMALL_rst_list.back().unpack(p);
//		p.generate_read_string(wb_idx, 150);
//		p.print_detail();
	}

	void store_result_gap_alignmet_NM_BIG(int32_t wb_ID, int aln_position_in_WB, uint8_t* read_string, int read_len){//, MM_idx_loader * wb_idx){

		gap_NM_BIG_rst_list.emplace_back();
		gap_NM_BIG_rst_list.back().set(wb_ID, aln_position_in_WB, read_string, read_len);

//		//debug:
//		Single_read_aln_rst_unpack p;
//		p.clear();
//		result_buff_single_sample_in_1M_block[buff_ID].gap_NM_BIG_rst_list.back().unpack(p);
//		p.generate_read_string(wb_idx, 150);
//		p.print_detail();
	}
//
};


#define EDIT_MAX_DIS_SMALL 4
#define EDIT_MAX_DIS 11
#define BAND_SIZE 5
#define ERROR_SCORE -10000

#define EDIT_DIS_D_DIAG 0
#define EDIT_DIS_D_UP 1
#define EDIT_DIS_D_LEFT 2
#define EDIT_DIS_D_DIAG_SNP 3

class SIMPLE_EDIT_DISTANCE_handler{
	//BASIC:
	int read_len;
	int **M;
	int **D;
	//INPUT:
	//OUTPUT:
	int max_s;
public:
	//OUTPUT:
	std::vector<uint16_t>change_detail_list;
	int max_pos_j;
	int ref_begin_pos_adjust;//skip some base in the begin of ref

//	void print(FILE * o){
//		fprintf(max_pos_j, );
//	}

	void set_edit_distance_error(){
		max_s = ERROR_SCORE;
	}

	void set_edit_distance_zero(){
		max_s = 0;
	}

	bool edit_dis_over_max(){
		return max_s < -EDIT_MAX_DIS;
	}

	int get_edit_distance(){
		return -max_s;
	}
	void init(int read_len){
		this->read_len = read_len;
		// M = (int **)xmalloc(sizeof(int *)*(read_len + 1));
		// D = (int **)xmalloc(sizeof(int *)*(read_len + 1));
		// for (int i = 0; i < read_len + 1; i++){
		// 	M[i] = (int *)xmalloc(sizeof(int)*(read_len + 1));
		// 	D[i] = (int *)xmalloc(sizeof(int)*(read_len + 1));
		// }
		M = (int **)xcalloc((read_len + 1), sizeof(int *));
		D = (int **)xcalloc((read_len + 1), sizeof(int *));
		for (int i = 0; i < read_len + 1; i++){
			M[i] = (int *)xcalloc(read_len+1, sizeof(int));
			D[i] = (int *)xcalloc(read_len +1, sizeof(int));
		}
		//init: S1
		M[0][0] = 0;	//init matrix
		D[0][0] = 0;	//init matrix
		//init: S2
		for (int i = 1; i < BAND_SIZE + 1; i++){
			M[i][0] = M[i-1][0] + (-1);
			D[i][0] = EDIT_DIS_D_UP;
		}
		for (int j = 1; j < BAND_SIZE + 1; j++){
			M[0][j] = M[0][j-1] + (-1);
			D[0][j] = EDIT_DIS_D_LEFT;
		}
		//set band
		//init: S3
		for (int i = BAND_SIZE + 1; i < read_len + 1; i++)
			 M[i - BAND_SIZE - 1][i] = ERROR_SCORE;
		for (int i = BAND_SIZE + 1; i < read_len + 1; i++)
			 M[i][i - BAND_SIZE - 1] = ERROR_SCORE;
	}

	void destroy(){
		for (int i = 0; i < read_len + 1; i++){
			free(M[i]);
			free(D[i]);
		}
		free(M);
		free(D);
	}

	void show_M(){
		for (int i = 0; i < read_len + 1; i++){
			for (int j = 0; j < read_len + 1; j++){
				fprintf(stderr, "%10d ", M[i][j]);
			}
			fprintf(stderr, "\n ");
		}
	}

	void show_D(){
		for (int i = 0; i < read_len + 1; i++){
			for (int j = 0; j < read_len + 1; j++){
				fprintf(stderr, "%10d ", D[i][j]);
			}
			fprintf(stderr, "\n ");
		}
	}

	void test1(){
		//debug code:
		//read 79
		uint8_t q_o[150] = {0, 2, 2, 0, 2, 3, 1, 0, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 2, 0, 2, 2, 0, 1, 0, 2, 2, 0, 0, 2, 0, 0, 1, 0, 2, 2, 0, 2, 2, 1, 3, 1, 0, 0, 3, 2, 0, 0, 0, 1, 3, 3, 3, 0, 1, 1, 0, 2, 2, 0, 1, 0, 1, 1, 3, 1, 3, 2, 0, 0, 2, 3, 2, 3, 2, 2, 0, 0, 1, 3, 2, 3, 2, 3, 3, 3, 1, 1, 0, 2, 2, 0, 0, 0, 2, 1, 1, 0, 2, 0, 0, 0, 3, 1, 3, 2, 3, 3, 2, 2, 2, 2, 3, 3, 1, 1, 3, 3, 3, 1, 0, 2, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 3, 3, 2, 5, 5, 5,};
		uint8_t t_o[150] = {0, 2, 2, 0, 2, 3, 1, 0, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 2, 0, 2, 2, 0, 1, 0, 2, 2, 0, 0, 2, 0, 0, 1, 0, 2, 2, 0, 2, 2, 1, 3, 1, 0, 0, 3, 2, 0, 0, 0, 1, 3, 3, 3, 0, 1, 1, 0, 2, 2, 0, 1, 0, 1, 1, 3, 1, 3, 2, 0, 0, 2, 3, 1, 3, 2, 2, 0, 0, 1, 3, 2, 3, 2, 3, 3, 3, 1, 1, 0, 2, 2, 0, 0, 0, 2, 1, 1, 0, 2, 0, 0, 0, 3, 1, 3, 2, 3, 3, 2, 2, 2, 2, 3, 3, 1, 1, 3, 3, 3, 1, 0, 2, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 3, 3, 1, 1, 3,};
		bool is_reverse = false;
		int q_len = 150;
		int t_len = 150;
		uint8_t *q = q_o;// + q_len - 1;
		uint8_t *t = t_o;// + t_len - 1;
		if(is_reverse){
			q = q_o + q_len - 1;
			t = t_o + t_len - 1;
		}
		//uint8_t *q = q_o + q_len - 1;
		//uint8_t *t = t_o + t_len - 1;
		count_the_edit_distance(q, t, t_len, is_reverse,true);
		//show_M();
		//show_D();
	}

	static void string_change(char * from, uint8_t *to, int len){
		for(int i = 0; i < len; i++){
			switch(from[i]){
			case 'A': to[i] = 0; break;
			case 'C': to[i] = 1; break;
			case 'G': to[i] = 2; break;
			case 'T': to[i] = 3; break;
			case 'N': to[i] = 5; break;
			}
		}
	}

//	void test4(MM_idx_loader *wb_idx){
//			//debug code:
//			//read 79
//			char q_c[151] = "GGAAAGGCATTTGAGAATTTTCTAATGTTTCGTCATTTGAAAAAAAAAAAAATGACTTTTGGCTACAACTTAATGCAACCAGTGACTTATTTCCCAACTCCCTATACCCCTCCTTGAATCTTGTATTCCTAAATTTAAAATTTCANNNNN";
//			char t_c[151] = "GGAAAGGCATTTGAGAATTTTCTAATGTTTCCATTTGAAAAAAAAAAAAATGACTTTTGGCTACAACTTAATGCAACCTGTGACTTATCTCCCAACTCCCTATACCCCTCCTTGAATCTTGTATTCCTAAATTTAAAATTTCAGGCTAAA";
//
//			uint8_t q_o[150];
//			uint8_t t_o[150];
//
//			bool is_reverse = false;
//			int q_len = 150;
//			int t_len = 150;
//			uint8_t *q = q_o;// + q_len - 1;
//			uint8_t *t = t_o;// + t_len - 1;
//			string_change(q_c, q_o, q_len);
//			string_change(t_c, t_o, t_len);
//			//uint8_t *q = q_o + q_len - 1;
//			//uint8_t *t = t_o + t_len - 1;
//			count_the_edit_distance(q, t, t_len, is_reverse);
//			//show_M();
//			//show_D();
//			//test the final results
//			Single_read_aln_rst_t3_80_gap r;
//			//r.set(546870, 20, change_detail_list);
//
//			Single_read_aln_rst_unpack p;
//			p.clear();
//			r.unpack(p, q_len);
//			p.generate_read_string(wb_idx, read_len);
//			p.print_detail(stderr);
//		}

	void test2(){
		//read 169
		uint8_t q[150] = {2, 3, 1, 0, 1, 1, 1, 1, 0, 1, 3, 2, 1, 0, 1, 3, 1, 3, 0, 2, 1, 1, 3, 2, 2, 2, 1, 0, 0, 1, 0, 2, 0, 2, 3, 2, 0, 2, 0, 2, 1, 1, 3, 2, 3, 1, 3, 3, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 3, 3, 0, 1, 0, 0, 0, 0, 2, 1, 0, 2, 2, 1, 0, 1, 0, 0, 1, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 2, 2, 3, 1, 1, 1, 0, 1, 0, 3, 1, 0, 3, 2, 2, 1, 2, 0, 0, 1, 3, 0, 0, 0, 2, 1, 2, 2, 2, 2, 0, 3, 3, 3, 0, 3, 3, 1, 3, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,};
		uint8_t t[150] = {2, 3, 1, 0, 1, 1, 1, 1, 0, 1, 3, 2, 1, 0, 1, 3, 1, 3, 0, 2, 1, 1, 3, 2, 2, 2, 1, 0, 0, 1, 0, 2, 0, 2, 3, 2, 0, 2, 0, 2, 1, 1, 3, 2, 3, 1, 3, 3, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 3, 3, 0, 1, 0, 0, 0, 0, 2, 1, 0, 2, 2, 1, 0, 1, 0, 0, 1, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 2, 2, 3, 1, 1, 1, 0, 1, 0, 3, 1, 0, 3, 2, 2, 1, 0, 0, 0, 1, 3, 0, 0, 0, 2, 1, 2, 2, 2, 2, 0, 3, 3, 3, 0, 3, 3, 1, 3, 2, 3, 1, 1, 1, 3, 0, 3, 1, 1, 3, 1, 2, 2, 0, 2, 3, 1, 0, 2, 2, 1,};

		bool is_reverse = false;
		int q_len = 150;
		int t_len = 150;
		count_the_edit_distance(q, t, t_len, is_reverse, true);
		//show_M();
		//show_D();
	}

	void test3(){
		//debug code:
		uint8_t q_o[150] = {0, 0, 0, 1, 3, 0, 0, 0, 2, 1, 2, 2, 2, 2, 0, 3, 3, 3, 0, 3, 3, 1, 3, 2, 3, 1, 1, 1, 3, 0, 3, 1, 1, 3, 1, 2, 2, 0, 2, 3, 3, 0, 2, 2, 1, 3, 0, 2, 0, 0, 1, 0, 0, 1, 3, 0, 3, 0, 2, 0, 3, 2, 3, 1, 3, 1, 0, 0, 2, 0, 2, 2, 3, 2, 2, 0, 1, 2, 0, 1, 3, 3, 3, 1, 3, 2, 3, 2, 3, 3, 1, 0, 3, 3, 2, 3, 1, 0, 3, 1, 3, 3, 0, 1, 3, 0, 2, 3, 2, 3, 2, 2, 1, 1, 2, 3, 3, 2, 0, 0, 2, 2, 0, 0, 1, 3, 1, 3, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, };
		uint8_t t_o[150] = {1, 0, 0, 0, 1, 3, 0, 0, 0, 2, 1, 2, 2, 2, 2, 0, 3, 3, 3, 0, 3, 3, 1, 3, 2, 3, 1, 1, 1, 3, 0, 3, 1, 1, 3, 1, 2, 2, 0, 2, 3, 1, 0, 2, 2, 1, 3, 0, 2, 0, 0, 1, 0, 0, 3, 0, 3, 0, 2, 0, 3, 2, 3, 1, 3, 1, 0, 0, 2, 0, 2, 2, 3, 2, 2, 0, 1, 2, 0, 1, 3, 3, 3, 1, 3, 2, 3, 2, 3, 3, 1, 0, 3, 3, 2, 3, 1, 0, 3, 1, 3, 3, 0, 1, 3, 0, 2, 3, 2, 3, 2, 2, 1, 1, 2, 3, 3, 2, 0, 0, 2, 2, 0, 0, 1, 3, 1, 3, 2, 2, 1, 3, 3, 1, 1, 0, 2, 3, 3, 1, 1, 0, 2, 2, 0, 0, 1, 3, 1, 1, };
		bool is_reverse = true;
		int q_len = 150;
		int t_len = 150;
		uint8_t *q = q_o + q_len - 1;
		uint8_t *t = t_o + t_len - 1;
		count_the_edit_distance(q, t, t_len, is_reverse, true);
		//simple_edit_distance(q, q_len, t, t_len, is_reverse);
		//show_M();
		//show_D();
	}

//describe of change detail;
//      000:0 INS: A
//      001:1 INS: C
//      010:2 INS: G
//      011:3 INS: T
//      100:4 DEL
//      101:5 SNP: + 1
//      102:6 SNP: + 2
//      103:7 SNP: + 3
	static uint16_t get_change_INS(int q_pos, uint8_t change_to_base){
		uint16_t change_detail = change_to_base;
		xassert(q_pos < 0xff, "MAX read length is 254");
		change_detail += (q_pos << 3);
		return change_detail;
	}

	static uint16_t get_change_DEL(int q_pos){
		uint16_t change_detail = 4;
		xassert(q_pos < 0xff, "MAX read length is 254");
		change_detail += (q_pos << 3);
		return change_detail;
	}
	//the query is to; the target is from
	static uint16_t get_change_SNP(int q_pos, uint8_t change_to_base, uint8_t change_from_base){
		uint16_t change_detail = 0;
		if(change_from_base < change_to_base){
			change_detail = (change_to_base - change_from_base + 4);
		}else if(change_from_base > change_to_base){
			change_detail = (4 + change_to_base - change_from_base + 4);
		}else{
			xassert(0, "ERROR, SNP without change.");
		}
		xassert(q_pos < 0xff, "MAX read pos is 254");
		change_detail += (q_pos << 3);
		return change_detail;
	}

	void simple_edit_distance(uint8_t *q, int q_len, uint8_t *t, int t_len)
	{
//Answer:
// t: TAAAAG
// q: T_AAAG
//		M _ T A A A G
//		_ 0 1 2 3 4 5
//		T 1 0 1 2 3 4
//		A 2 1 0 1 2 3
//		A 3 2 1 0 1 2
//		A 4 3 2 1 0 1
//		A 5 4 3 2 1 1
//		G 6 5 4 3 2 1
//
//		D _ T A A A G
//		_ _ L L L L L
//		T U D L L L L
//		A U U D D/L D/L L
//		A U U D/U D D/L L
//		A U U D/U D/U D L
//		A U U D/U D/U D/U D
//		G U U U U U D
//
		//alignment
		for (int i = 1; i < q_len + 1; i++) {
			for (int j = i - BAND_SIZE ; j < i + BAND_SIZE + 1; j++) {
				if(j <= 0) continue;
				if(j > t_len) continue;
				bool is_match = (q[i-1] == t[j-1]);
				int scoreDiag = M[i-1][j-1] + ((is_match)?0:(-1));
				int scoreUp   = M[i-1][j] - 1;
				int scoreLeft = M[i][j-1] - 1;
				M[i][j] = MAX(scoreDiag, scoreUp);
				M[i][j] = MAX(M[i][j], scoreLeft);
				if(M[i][j] == scoreDiag){
					if(is_match)
						D[i][j] = EDIT_DIS_D_DIAG;
					else
						D[i][j] = EDIT_DIS_D_DIAG_SNP;
				}else if(M[i][j] == scoreUp){
					D[i][j] = EDIT_DIS_D_UP;
				}else{
					D[i][j] = EDIT_DIS_D_LEFT;
				}
			}
		}

		//search the max score when reach the end of query
		max_s = ERROR_SCORE; max_pos_j = 0;
		for (int j = q_len - BAND_SIZE; j < q_len + BAND_SIZE + 1 && j < t_len + 1; j++){
			if(max_s < M[q_len][j]){
				max_s = M[q_len][j];
				max_pos_j = j;
			}
		}

		//search reverse
		//trace back and generate change detail:
		int i_pos = q_len;
		int j_pos = max_pos_j;
		int read_not_enough = 0;
		change_detail_list.clear();

		while(j_pos >= 1){
			if(i_pos < 1){//read is not enough
				if(false) fprintf(stderr, "i_pos %d , j_pos %d \t @ Read not enough CLIP\n",i_pos, j_pos);
				read_not_enough ++;
				j_pos --;
				continue;
			}
			//WARNING:: The position for INS and SNP is i_pos - 1; the position for DEL is i_pos.
			//DONT change the number.
			switch(D[i_pos][j_pos]){
			case EDIT_DIS_D_DIAG:
				i_pos --; j_pos --;
				break;
			case EDIT_DIS_D_UP: //INS
				change_detail_list.emplace_back(get_change_INS(i_pos - 1, q[i_pos - 1]));
				if(false) fprintf(stderr, "i_pos %d , j_pos %d \t with base %c @ INS\n",i_pos - 1, j_pos - 1, "ACGT"[q[i_pos - 1]]);
				i_pos --;
				break;
			case EDIT_DIS_D_LEFT: //DEL
				change_detail_list.emplace_back(get_change_DEL(i_pos));
				if(false) fprintf(stderr, "i_pos %d , j_pos %d \t @ DEL\n",i_pos -  1, j_pos - 1);
				j_pos --;
				break;
			case EDIT_DIS_D_DIAG_SNP:
				change_detail_list.emplace_back(get_change_SNP(i_pos - 1, q[i_pos - 1], t[j_pos - 1]));
				if(false) fprintf(stderr, "i_pos %d , j_pos %d \t with base %c @ SNP\n",i_pos - 1, j_pos - 1, "ACGT"[q[i_pos - 1]]);
				i_pos --; j_pos --;
				break;
			}
		}

		if(false) fprintf(stderr, "i_pos %d j_pos %d \n", i_pos,j_pos);

		//std::reverse(change_detail_list.begin(), change_detail_list.end());
		max_s += read_not_enough;
		//Note: when clipping alignment of ref occurred, adjust the final reference mapping position
		ref_begin_pos_adjust += read_not_enough + (j_pos - i_pos);
	}

	//P1: this function search the hamming distance between the read and a SUB-STRING of ref
	//return the position with most small distance and the distance
	//besides, when the most small distance is over MAX_DIS, it return (MAX_DIS+1) and the position is set to -1.
	//P2: about the storing of change detail:
	//	each change_detail occupied 12bit, the 12 bit stored in the lower 12 bit of uint16_t,
	//which is :
		//2 bit, check bit. set to 0 when the original base is C or T, set to 1 when original base is A or G.
		//8 bit: the position of read(max length of read is 256);
		//2 bit: base change to;
	//for example, when REF is A and ALT is C (A ---->> C), the check bit is 00 and the changed-to base is 01
	/*input stored in binary: ACGT is 0/1/2/3*/
	int count_the_edit_distance(uint8_t *q, uint8_t * t, int t_len, bool is_reverse, bool show_log){
		//simple get distance
		change_detail_list.clear();
		max_pos_j = 0;
		max_s = 0;
		int skip_base = 0;
		// {
		// 	ref_begin_pos_adjust = 0;
		// 	for(int i = 0; i < read_len; i++)
		// 		if(q[i] >= 4)
		// 			skip_base++;
		// 	if(is_reverse){
		// 		q += skip_base;
		// 		t += skip_base;
		// 		ref_begin_pos_adjust += skip_base;
		// 	}
		// }
		//remove the clip N of reads; some base of the reads is UNKNOWN
		int cur_read_len = read_len - skip_base;
		int cur_ref_len = t_len - skip_base;
		//store in reverse
		
		// for(int i = cur_read_len - 1; i >= 0; i--){
		// 	if(q[i]>=4) continue;
		// 	uint8_t q_k = q[i];
		// 	if(is_reverse){
		// 		int k = cur_read_len-1-i;
		// 		// add by zhanganqi
		// 		if(q[k]>=4){
		// 			continue;
		// 		}
		// 		q_k = 3-q[k];
				
		// 	}
		// 	if(q[i] < 5 && t[i] != q_k){
		// 		//fprintf(stderr, "aaaaaaaaaaaaaa%d,%d\n",i,k);
		// 		max_s --;
		// 		if(max_s < -EDIT_MAX_DIS)
		// 			break;
		// 		change_detail_list.emplace_back(get_change_SNP(i, q_k, t[i]));
		// 	}
		// }
		for(int i = cur_read_len - 1; i >= 0; i--){
			//zhanganqi q[i] 5？
			if(q[i] < 4 && t[i] != q[i]){
				//fprintf(stderr, "aaaaaaaaaaaaaa%d,%d\n",i,k);
				max_s --;
				if(max_s < -EDIT_MAX_DIS)
					break;
				change_detail_list.emplace_back(get_change_SNP(i, q[i], t[i]));
			}
		}
		 if(max_s < -EDIT_MAX_DIS){
			return 0;
		 }
		 return 1;
		//when EDIT distance is too high
		// if(max_s < -EDIT_MAX_DIS){
		// 	if(false){
		// 		fprintf(stderr, "hamming_distance %d\n", max_s);
		// 		for(int i = 0; i < read_len; i++) {fprintf(stderr, "%d, ", q[i]);}	fprintf(stderr, "\n");
		// 		for(int i = 0; i < t_len; i++)    {fprintf(stderr, "%d, ", t[i]);}	fprintf(stderr, "\n");
		// 	}
		// 	simple_edit_distance(q, cur_read_len, t, cur_ref_len);
		// 	if(false && max_s < -EDIT_MAX_DIS){
		// 		fprintf(stderr, "EDIT_distance %d\n", max_s);
		// 	}
		// }
	}
};

}
#endif /* SRC_BWT_ONLINE_PROCESS2_SA_HPP_ */
