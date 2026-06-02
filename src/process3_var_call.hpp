/*
 * var_call.hpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

#ifndef SRC_BWT_ONLINE_PROCESS3_VAR_CALL_HPP_
#define SRC_BWT_ONLINE_PROCESS3_VAR_CALL_HPP_



//for WB and POS
#define POS_OFFSET 22
#define POS_MASK 0x3fffff
#define WB_POS_MAX ((POS_MASK)-0x3ff)
#define REF_POS_MAX (POS_MASK)
#include "aln_para.hpp"
namespace Aln_online{

//***************************************************ALN buff pool ****************************************************//

//***************************************************SA process buff pool ****************************************************//

//***************************************************Pipeline 3：variant calling buff pool ****************************************************//

#define READ_RST_TYPE_E1 1
#define READ_RST_TYPE_G2 2
#define READ_RST_TYPE_G3 3
#define READ_RST_TYPE_G4 4
#define READ_RST_TYPE_G5 5
class WB_CALLER{

private:
	//
	int window_basic_bg;
	int window_ed;
	//index:
	ALN_ONLINE::MM_idx_loader * wb_idx;
	std::vector<Window_t> *window_info;
	int gap_aln_default_read_len;
	//output
	Window_block_counter_Known_Var *exact_counter;
	Nearby_WB_Signal_handler_Novel *gap_counter;
	Var_generater_all_type *var_generater;

	int known_var_ID;
	int novel_var_ID;

	int chip_id;
	int sample_id;

	OL_PAR *p;
	FILE*log_f;
public:
	void init(ALN_ONLINE::Simple_ref_handler *ref_h, FILE*vcf_out_log, uint64_t ave_read_depth,
			FILE * vcf_out, ALN_ONLINE::MM_idx_loader * wb_idx,	std::vector<Window_t> *window_info, int global_sample_ID, OL_PAR *p){
		this->p = p;
		this->log_f = vcf_out_log;
		//sample_id = get_original_sample_id(chip_id, global_sample_ID);
		sample_id =1;
		char *sample_name_string = (char *)xmalloc(1024);
		sprintf(sample_name_string, "C%d_S%d", chip_id, sample_id);

		window_basic_bg = 0;
		window_ed = 0;

		this->wb_idx = wb_idx;
		this->window_info = window_info;
		//for gap alignment, the default read length is true read length minus the seed length;
		//because we don`t know what`s the true base in the last 25 bits
		int gap_aln_seed_length = 25;
		gap_aln_default_read_len = 150 - gap_aln_seed_length;//todo::

		//exact counter
		exact_counter = (Window_block_counter_Known_Var*) calloc(1, sizeof(Window_block_counter_Known_Var));
		exact_counter->init(*window_info, vcf_out_log);
		exact_counter->ref_h = ref_h;
		//gap counter
		gap_counter = (Nearby_WB_Signal_handler_Novel*) calloc(1, sizeof(Nearby_WB_Signal_handler_Novel));
		gap_counter->init(this->wb_idx, ref_h, vcf_out_log, p->show_assembly);
		//
		var_generater = (Var_generater_all_type*) calloc(1, sizeof(Var_generater_all_type));
		var_generater->init(ave_read_depth, vcf_out, vcf_out_log, ref_h, sample_name_string);

		known_var_ID = 0;
		novel_var_ID = 0;
	}
	void show_statistics(){
		gap_counter->show_statistics();
		exact_counter->show_statistics();
	}

	void destroy(){
		exact_counter->destroy();
		if(exact_counter != NULL){ free(exact_counter); exact_counter = NULL;}
		gap_counter->destroy();
		if(gap_counter != NULL){ free(gap_counter); gap_counter = NULL;}
		var_generater->destroy();
		if(var_generater != NULL){ free(var_generater); var_generater = NULL;}
	}
private:
	//input
	//for exact t1, global
	int e_t1_i;// = 0;
	int e_t1_total;
	Single_read_aln_rst_t1_48_exact * e_t1_l;
	//for exact t1, last data
	int e_t1_last_wb_idx;
	int e_t1_last_pos;

	//for gap t2, global
	int g_t2_i;// = 0;
	int g_t2_total;
	Single_read_aln_rst_t2_48_gap * g_t2_l;
	//for gap t2, last data
	Single_read_aln_rst_unpack g_t2_last;

	//for gap t3, global
	int g_t3_i;// = 0;
	int g_t3_total;
	Single_read_aln_rst_t3_80_gap * g_t3_l;
	//for gap t3, last data
	Single_read_aln_rst_unpack g_t3_last;

	//for gap t4, global
	int g_t4_i;// = 0;
	int g_t4_total;
	Single_read_aln_rst_t4_160_gap * g_t4_l;
	//for gap t3, last data
	Single_read_aln_rst_unpack g_t4_last;


	//for gap t4, global
	int g_t5_i;// = 0;
	int g_t5_total;
	Single_read_aln_rst_t5_cpx * g_t5_l;
	//for gap t4, last data
	Single_read_aln_rst_unpack g_t5_last;

public:
	void init_input(
			int e_t1_total, Single_read_aln_rst_t1_48_exact * e_t1_l,
			int g_t2_total, Single_read_aln_rst_t2_48_gap * g_t2_l,
			int g_t3_total, Single_read_aln_rst_t3_80_gap * g_t3_l,
			int g_t4_total, Single_read_aln_rst_t4_160_gap * g_t4_l,
			int g_t5_total, Single_read_aln_rst_t5_cpx * g_t5_l){
		e_t1_i = 0;
		this->e_t1_total = e_t1_total;
		this->e_t1_l = e_t1_l;
		//for exact t1, last data
		e_t1_last_wb_idx = -1;
		e_t1_last_pos = 0;

		//for gap t2, global
		g_t2_i = 0;
		this->g_t2_total = g_t2_total;
		this->g_t2_l = g_t2_l;
		g_t2_last.wb_ID = -1;

		//for gap t3, global
		g_t3_i = 0;
		this->g_t3_total = g_t3_total;
		this->g_t3_l = g_t3_l;
		g_t3_last.wb_ID = -1;

		//for gap t4, global
		g_t4_i = 0;
		this->g_t4_total = g_t4_total;
		this->g_t4_l = g_t4_l;
		g_t4_last.wb_ID = -1;

		//for gap t4, global
		g_t5_i = 0;
		this->g_t5_total = g_t5_total;
		this->g_t5_l = g_t5_l;
		g_t5_last.wb_ID = -1;
	}

private:
	void add_signal_e_t1(int window_id, int window_basic, int pos){
		//for reference signals
		if(pos >= WB_POS_MAX){
			pos -= WB_POS_MAX;
			pos += (*window_info)[window_id].total_hap_length - 300;//the reads aligned to the reference is stored after all other haplotype
		}
		int32_t st_pos = pos;
		int32_t ed_pos = st_pos + 150;
		//window_info[window_id].with_in_repeat_check(st_pos, ed_pos);
		//if(ed_pos > st_pos + 16){//skip too short reads
		//debug code:
		if(false){
			fprintf(stderr, "add_signal_e_t1  window_id %d ,window_basic %d, pos %d \n", window_id, window_basic, pos);
		}
		exact_counter->add_signal(window_id - window_basic, pos, ed_pos - st_pos);
	}

	void add_all_signal_exact_rst(){
		//try last:
		if(e_t1_last_wb_idx >= window_ed){
			return;
		}else{
			if(e_t1_last_wb_idx != -1){
				xassert(e_t1_last_wb_idx >= window_basic_bg, "");
				add_signal_e_t1(e_t1_last_wb_idx, window_basic_bg, e_t1_last_pos);
				e_t1_last_wb_idx = -1;//clear the result after use the data
			}
		}
		for(; e_t1_i < e_t1_total;){
			e_t1_l[e_t1_i++].restore(e_t1_last_wb_idx, e_t1_last_pos);
			//new a window block counter
			if (e_t1_last_wb_idx >= window_ed){
				break;
			}
			add_signal_e_t1(e_t1_last_wb_idx, window_basic_bg, e_t1_last_pos);
			e_t1_last_wb_idx = -1;//clear the result after use the data
		}
	}

	void add_all_signal_g_rst(int store_type){
		//type:
		Single_read_aln_rst_unpack * g_r = &g_t2_last;
		int *g_tn_i = &g_t2_i;
		int g_tn_total = g_t2_total;
		switch(store_type){
			case READ_RST_TYPE_G3:g_r = &g_t3_last; g_tn_i = &g_t3_i;g_tn_total = g_t3_total; break;
			case READ_RST_TYPE_G4:g_r = &g_t4_last; g_tn_i = &g_t4_i;g_tn_total = g_t4_total; break;
			case READ_RST_TYPE_G5:g_r = &g_t5_last; g_tn_i = &g_t5_i;g_tn_total = g_t5_total; break;
		}

		//try last:
		if(g_r->wb_ID >= window_ed){
			return;
		}else{
			if(g_r->wb_ID != -1){
				gap_counter->add_signal(*g_r, window_basic_bg);
				g_r->wb_ID = -1;//clear the result after use the data
			}
		}
		while((*g_tn_i) < g_tn_total){
			g_r->clear();
			switch(store_type){
			case 2: g_t2_l[*g_tn_i].unpack(*g_r, gap_aln_default_read_len); break;
			case 3: g_t3_l[*g_tn_i].unpack(*g_r, gap_aln_default_read_len); break;
			case 4: g_t4_l[*g_tn_i].unpack(*g_r, gap_aln_default_read_len); break;
			case 5: g_t5_l[*g_tn_i].unpack(*g_r); break;
			}
			(*g_tn_i)++;
			//	p.generate_read_string(wb_idx, read_len);
			//	p.print_detail();
			//new a window block counter
			if (g_r->wb_ID >= window_ed){
				break;
			}
			gap_counter->add_signal(*g_r, window_basic_bg);
			g_r->wb_ID = -1;//clear the result after use the data
		}
	}

	void set_wb_region(	int window_basic_bg){
		this->window_basic_bg = window_basic_bg;
		this->window_ed = window_basic_bg + 4;
	}

	void add_all_signal_all_type(){
		//adding signals for exact results
		add_all_signal_exact_rst();
		//adding signals for gap results
		gap_counter->clear(window_basic_bg);
		add_all_signal_g_rst(READ_RST_TYPE_G2);
		add_all_signal_g_rst(READ_RST_TYPE_G3);
		add_all_signal_g_rst(READ_RST_TYPE_G4);
		add_all_signal_g_rst(READ_RST_TYPE_G5);
	}

public:
	void get_signals_and_call_variants(int window_basic){
		//S1: set work region
		set_wb_region(window_basic);
		exact_counter->set_new_wb(*window_info, window_basic);
		//S2: get signals
		add_all_signal_all_type();
		if(window_basic < 10)//skip the first WB, only adding signals
			return;
		//S3: generate candidates
		//calling:
		//init 4 window block per time, not only one, because wb ID is not sorted
		//call variant for all novel variants
		gap_counter->generate_novel_vars_candidate_in_4_nearby_WB();
		// call the KNOWN variant and get the final results
		exact_counter->generate_known_var_candidate(*window_info);
		//S4: call variants
		{
			uint32_t st_wb_ID = window_basic_bg;
			for(int i = 0; i < 4; i++){
				//debug:
				bool with_var = var_generater->vcf_generatig(
						*exact_counter->get_pos_depth_by_id(i - 1), *exact_counter->get_pos_depth_by_id(i),
						*window_info, st_wb_ID + i - 1, st_wb_ID + i,
						*exact_counter->get_process_block_by_id(i - 1), *exact_counter->get_process_block_by_id(i),
						&(gap_counter->result_4[i]) );
				if(with_var)
					fprintf(log_f, "Above VAR(s) is from: WB %d@ [%d:%d~%d]\n", st_wb_ID + i, window_info[0][st_wb_ID + i].chr_ID, window_info[0][st_wb_ID + i].st_pos, window_info[0][st_wb_ID + i].st_pos + 150 - 1);
			}
		}

		{//after process, store "cur" to "pre"
			gap_counter->store_cur_to_old();
			exact_counter->final_process_after_var_calling();
		}
	}

	bool all_read_is_processed(){
		return (e_t1_i >= e_t1_total &&
				g_t2_i >= g_t2_total &&
				g_t3_i >= g_t3_total &&
				g_t4_i >= g_t4_total );
	}
};




#endif /* SRC_BWT_ONLINE_PROCESS3_VAR_CALL_HPP_ */
