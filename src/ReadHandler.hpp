/*
 * BamHandler.hpp
 *
 *  Created on: 2020年4月28日
 *      Author: fenghe
 */

#ifndef BAMHANDLER_HPP_
#define BAMHANDLER_HPP_

extern "C"
{
	#include "clib/bam_file.h"
	#include "clib/vcf_file.h"
extern uint64_t kmerMask[33];
}

#include"RefHandler.hpp"
#include<vector>
#include<algorithm>

struct READ_record{
	READ_record(uint16_t flag_, uint8_t direction_, uint64_t position_, int read_l_, int cigar_l_,
			int read_index_, int quality_index_, int cigar_index_,	int soft_left_, int soft_right_, int NM_NUM_) :
			flag(flag_), direction(direction_), position(position_), read_l(read_l_), cigar_l(cigar_l_),
			read_index(read_index_), quality_index(quality_index_), cigar_index(cigar_index_),
			soft_left(soft_left_), soft_right(soft_right_), NM_NUM(NM_NUM_) {
	}
	//basic
	uint16_t 	flag;
	uint8_t 	direction; // the direction of the read, and will be direction of mate when reads are unmapped
	int		 	position; // the position of the read, and will be position of mate when reads are unmapped
	int 		read_l;
	int 		cigar_l;
	int 	 	read_index;
	int 		quality_index;
	int  		cigar_index;
	int soft_left;
	int soft_right;
	int NM_NUM;

	//std::sort(l.begin(), l.end(), SVE::cmp_by_position);
	static inline int cmp_position(const READ_record &a, const READ_record &b){
		return a.position < b.position;
	}

	void show(){
		fprintf(stderr, "[%d %d %d]\t", flag, direction, position);
		fprintf(stderr, "[%d %d %d]\t", read_l, cigar_l, read_index);
		fprintf(stderr, "[%d %d]\t", quality_index, cigar_index);
		fprintf(stderr, "[%d %d %d]\n", soft_left, soft_right, NM_NUM);
	}
};

namespace Read_type{
	enum T {SR, DR, UM, TL, unknown};
}
//used to store trans-location read pairs:
struct trans_Read_item{
	int32_t tid;
	int32_t mtid;
	int32_t pos;
	int32_t mpos;
	uint8_t flag;

	trans_Read_item(	int32_t tid_, int32_t mtid_, int32_t pos_, int32_t mpos_, uint8_t flag_){ tid = tid_, mtid = mtid_, pos = pos_, mpos = mpos_, flag = flag_; }
	trans_Read_item(	bam1_t *br){ tid = br->core.tid, mtid = br->core.mtid, pos = br->core.pos, mpos = br->core.mpos, flag = br->core.flag; }
	//copy:
	trans_Read_item(const trans_Read_item &b){ memcpy(this, &b, sizeof(trans_Read_item)); }
	//sort mtid:
	static inline int cmp_by_mtid(const trans_Read_item &a, const trans_Read_item &b){
		if(a.mtid == b.mtid)	return a.mpos < b.mpos;
		else							return a.mtid < b.mtid;
	}

	void prinf(FILE * log){
		fprintf(stderr, "[cur:%d:%d ; mate: %d:%d]\n", tid, pos, mtid, mpos );
	}
};


struct BAM_handler{
public:
	int chr_ID;
	int pos;
	int window_id;
	uint localoffset;
	//typedef std::vector<READ_record> READ_LIST;
	//void init(char * bamFileName, char *realignment_bamFileName, char *TL_read_Filename, char * ref_file_name){
	void init(){
		chr_ID = 0;
		pos = 0;
		window_id = 0;
		localoffset = 0;
		// const int MAX_READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
		// sr_read_index = new int [MAX_READ_INDEX_SIZE];
		// um_read_index = new int [MAX_READ_INDEX_SIZE];
		// storeReadBuff = new uint8_t[5000];
		// kv_init(base);
		// kv_init(quality);
		kv_init(cigar);
		//bam_file_open(bamFileName, ref_file_name, NULL, &file);
	}

	void destroy(){

		// delete[]sr_read_index;
		// delete[]um_read_index;
		// delete[]storeReadBuff;
		// bam_file_close(&file);
		// if(re_alignment_file._hdr != NULL)
		// 	bam_file_close(&re_alignment_file);
	}
	bool clip_low_quality_Filter(bam1_t *br, int read_len, int soft_left, int soft_right);
	bool storeReadSR(bam1_t *br, int soft_left, int soft_right, int gap_mismatch_inside);// store normal SR signals
	void storeReadUM(bam1_t *br, uint8_t *query);
	//void storeReadCore(	Read_type::T t, int readLen, uint8_t *seq, uint8_t * qual,	int cigarLen,
//			uint32_t* bam_cigar, uint16_t flag, int32_t pos,	int soft_left, int soft_right, int gap_mismatch_inside);
	void storeReadCore1(int readLen, uint8_t *seq, uint8_t * qual,	int cigarLen,
			uint32_t* bam_cigar, uint16_t flag, int32_t pos);
	void storeReadCore(bam1_t *br);
	// void getReadStr(uint8_t * to, READ_record& c_r, int bg, int len){
	// 	memcpy(to, base.a + c_r.read_index + bg, len);
	// }
	// inline uint8_t* getReadStr(READ_record& c_r){ return base.a + c_r.read_index;}
	// inline uint8_t* getQualStr(READ_record& c_r){ return quality.a + c_r.quality_index;}

	bool isSameHeader(const BAM_handler &B) const;
	static bool pass_compact_filter(uint8_t *s, const int len);
	//bool clip_AAA_TailFilter(uint8_t *query, const int read_len);//when a clip string clip to an "AAAAAAAA..." tail, return true
	bool clip_AAA_TailFilter(uint8_t *read_bin, const int read_len, int soft_left, int soft_right);

	void clear(){
		// read_list.clear();
		// um_read_list.clear();
		chr_ID = 0;
		pos = 0;
		window_id = 0;
		localoffset = 0;
		// base.n = 0;
		// quality.n = 0;
		cigar.n = 0;
		//index_already_built = false;
	}

	// void build_read_index_SINGLE(READ_LIST &l, int *index, int region_st_index){
	// 	for(int i = 0; i < READ_INDEX_SIZE; i++)
	// 		index[i] = -1;
	// 	int read_ID = 0;
	// 	for(auto r: l){
	// 		int c_index = (r.position - region_st_index)/ 100;
	// 		if(c_index >= 0 && c_index < READ_INDEX_SIZE && index[c_index] == -1)
	// 			index[c_index] = read_ID;
	// 		read_ID++;
	// 	}
	// }

	// inline void build_read_index(int region_st_index){
	// 	if(index_already_built) return;
	// 	build_read_index_SINGLE(read_list, sr_read_index, region_st_index);
	// 	build_read_index_SINGLE(um_read_list, um_read_index, region_st_index);
	// 	index_already_built = true;
	// }

	
	// READ_LIST read_list;
	// READ_LIST um_read_list;

	//TransReadLoader tr_loader;
	//std::vector<trans_Read_item> trans_read_list;//used to store trans-location read pairs:

	// Bam_file file;
	// Bam_file TL_file;
	// Bam_file re_alignment_file;

private:
	BAM_handler(const BAM_handler &b);//can`t be copy

	//read lists and index
	//bool index_already_built = false;
	const static int READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
	// int * sr_read_index;
	// int * um_read_index;
	// uint8_t * storeReadBuff;
public:
	//buffs
	// kvec_T(uint8_t, BASE_STR);//store read string
	// kvec_T(uint8_t, QUALITY_STR);

	// BASE_STR base;
	// QUALITY_STR quality;
	path_t cigar;
};

#endif /* BAMHANDLER_HPP_ */
