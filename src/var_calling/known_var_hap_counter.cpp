/*
 * hap_count.c
 *
 *  Created on: 2022年6月3日
 *      Author: zyx
 */


#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "known_var_hap_counter.hpp"

void Hap_counter::add_helper(int head_tail, uint64_t pos, uint64_t *add_pos, uint64_t *save_pos){
	if((*add_pos & HIGH_DIGIT_CHECK)!=0 ){
		if(head_tail==1){
			for(uint64_t i=0; i<8-pos; i++){
				uint64_t a=0x80;
				uint64_t b=0x01;
				uint64_t check_tmp = (*add_pos & ( a <<(i*8)));
				if(check_tmp!=0){
					uint64_t reset_num  = ( b <<(i*8));
					*save_pos = *save_pos+reset_num;
					*add_pos = (*add_pos - (a << (i*8)));
				}
			}
		}
		else{
			for(uint i=0; i<8-pos; i++){
				uint64_t check_tmp = (*add_pos &(CHECK_HEAD>>(i*8)));
				if(check_tmp!=0){
					uint64_t reset_num  = (CHECK_HEAD>>(i*8));
					*save_pos = *save_pos+reset_num;
					*add_pos = (*add_pos - (CHECK_HEAD >> (i*8)));
				}
			}
		}
	}
	if(head_tail==1){
		*add_pos = (*add_pos)+(ADD_BASE>>(pos*8));
	}else{
		*add_pos = (*add_pos)+(ADD_BASE<<((8-pos)*8));
	}
}
