/*
 * var_map.cpp
 *
 *  Created on: 2023年2月20日
 *      Author: fenghe
 */

#include "var_map.hpp"
//int MAX_LOW:: when the depth is <= MAX_LOW, the VAR will be set as LOW
void get_GT_type(uint var_depth, uint  ref_depth, int &GT_TYPE, int &PASS_TYPE, int MAX_LOW){
	GT_TYPE = 0;// 0: 0/0; 1: 0/1; 2: 1/1
	PASS_TYPE = 1; // 0: NOT_PASS; 1: PASS; 2:LOW_COV
	//basic condition: low VAR depth;
	if((int)var_depth < MAX_LOW){
		PASS_TYPE = 2; //LOW_COV
	}

	if(ref_depth > 8 * var_depth){
		PASS_TYPE = 2;
		GT_TYPE = 0;
	}else if(ref_depth == 1 || ref_depth * 8 < var_depth)
		GT_TYPE = 2;
	else
		GT_TYPE = 1;
}



