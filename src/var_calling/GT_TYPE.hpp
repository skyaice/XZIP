/*
 * GT_TYPE.hpp
 *
 *  Created on: 2023年3月20日
 *      Author: fenghe
 */

#ifndef SRC_BWT_IDX_GT_TYPE_HPP_
#define SRC_BWT_IDX_GT_TYPE_HPP_

#include<stdlib.h>

void get_GT_type_KNOWN(uint var_depth, uint  ref_depth, int &GT_TYPE, int &PASS_TYPE, int MAX_LOW);
void get_GT_type_NOVEL(uint var_depth, uint  ref_depth, int &GT_TYPE, int &PASS_TYPE, int MAX_LOW);


#endif /* SRC_BWT_IDX_GT_TYPE_HPP_ */
