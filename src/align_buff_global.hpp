/*
 * align_buff_global.hpp
 *
 *  Created on: 2022年5月9日
 *      Author: fenghe
 */

#ifndef ALIGN_BUFF_GLOBAL_HPP_
#define ALIGN_BUFF_GLOBAL_HPP_


class ALIGN_BUFF_global{

public:
  //数据项：
//  TYPE1 *data1;
//  TYPE2 *data2;
//  ......
//  TYPEN *dataN;


  //初始化函数
  //如果该函数需要参数，联系我提供接口
  void init(int thread_idx_number){
//      data1 = malloc();
//      .....
//      .....
  }

  //销毁函数（可选，不写也行）
  void destroy(){
//      free(data1);

  }
};

#endif /* ALIGN_BUFF_GLOBAL_HPP_ */
