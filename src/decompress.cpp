#include <getopt.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <sstream>
#include <queue>
#include <inttypes.h>
#include <random>
#include "BWT_aln.hpp"
#include "CPPLIB/tools.hpp"
extern "C"{
#include "clib/kthread.h"
#include "clib/utils.h"
#include "clib/desc.h"
}
#include "occ_RST_DEF.hpp"
#include "./BWT_idx/var_map.hpp"


#include <set>

#include <math.h>

using namespace BWT_aln;


#define PIPELINE_T_NUM 3//one for reading; one for classifying; one for writing
#define STEP_NUM PIPELINE_T_NUM
#define N_NEEDED_B 2000000 //2M read pair per time
#define N_NEEDED N_NEEDED_B+100
#define WINDOW_ID_MAX 20588639
uint32_t hap_id_list[WINDOW_ID_MAX];
uint64_t read_sum = 0;
uint32_t quick_pow_32(int di, int zhi) {
	uint32_t sum = 1;
	while (zhi!= 0) {
		if (zhi & 1) {    //这里是zhi%2==1的意思
			sum *= di;
		}
		zhi >>= 1;    //这里是zhi/=2的意思，即向右移动一位在二进制中即除以2
		di = di * di;
	}
	return sum;
}
void read_hap_id_list()
{
    std::string filename = "/home/user/yexiang/hap_id_list.txt";
    std::ifstream in_file(filename, std::ios::binary);
    in_file.read(reinterpret_cast<char*>(hap_id_list), 
                 WINDOW_ID_MAX * sizeof(uint32_t));
}
uint64_t reverse_dev_num(std::string dev_string)
{
    uint64_t final = 0;
    int length = dev_string.size();
    std::string num_string;
    for(int i = 0; i < length; i++)
    {
        if(dev_string[i] == 'A')num_string += "00";
        else if(dev_string[i] == 'C')num_string += "01";
        else if(dev_string[i] == 'G')num_string += "10";
        else if(dev_string[i] == 'T')num_string += "11";
    }
    int num_length = num_string.size();
    //std::cout<<"num string is"<<num_string<<std::endl;
    for(int i = num_length - 1;i >= 0; i--)
    {
        if(num_string[i] == '1')final += quick_pow_32(2, num_length - i - 1);
    }
    return final;
}
uint64_t reverse_offset_num(std::string offset_string)
{
    //CAATA 0100001100 = 268
    uint64_t final = 0;
    int length = offset_string.size();
    std::string num_string;
    for(int i = 0; i < length; i++)
    {
        if(offset_string[i] == 'A')num_string += "00";
        else if(offset_string[i] == 'C')num_string += "01";
        else if(offset_string[i] == 'G')num_string += "10";
        else if(offset_string[i] == 'T')num_string += "11";
    }
    //std::cout<<"num string is"<<num_string<<std::endl;
    int num_length = num_string.size();
    int base = 1;
    for(int i = num_length - 1;i >= 0; i--)
    {
        if(num_string[i] == '1')final += base;
        base *= 2;
    }
    return final;
}
std::string convertBitsToString(const std::vector<unsigned char>& binary_data, size_t bit_length) {
    std::string result;
    result.reserve(bit_length);
    
    for (size_t i = 0; i < bit_length; ++i) {
        size_t byte_index = i / 8;
        size_t bit_index = 7 - (i % 8);  // 匹配写入时的位序（高位在前）
        
        bool bit_value = (binary_data[byte_index] >> bit_index) & 1;
        result.push_back(bit_value ? '1' : '0');
    }
    
    if (result.length() != bit_length) {
        throw std::runtime_error("比特转换错误: 长度不匹配");
    }
    
    return result;
}
bool read_compressed_block(Compress_final_block_with_huffman_table& blocks_with_huffman_table, std::ifstream& in_file) 
{
     // 读取Final_blocks数量
     
     size_t block_count;
     in_file.read(reinterpret_cast<char*>(&block_count), sizeof(size_t));
     blocks_with_huffman_table.Final_blocks.resize(block_count);
     std::cout << "读取区块数量: " << block_count << std::endl;
 
     // 读取每个压缩块
     for (auto& block : blocks_with_huffman_table.Final_blocks) 
     {
        size_t code_len = 0;
        in_file.read(reinterpret_cast<char*>(&code_len), sizeof(size_t));
        size_t byte_count = (code_len + 7) / 8;
        std::vector<unsigned char> binary_data(byte_count);
        in_file.read(reinterpret_cast<char*>(binary_data.data()), byte_count);
        block.huffman_code = convertBitsToString(binary_data, code_len);
        if (block.huffman_code.size() != code_len) {
            throw std::runtime_error("哈夫曼编码长度不匹配: " + 
                                     std::to_string(block.huffman_code.size()) +
                                     " vs " + std::to_string(code_len));
        }
		// 2. 读取qname元数据
   		 in_file.read(reinterpret_cast<char*>(&block.lqname), sizeof(int));

   		 // 3. 读取qname原始数据
    	 in_file.read(block.qname, block.lqname);

    	// 4. 读取real_bio_string
    	size_t bio_len;
    	in_file.read(reinterpret_cast<char*>(&bio_len), sizeof(size_t));
    	block.real_bio_string.resize(bio_len);
    	in_file.read(&block.real_bio_string[0], bio_len);

    	// 5. 读取8个uint64_t字段
    	in_file.read(reinterpret_cast<char*>(&block.window_id), sizeof(uint64_t));
   	 	in_file.read(reinterpret_cast<char*>(&block.window_id2), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hap_id), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hap_id2), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hap_offset), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hap_offset2), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hapid1), sizeof(uint64_t));
    	in_file.read(reinterpret_cast<char*>(&block.hapid2), sizeof(uint64_t));
     }
     std::cout << "开始读取哈夫曼表" << std::endl;
    
    // 3. 读取哈夫曼表
    // 使用eof()检查是否还有更多数据，但需谨慎处理
    while (in_file.peek() != EOF) {
        char key = in_file.get();
        
        if (in_file.eof()) break; // 检查是否到达文件末尾
        
        uint16_t code_len = 0;
        in_file.read(reinterpret_cast<char*>(&code_len), sizeof(uint16_t));
        
        if (in_file.eof()) break; // 再次检查
        
        std::string code(code_len, '\0');
        in_file.read(&code[0], code_len);
        
        if (in_file.gcount() != static_cast<std::streamsize>(code_len)) {
            throw std::runtime_error("哈夫曼编码读取不完整");
        }
        
        blocks_with_huffman_table.huffmanCode[key] = code;
        std::cout << key << ":" << code << std::endl;
    }
    
     return true;
}
bool load_from_file(const std::string& filename, Compress_final_block_with_huffman_table& block) 
{
    std::ifstream in_file(filename, std::ios::binary);
    if (!in_file) return false;
    // 读取结构体数据
    if (!read_compressed_block(block, in_file)) 
    {
        std::cerr << "读取数据失败" << std::endl;
        return false;
    }

    return true;
}
std::string decode_huffman(const std::string& encoded_str,const std::unordered_map<std::string, char>& reverse_map) 
{
    std::string result;
    std::string current_code;

    for (char bit : encoded_str) {
        current_code += bit;
        
        // 正确用法：检查迭代器有效性
        const auto it = reverse_map.find(current_code);
        if (it != reverse_map.end()) { 
            result += it->second;
            current_code.clear();
        }
    }

    
    if (!current_code.empty())
    {
        throw std::runtime_error("解码失败：存在未匹配的编码前缀 - " + current_code);
    }

    return result;
}
void bio_string_to_num(std::vector<Compress_bio_string_block> bio_string_blocks, uint32_t hap_id_list[WINDOW_ID_MAX], uint32_t *pre_hap_id)
 {
	 std::vector<Compress_block_store> reverse_store_blocks;
	 Compress_block_store reverse_store_block;
	 for (const Compress_bio_string_block& block : bio_string_blocks)//提取bio_string中的相关字符串
	 {
		 
		 int length = block.bio_string.size();
		 int pos = 0;
		 int flag = 1;
		 std::string reverse_string;
		 reverse_string.clear();
		 while(pos < length)
		 {
			 if(pos == length - 1)
			 {
				 memset(reverse_store_block.qname, 0, sizeof(reverse_store_block.qname));
				 if(block.bio_string[pos] == 'A')reverse_store_block.add_flag = 0;
				 else reverse_store_block.add_flag = 1;
				 int length_qname = block.lqname;
				 memcpy(reverse_store_block.qname, block.qname, length_qname);
				 reverse_store_block.lqname = block.lqname;
				 reverse_store_block.real_window_id = block.window_id;
				 reverse_store_block.real_window_id2 = block.window_id2;
				 reverse_store_block.real_hap_id = block.hap_id;
				 reverse_store_block.real_hap_id2 = block.hap_id2;
				 reverse_store_block.real_hap_offset = block.hap_offset;
				 reverse_store_block.real_hap_offset2 = block.hap_offset2;
				 reverse_store_block.real_hapid1 = block.hapid1;
				 reverse_store_block.real_hapid2 = block.hapid2;
				 reverse_store_blocks.push_back(reverse_store_block);
				 break;
			 }
			 if(block.bio_string[pos] == '$')
			 {
				 if(flag == 1)
				 {
					 reverse_store_block.hapdev_num_ACGT = reverse_string;
					 reverse_string.clear();
					 pos++;
					 flag++;
				 }
				 else if(flag == 2)
				 {
					 reverse_store_block.hapoffset_num_ACGT = reverse_string;
					 reverse_string.clear();
					 pos++;
					 flag++;
				 }
				 else if(flag == 3)
				 {
					 reverse_store_block.hapdev_num_ACGT2 = reverse_string;
					 reverse_string.clear();
					 pos++;
					 flag++;
				 }
				 else if(flag == 4)
				 {
					 reverse_store_block.hapoffset_num_ACGT2 = reverse_string;
					 reverse_string.clear();
					 pos++;
					 flag++;
				 }
				 else if(flag == 5)
				 {
					 for(int i = 0; i < reverse_string.length(); i++)
					 {
						reverse_store_block.head_seq[i] = reverse_string[i];
					 }
					 reverse_string.clear();
					 pos++;
					 flag++;
				 }
				 else if(flag == 6)
				 {
					for(int i = 0; i < reverse_string.length(); i++)
					{
					   reverse_store_block.tail_seq[i] = reverse_string[i];
					}
					reverse_string.clear();
					pos++;
					flag++;
				 }
				 else if(flag == 7)
				 {
					for(int i = 0; i < reverse_string.length(); i++)
					{
					   reverse_store_block.head_seq2[i] = reverse_string[i];
					}
					reverse_string.clear();
					pos++;
					flag++;
				 }
				 else if(flag == 8)
				 {
					for(int i = 0; i < reverse_string.length(); i++)
					{
					   reverse_store_block.tail_seq2[i] = reverse_string[i];
					}
					reverse_string.clear();
					pos++;
					flag++;
				 }
			 }
			 else
			 {
				 reverse_string += block.bio_string[pos];
				 pos++;
			 }
		 }
		 
	 }
 
	 for (Compress_block_store& block : reverse_store_blocks)//进一步恢复compress_block_store相关信息
	 {
		 
		 block.hap_dev = reverse_dev_num(block.hapdev_num_ACGT);
		 block.read_pair_hapid_dev = reverse_dev_num(block.hapdev_num_ACGT2);
		 block.hap_offset = reverse_offset_num(block.hapoffset_num_ACGT);
		 block.hap_offset2 = reverse_offset_num(block.hapoffset_num_ACGT2);
	 }
	 
	 //debug
	 
	 int cnt = 0;
	 cnt = 0;
	 for (Compress_block_store& block : reverse_store_blocks)
	 {
		 Compress_block reverse_block;
		 uint32_t hapid1 = *pre_hap_id + block.hap_dev;
		 uint32_t hapid2;
		 if(block.add_flag == 0)hapid2 = hapid1 - block.read_pair_hapid_dev;
		 else hapid2 = hapid1 + block.read_pair_hapid_dev;
		 reverse_block.hap_offset = block.hap_offset;
		 reverse_block.hap_offset2 = block.hap_offset2;
		 *pre_hap_id = hapid1;


		 //二分查找
		 int l = 0 ;
		 int r = WINDOW_ID_MAX;
		 while (l < r)
		 {
			 int mid = l + r >> 1;	//(l+r)/2
			 if (hapid1 < hap_id_list[mid])  r = mid;    // check()判断mid是否满足性质
			 else l = mid + 1;
		 }
		 reverse_block.hap_id = hapid1 - hap_id_list[l - 1];
		 reverse_block.window_id = l - 1;

		 l = 0 ;
		 r = WINDOW_ID_MAX;
		 while (l < r)
		 {
			 int mid = l + r >> 1;	//(l+r)/2
			 if (hapid2 < hap_id_list[mid])  r = mid;    // check()判断mid是否满足性质
			 else l = mid + 1;
		 }
		 reverse_block.hap_id2 = hapid2 - hap_id_list[l - 1];
		 reverse_block.window_id2 = l - 1;

		 memcpy(reverse_block.qname, block.qname, block.lqname);
		 if(cnt < 100)
		 {
			 std::cout<<"-----真实信息-----"<<std::endl;
			 std::cout<<block.qname<<std::endl;
			 std::cout<<block.real_window_id<<std::endl;
			 std::cout<<block.real_window_id2<<std::endl;
			 std::cout<<block.real_hap_id<<std::endl;
			 std::cout<<block.real_hap_id2<<std::endl;
			 std::cout<<block.real_hap_offset<<std::endl;
			 std::cout<<block.real_hap_offset2<<std::endl;
			 std::cout<<block.real_hapid1<<std::endl;
			 std::cout<<block.real_hapid2<<std::endl;

			 std::cout<<"-----解码信息-----"<<std::endl;
			 std::cout<<reverse_block.qname<<std::endl;
			 std::cout<<reverse_block.window_id<<std::endl;
			 std::cout<<reverse_block.window_id2<<std::endl;
			 std::cout<<reverse_block.hap_id<<std::endl;
			 std::cout<<reverse_block.hap_id2<<std::endl;
			 std::cout<<reverse_block.hap_offset<<std::endl;
			 std::cout<<reverse_block.hap_offset2<<std::endl;
			 std::cout<<"-----"<<std::endl;
		 }
		 cnt++;
	 }
 }

bool read_from_bin(const std::string& filename, uint32_t hap_id_list[WINDOW_ID_MAX], uint32_t *pre_hap_id) 
 {
	 Compress_final_block_with_huffman_table loaded;
	 std::vector<Compress_bio_string_block> bio_string_blocks;
	 std::vector<Compress_block> origin_blocks;
	 if (load_from_file(filename, loaded)) 
	 {
		 std::vector<Compress_final_block> blocks = loaded.Final_blocks;
		 std::unordered_map<char, std::string> huffmantable = loaded.huffmanCode;
		 std::unordered_map<std::string, char> reverse_map;
		 for (const auto& pair : huffmantable) {
			 reverse_map[pair.second] = pair.first;
		 }
		 int cnt = 0;
		 for (auto& block : blocks)
		 {
			 Compress_bio_string_block bio_string_block;
			 std::string huffmancode = block.huffman_code;
			 std::string decoded = decode_huffman(block.huffman_code, reverse_map);
			 bio_string_block.bio_string = decoded;
			 if(cnt==0)
			 {
				std::cout<<block.qname<<std::endl;
				std::cout<<decoded<<std::endl;
				std::cout<<block.real_bio_string<<std::endl;
			 }
			 cnt++;
			 memcpy(bio_string_block.qname, block.qname, block.lqname);
			 bio_string_block.lqname = block.lqname;
			 bio_string_block.window_id = block.window_id;
			 bio_string_block.window_id2 = block.window_id2;
			 bio_string_block.hap_id = block.hap_id;
			 bio_string_block.hap_id2 = block.hap_id2;
			 bio_string_block.hap_offset = block.hap_offset;
			 bio_string_block.hap_offset2 = block.hap_offset2;
			 bio_string_block.hapid1 = block.hapid1;
			 bio_string_block.hapid2 = block.hapid2;
			 bio_string_blocks.push_back(bio_string_block);
		 }
		 bio_string_to_num(bio_string_blocks, hap_id_list, pre_hap_id);
	 }
 }
void go_to_decompress(uint32_t hap_id_list[WINDOW_ID_MAX])
{
     std::string filename;
	 uint32_t pre_hap_id = 0;
	 for(int i = 0;i<=0;i++)
	 {
		 double start_time = clock();
		 filename = "/home/user/yexiang/final_store_bgi/store." + std::to_string(i) +".bin";
		 read_from_bin(filename, hap_id_list, &pre_hap_id);
		 std::cout<<"successfully decompress "<< i <<".bin"<<std::endl;
		 double end_time = clock();
		 std::cout<<"spend "<<end_time - start_time<<std::endl;
	 }
}
void DECOMPRESS_MAIN::run(int argc, char* argv[])
{
	std::cout<<"成功到这里"<<std::endl;
	 double cpu_time = cputime();
     CLASSIFY_SHARE_DATA *share = NULL;
	 //share data
	 share = (CLASSIFY_SHARE_DATA*)xcalloc(1, sizeof(CLASSIFY_SHARE_DATA));
 
	 //loading parameters
	 share->o = (OL_PAR *)xcalloc(1, sizeof(OL_PAR));//mapping option
	 if(share->o->get_option(argc, argv) != 0)
		 return;
	 const char* statics_file_names[share->o->thread_n];
	 FILE* statics_file_handler[share->o->thread_n];
 
	 //open SAM/BAM files
	 char * input_bam_fn = share->o->read_bam;
	 htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	 bam_hdr_t *header = sam_hdr_read(input_file);
	 share->header = header;
 
	 char *work_dir = share->o->work_dir;

	 std::string out_file_name = std::string(work_dir) + "/fq.bin";
	 share->sam_out_fp = fopen(out_file_name.c_str(), "w");
	 fprintf(stderr, "BEGIN LOADING INDEX\n");
	 share->idx = (Aln_online::IDX_loader *)xcalloc(1, sizeof(Aln_online::IDX_loader));
	 fprintf(stderr, "Restore ALL INDEX\n");
	 share->idx->restore_idx(share->o);
	 fprintf(stderr, "load index CPU: %.3f sec\n", cputime() - cpu_time);
     read_hap_id_list();
     go_to_decompress(hap_id_list);
     return;
}
