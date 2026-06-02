/*
 * tools.hpp
 *
 *  Created on: 2021年10月21日
 *      Author: fenghe
 */

#ifndef TOOLS_HPP_
#define TOOLS_HPP_

#define AF_LEVEL_number 9
#define MAX_LINE_LENGTH 10000000
#define MAX_LINE_ITEM_NUM 2000000

//void _err_fatal_simple_core(const char *func, const int line, const char *msg);
#define xassert(cond, msg) 		if ((cond) == 0) _err_fatal_simple_core(__func__,__LINE__, msg)

void split_string_append(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str);
void split_string(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str);
int get_chr_ID_from_string(std::string & chr_str);
bool endswith(const std::string&str, const std::string&suffix, int start, int end);
bool startswith(const std::string&str, const std::string&suffix, int start, int end);

int get_VAR_TYPE(std::string & REF, std::string & ALT, int max_index_len);
void get_SNP_change_str(std::string & REF, std::string & ALT, int AF_level, std::string &SNP_str, std::string &SNP_str_AF);
int get_SNP_change_int(std::string & REF, std::string & ALT);
bool isSNP_an_Ti(int var_type, std::string & REF, std::string & ALT);
int getIndel_len(std::string & REF, std::string & ALT);
void get_INDEL_change_str(std::string & REF, std::string & ALT, int AF_level, std::string &INDEL_str, std::string &INDEL_str_AF);

void load_strings_from_file(char * string_fn, std::vector<std::string> &v, uint64_t max_load_line);
void load_int_from_file(char * int_fn, std::vector<int> &v, uint64_t max_load_line);

void dump_int64_v_to_file(char * out_fn, std::vector<int64_t> &v);
void load_int64_v_from_file(char * int_fn, std::vector<int64_t> &v);

void CPP_vector_dump_bin(FILE * f_dump, void * data, uint64_t data_size);
uint64_t CPP_vector_load_bin(FILE * f_dump, void ** data, uint64_t type_size);
void set_magic_string(FILE * f_dump, const char * magic_string);
void check_magic_string(FILE * f_dump, const char * magic_string);


#endif /* TOOLS_HPP_ */
