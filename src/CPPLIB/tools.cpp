/*
 * tools.cpp
 *
 *  Created on: 2021年10月21日
 *      Author: fenghe
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include"tools.hpp"
#include <fstream>

extern "C"{
	#include "../clib/utils.h"
}


void split_string_append(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str){
	if(strlen(split_line) == 0) return;
	strcpy(temp, split_line);
	char * token_value = NULL;
	char *save_ptr = NULL;
	token_value = strtok_r(temp, split_str, &save_ptr); item_value.emplace_back(token_value);
	for(int item_idx = 1; item_idx < MAX_LINE_ITEM_NUM; item_idx++){
		token_value = strtok_r(NULL, split_str, &save_ptr);
		if(token_value == NULL)
			break;
		item_value.emplace_back(token_value);
	}
}

void split_string(std::vector<std::string> &item_value, char * temp, const char * split_line, const char *split_str){
	item_value.clear();
	split_string_append(item_value, temp, split_line, split_str);
}

int get_chr_ID_from_string(std::string & chr_str){
	int chrID = 0;
	if(chr_str.empty() || chr_str.size() > 5 ) //begin with 'chr', like chr1 - chr22 - chrM
		return -1;
	const char *chr_char = chr_str.c_str();
	if(chr_char[0] == 'c')//begin with 'chr', like chr1 - chr22 - chrM
		chr_char += 3;
	int ID_length = strlen(chr_char);
	if(ID_length < 1 || ID_length > 2)
		return -1;
	if(chr_char[0] <= '9' && chr_char[0] >= '1' && (ID_length == 1 || (chr_char[1] <= '9' && chr_char[1] >= '0')))
		chrID = atoi(chr_char) - 1;
	else if(chr_char[0] == 'X' || chr_char[0] == 'x')			chrID = 22;
	else if(chr_char[0] == 'Y' || chr_char[0] == 'y')			chrID = 23;
	else if(chr_char[0] == 'M' || chr_char[0] == 'm')			chrID = 24;
	else 														chrID = -1;
	return chrID;
}

#define ADJUST_INDICES(start, end, len)     \
    if (end > len)                          \
        end = len;                          \
    else if (end < 0) {                     \
        end += len;                         \
        if (end < 0)                        \
        end = 0;                            \
    }                                       \
    if (start < 0) {                        \
        start += len;                       \
        if (start < 0)                      \
        start = 0;                          \
    }

int _string_tailmatch(const std::string&self, const std::string&substr, int start, int end, int direction)
{
	int selflen = (int)self.size();
	int slen = (int)substr.size();
	const char* str = self.c_str();
	const char* sub = substr.c_str();
	ADJUST_INDICES(start, end, selflen);
	if (direction < 0){
		if (start + slen>selflen)				return 0;
	}
	else{
		if (end - start<slen || start>selflen)	return 0;
		if (end - slen > start)					start = end - slen;
	}
	if (end - start >= slen)					return !memcmp(str + start, sub, slen);
	return 0;
}

bool endswith(const std::string&str, const std::string&suffix, int start, int end)
{
	int result = _string_tailmatch(str, suffix, start, end, +1);
	return static_cast<bool>(result);
}

bool startswith(const std::string&str, const std::string&suffix, int start, int end)
{
	int result = _string_tailmatch(str, suffix, start, end, -1);
	return static_cast<bool>(result);
}

#define MAX_INDEL_LEN 50
//return 0 if SNP; return 1 if INDEL; return 2 if SV
int get_VAR_TYPE(std::string & REF, std::string & ALT, int max_index_len){
	if(REF.size() < 1){ fprintf(stderr, "FATAL ERROR, is_SNP???\n"); exit(-1);}
	if(ALT.size() < 1){ fprintf(stderr, "FATAL ERROR, is_SNP???\n"); exit(-1);}

	int indel_len = (REF.size() - ((REF[0] == '-')?1:0)) - (ALT.size() - ((ALT[0] == '-')?1:0));
	if(indel_len < 0) indel_len = - indel_len;
	if(indel_len > MAX_INDEL_LEN) 	return 2;//SV

	if(REF[0] == '-') return 1;//INDEL
	if(ALT[0] == '-') return 1;//INDEL

	if(REF.size() == 1 && ALT.size() == 1) 	return 0;//SNP
	if(REF.size() != ALT.size())			return 1;//INDEL
	int mis_numer = 0;
	int first_diff = -1;
	for(int i = 0; i < (int)REF.size(); i++){
		if(REF[i] != ALT[i]){
			mis_numer++;
			if(first_diff == -1) first_diff = i;
		}
	}
	if(mis_numer <= 1){
		if(first_diff == -1) first_diff = 0;
		char ref_char[2];ref_char[0] = REF[first_diff]; ref_char[1] = 0;
		char alt_char[2];alt_char[0] = ALT[first_diff]; alt_char[1] = 0;
		REF.clear(); REF.append(ref_char);
		ALT.clear(); ALT.append(alt_char);
		return 0;//SNP
	}

	return 1;//INDEL
}

int get_SNP_change_int(std::string & REF, std::string & ALT){
	int ref_int = -1; int alt_int = -1;
	if(REF[0] == 'A') ref_int = 0;
	else if(REF[0] == 'C') ref_int = 1;
	else if(REF[0] == 'G') ref_int = 2;
	else if(REF[0] == 'T') ref_int = 3;
	if(ALT[0] == 'A') alt_int = 0;
	else if(ALT[0] == 'C') alt_int = 1;
	else if(ALT[0] == 'G') alt_int = 2;
	else if(ALT[0] == 'T') alt_int = 3;

	if(ref_int == -1 || alt_int == -1){	fprintf(stderr, "fatal error, get_SNP_change_int Wrong change, and set change int to 0: @ REF:%s, ALT: %s", REF.c_str(), ALT.c_str());	return 0;}
	return ref_int*4 + alt_int;
}

//return true when it is an Ti, return false when it is an Tv
bool isSNP_an_Ti(int var_type, std::string & REF, std::string & ALT){
	if(var_type != 0){		fprintf(stderr, "fatal error, in isSNP_an_Ti");		exit(-1);	}
	if(REF.size() > 1 || ALT.size() > 1){
		fprintf(stderr, "fatal error, in isSNP_an_Ti REF.size() > 1 ");	exit(-1);
	}
	int ref_type = -1; int alt_type = -1;
	if(REF[0] == 'A' || REF[0] == 'G')			ref_type = 0;
	else if(REF[0] == 'C' || REF[0] == 'T')		ref_type = 1;
	else{		fprintf(stderr, "fatal error, in isSNP_an_Ti UNKOWN type REF, treat as false");		return false;}
	if(ALT[0] == 'A' || ALT[0] == 'G')			alt_type = 0;
	else if(ALT[0] == 'C' || ALT[0] == 'T')		alt_type = 1;
	else{		fprintf(stderr, "fatal error, in isSNP_an_Ti UNKOWN type ALT, treat as false");		return false;}

	if(alt_type == ref_type) return true;
	else return false;
}

int getIndel_len(std::string & REF, std::string & ALT){
	return ((ALT.size() - ((ALT[0] == '-')?1:0)) - (REF.size() - ((REF[0] == '-')?1:0)));
}

void load_strings_from_file(char * string_fn, std::vector<std::string> &v, uint64_t max_load_line){
	//get file names
	FILE * try_f = xopen(string_fn, "rb");
	fclose(try_f);
	char *temp = new char[MAX_LINE_LENGTH];//10M
	std::ifstream name_list_File(string_fn);
	uint64_t c_load = 0;
	while(true){
		if(c_load++ > max_load_line) break;
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(name_list_File.eof() && (temp[0]) == 0) break;
		v.emplace_back(temp);
		if(name_list_File.eof())	break;
	}
	name_list_File.close();
}

void load_int_from_file(char * int_fn, std::vector<int> &v, uint64_t max_load_line){
	FILE * try_f = xopen(int_fn, "rb");
	fclose(try_f);
	//get file names
	char *temp = new char[MAX_LINE_LENGTH];//10M
	std::ifstream name_list_File(int_fn);
	uint64_t c_load = 0;
	while(true){
		if(c_load++ > max_load_line) break;
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)	break;
		v.emplace_back(atoi(temp));
	}
	name_list_File.close();
}

void dump_int64_v_to_file(char * out_fn, std::vector<int64_t> &v){
	FILE * try_f = xopen(out_fn, "w");
	for(uint64_t d: v)
		fprintf(try_f, "%ld\n", d);
	fclose(try_f);
}

void load_int64_v_from_file(char * int_fn, std::vector<int64_t> &v){
	FILE * try_f = xopen(int_fn, "rb");
	fclose(try_f);
	//get file names
	char *temp = new char[MAX_LINE_LENGTH];//10M
	std::ifstream name_list_File(int_fn);
	while(true){
		name_list_File.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)	break;
		v.emplace_back(atol(temp));
	}
	name_list_File.close();
}

void CPP_vector_dump_bin(FILE * f_dump, void * data, uint64_t data_size){
	fwrite(&data_size, 1, 8, f_dump);
	fwrite(data, 1, data_size, f_dump);
}

//load from file, return the true load data size(in byte)
uint64_t CPP_vector_load_bin(FILE * f_dump, void ** data, uint64_t type_size){
	uint64_t data_size;
	xread(&data_size, 1, 8, f_dump);//equal to fread ....
	(*data) = (void* ) xmalloc (data_size);
	xread(*data, 1, data_size, f_dump);
	return data_size/type_size;
}

void set_magic_string(FILE * f_dump, const char * magic_string){
	uint32_t magic_len = strlen(magic_string);
	fwrite(magic_string, 1, magic_len, f_dump);
}

void check_magic_string(FILE * f_dump, const char * magic_string){
	uint32_t magic_len = strlen(magic_string);
	char magic_string_load[128];
	xread(magic_string_load, 1, magic_len, f_dump); //equal to fread ....
	//char warning_str[1024];
	//sprintf(warning_str, "Magic string check fail! magic_string %s :  magic_string_load %s", magic_string, magic_string_load);
	xassert(strncmp(magic_string, magic_string_load, magic_len) == 0, "Magic string check fail! "); //equal to assert ....
}

