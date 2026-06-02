#ifndef SAMPLE_HPP_
#define SAMPLE_HPP_

#include <zlib.h>
#include <string>
#include <vector>
#include <fstream>
#include "../CPPLIB//tools.hpp"

using std::string;
using std::vector;
#define MAX_SAMPLE_NUM 500
class Sample {
public:
    string fq1_path;
    string fq2_path;
    string sample_name;

    gzFile *fq1;
    gzFile *fq2;

    Sample( string fq1_path,  string fq2_path, string sample_name){
        this->fq1_path = fq1_path;
        this->fq2_path = fq2_path;
        this->sample_name = sample_name;
    }

    bool check_file(){
        char fn[1024];
		sprintf(fn, "%s/%s", fq1_path.c_str());
		if (access(fn, F_OK) == -1) {
			fprintf(stderr, "%s not exist!\n", fn);
            return false;
		}
        return true;
    }
    


};

class Sample_list{
    std::vector<Sample*> sample_list;
    
    int sample_num;
    char * fq_list_file_path;
    char *map_load_temp = new char[MAX_SAMPLE_NUM];//10M
    char *temp = (char *)malloc(10000); //10k
    vector<std::string> line_tmp;
    int load_fq_list(char *fn){
        std::ifstream map_f(fn);
        while(true){
            map_f.getline(map_load_temp, MAX_SAMPLE_NUM);
            if(map_f.eof() && (map_load_temp[0]) == 0) break;
            split_string(line_tmp, temp, map_load_temp , " ");
            Sample *sample_input = &Sample(line_tmp[0], line_tmp[1], line_tmp[2]);
            if(sample_input->check_file()){
                sample_list.emplace_back(sample_input);
            }

        }
    }

};

#endif