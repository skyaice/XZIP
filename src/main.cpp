//============================================================================
// Name        : annotationAnaysis.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <map>
//#include <set>
#include <unordered_set>
#include <string>
#include <algorithm>
extern "C"{
	#include "clib/utils.h"
	#include "clib/bam_file.h"
	#include "clib/desc.h"
	#include <zlib.h>
}

#include"BWT_idx/bwt.hpp"
#include"CPPLIB/tools.hpp"
//#include "SV_handler.hpp"
#include "index_building/haplotype_online.hpp"
#include "BWT_aln.hpp"

#include "BWT_idx/bwt_cpt_sa.hpp"
//#include "variant_caller_main.hpp"

#include"CPPLIB/JsonObject/CJsonObject.hpp"
#include"CPPLIB/get_option_cpp.hpp"





int BWT_aln_main(int argc, char *argv[]){
	fprintf(stderr, "\n\n V1.00\n\n");
	//get option
	BWT_aln::BWT_CLASSIFY_MAIN cm;
	cm.init_run(argc, argv);
	return 0;
}

int print_complementary_seq(int argc, char *argv[]){
	char_seq_reverse(strlen(argv[optind]), argv[optind], true);
	fprintf(stderr, "%s\n", argv[optind]);
	char_seq_reverse(strlen(argv[optind]), argv[optind], true);
	char_seq_reverse(strlen(argv[optind]), argv[optind], false);
	fprintf(stderr, "%s\n", argv[optind]);

	return 0;
}

int get_ref_region(int argc, char *argv[]){

	//load reference string;
	faidx_t *fai = fai_load(argv[optind]);//load reference
	int true_region_load_len = 0;
	char *ref_seq = fai_fetch(fai, argv[optind + 1], &true_region_load_len);
	fprintf(stdout, "%s\n", ref_seq);

	fai_destroy(fai);

	return 0;
}

int bwt_sa_compact(int argc, char *argv[]){
	BWT_SA_compactor sc;
	sc.generate_SA_index_from_BWA_SA(argc, argv);
	return 0;
}

//int variant_calling(int argc, char *argv[]){
//        fprintf(stderr, "variant counting start \n");
//        VAR_CALLING::Variant_Caller_main sh;
//        sh.run(argc, argv);
//        return 0;
//}

int decompress(int argc, char *argv[])
{
	BWT_aln::DECOMPRESS_MAIN dc;
	dc.run(argc, argv);
	return 0;
}

int main(int argc, char *argv[])
{
	//little endian check
	union W{ int a; char b; }c;	c.a = 1; xassert(c.b == 1, "System must be little endian");

	fprintf(stderr, "\n\n main version V1.00\n\n");
	COMMAND_HANDLER ch;

	ch.add_function("decompress", "decompress binary file", decompress);
	ch.add_help_msg_back("Usage: 'decompress the file");


	ch.add_function("bwt_aln", "align and variants calling using BWT", BWT_aln_main);
	ch.add_help_msg_back("Used for BGI T20");



	ch.add_help_msg_back("\n············Variants calling functions·················\n\n");

//    ch.add_function("variant_calling", "variants count from wb block", variant_calling);
//    ch.add_help_msg_back("Usage: 'variant_count [wb_block_path]'");

	ch.add_help_msg_back("\n············Tools·················\n\n");
	ch.add_function("reverse", "print complementary seq", print_complementary_seq);

	ch.add_function("ref_region", "print seq in reference region ", get_ref_region);
	ch.add_help_msg_back("Usage: 'ref_region [ref.fa] [chr1: 10000~10100]'");



	return ch.run(argc, argv);
}
