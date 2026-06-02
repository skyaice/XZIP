
/*
 * var_call.cpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

extern "C"
{
#include "clib/ksort.h"
}

#include "var_calling/known_var_candidate_generator.hpp"
#include "var_calling/novel_var_candidate_generator.hpp"
#include "var_calling/var_caller_all_type.hpp"
#include "aln_RST_DEF.hpp"
#include "process3_var_call.hpp"

#include "BWT_aln.hpp"

namespace Aln_online{

static void show_simple_vcf_header(FILE* output){
	fprintf(output,
			"##fileformat=VCFv4.2\n"
			"##contig=<ID=chr1,length=248956422>\n"
			"##contig=<ID=chr2,length=242193529>\n"
			"##contig=<ID=chr3,length=198295559>\n"
			"##contig=<ID=chr4,length=190214555>\n"
			"##contig=<ID=chr5,length=181538259>\n"
			"##contig=<ID=chr6,length=170805979>\n"
			"##contig=<ID=chr7,length=159345973>\n"
			"##contig=<ID=chr8,length=145138636>\n"
			"##contig=<ID=chr9,length=138394717>\n"
			"##contig=<ID=chr10,length=133797422>\n"
			"##contig=<ID=chr11,length=135086622>\n"
			"##contig=<ID=chr12,length=133275309>\n"
			"##contig=<ID=chr13,length=114364328>\n"
			"##contig=<ID=chr14,length=107043718>\n"
			"##contig=<ID=chr15,length=101991189>\n"
			"##contig=<ID=chr16,length=90338345>\n"
			"##contig=<ID=chr17,length=83257441>\n"
			"##contig=<ID=chr18,length=80373285>\n"
			"##contig=<ID=chr19,length=58617616>\n"
			"##contig=<ID=chr20,length=64444167>\n"
			"##contig=<ID=chr21,length=46709983>\n"
			"##contig=<ID=chr22,length=50818468>\n"
			"##contig=<ID=chrX,length=156040895>\n"
			"##contig=<ID=chrY,length=57227415>\n"
			"##contig=<ID=chrM,length=16569>\n"
			"##FILTER=<ID=LOW_COV,Description=\"Low coverage\">\n"
			"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant TYPE\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus Genotype across all datasets with called genotype\">\n"
			"##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reads that support REF\">\n"
			"##FORMAT=<ID=DA,Number=1,Type=Integer,Description=\"Number of reads that support ALT\">\n"
			"##FORMAT=<ID=NK,Number=1,Type=String,Description=\"(Novel or Known?)How the variants is called.\">\n"
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INTEGRATION\n"
			"");
}
}
}


