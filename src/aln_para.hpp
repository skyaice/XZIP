/*
 * Online_para.hpp
 *
 *  Created on: 2022年12月2日
 *      Author: fenghe
 */

 #ifndef SRC_BWT_ONLINE_ALN_PARA_HPP_
 #define SRC_BWT_ONLINE_ALN_PARA_HPP_
 #include "CPPLIB/get_option_cpp.hpp"
 #include <string>
 #include <cstring>
 #include <cstdlib>
 #include <iostream>
 
 static std::string remove_trailing_slash(std::string path) {
	 while (path.size() > 1 && (path.back() == '/' || path.back() == '\\')) {
		 path.pop_back();
	 }
	 return path;
 }
 
 static std::string ensure_trailing_slash(std::string path) {
	 if (!path.empty() && path.back() != '/' && path.back() != '\\') {
		 path += "/";
	 }
	 return path;
 }
 
 static std::string join_path(const std::string& parent, const std::string& child) {
	 if (parent.empty()) return child;
 
	 if (parent.back() == '/' || parent.back() == '\\') {
		 return parent + child;
	 }
 
	 return parent + "/" + child;
 }
 bool create_directory(const std::string& dir_path);
 struct OL_PAR{
	 //basic option
	 int thread_n;//the max thread number of execute sub thread
	 //
	 int min_base_qual;
	 int read_length;
	 int max_use_read;
 
	 int webp_lossless;
	 int webp_quality;
	 int webp_method;
	 int webp_width;
	 int webp_height;
	 int zstd_level;
	 int paired_end;
 
	 char * idx_dir;
	 char * read_bam;
	 //char * read_fastq2;
	 char * work_dir;
	 char * pos_dir;
	 char * quality_score_dir;
	 char * webp_dir;
	 char * fastq_output_dir;
	 char * fasta_path;
	 char * unmapped_read1_file;
	 char * unmapped_read2_file;
 
 
	 char *sample_name;
	 char *readset;
	 char *pl;
	 char *cn;
	 char *rg;
 
	 int show_assembly = 0;
	 int exact_statics = 0;
	 int lowseq = 0;
 
	 int seeding_level;
 
	 
	 int get_option(int argc, char *argv[]) {
		 options_list l;
	 
		 l.show_command(stderr, argc + 1, argv - 1);
	 
		 l.add_title_string("\n");
		 l.add_title_string("  Usage:     ");
		 l.add_title_string(PACKAGE_NAME);
		 l.add_title_string(" bwt_aln [Options] <idx_dir> <read_bam> <work_dir> <fasta_path> [rg]\n");
		 l.add_title_string("             or use --idx-dir --read-bam --work-dir --fasta-path [--rg]\n");
	 
		 l.add_title_string("  Basic:\n");
		 l.add_title_string("    <idx_dir>      Folder  Directory used for storing index.\n");
		 l.add_title_string("    <read_bam>     File    Input BAM file.\n");
		 l.add_title_string("    <work_dir>     Folder  Working directory. Subfolders will be created automatically.\n");
		 l.add_title_string("    <fasta_path>   File    Reference FASTA file.\n");
		 l.add_title_string("    [rg]           String  Optional read group.\n");
		 l.add_title_string("    Long path options accept both '-' and '_' spellings, e.g. --idx-dir/--idx_dir.\n");
		 l.add_title_string("    Read length accepts --read_length, --read-length, --rlen, or -n.\n");
		 l.add_title_string("\n");
	 
		 l.add_title_string("  Quality-score compression options:\n");
		 l.add_title_string("    -W / --webp_lossless / --quality_lossless\n");
		 l.add_title_string("                         1=zstd lossless quality scores, 0=WebP lossy. Default: 1\n");
		 l.add_title_string("    -Q / --webp_quality     WebP lossy quality, 0-100. Used only when -W 0. Default: 80\n");
		 l.add_title_string("    -M / --webp_method      WebP method, 0-6. Used only when -W 0. Default: 6\n");
		 l.add_title_string("    -x / --webp_width       WebP tile width, 1-16383. Rounded down to read_length multiple. Default: 4096\n");
		 l.add_title_string("    -y / --webp_height      WebP tile height, 1-16383. Default: 4096\n");
		 l.add_title_string("    -Z / --zstd_level       zstd compression level, 1-19. Default: 19\n");
		 l.add_title_string("    --paired_end            1=paired-end BAM, 0=single-end BAM. Default: 1\n");
		 l.add_title_string("\n");
	 
		 /*
			 Basic options
		 */
		 l.add_option("read_length", 'n', "read length", true, 100);
		 l.add_alias_back("rlen");
		 l.set_arg_pointer_back((void *)&read_length);
	 
		 l.add_option("max_use_read", 'r', "max_use_read", true, MAX_int32t);
		 l.set_arg_pointer_back((void *)&max_use_read);
	 
		 l.add_option("seeding_level", 'l', "seeding_level", true, 1);
		 l.set_arg_pointer_back((void *)&seeding_level);
	 
		 l.add_option("thread", 't', "Number of threads", true, 4);
		 l.set_arg_pointer_back((void *)&thread_n);
 
		 l.add_option("paired_end", 'e', "[0/1] Input BAM layout: 1=paired-end, 0=single-end", true, 1);
		 l.set_arg_pointer_back((void *)&paired_end);
 
		 const char *opt_idx_dir = NULL;
		 const char *opt_read_bam = NULL;
		 const char *opt_work_dir = NULL;
		 const char *opt_fasta_path = NULL;
		 const char *opt_rg = NULL;
 
		 l.add_option("idx_dir", 'i', "Index directory", false, "");
		 l.set_arg_pointer_back((void *)&opt_idx_dir);
 
		 l.add_option("read_bam", 'b', "Input BAM file", false, "");
		 l.set_arg_pointer_back((void *)&opt_read_bam);
 
		 l.add_option("work_dir", 'd', "Working directory", false, "");
		 l.set_arg_pointer_back((void *)&opt_work_dir);
 
		 l.add_option("fasta_path", 'f', "Reference FASTA file", false, "");
		 l.set_arg_pointer_back((void *)&opt_fasta_path);
 
		 l.add_option("rg", 'G', "Optional read group", false, "");
		 l.set_arg_pointer_back((void *)&opt_rg);
	 
		 l.add_option("ASS_detail", 'A', "[Bool] Whether show the detail of assembly in VC.");
		 l.set_arg_pointer_back((void *)&show_assembly);
	 
		 l.add_option("Exact_Gap_statistic", 'g', "[Bool] statistic exact aln num and gap aln num");
		 l.set_arg_pointer_back((void *)&exact_statics);
	 
		 l.add_option("sample_name", 's', "sample name", false, "");
		 l.set_arg_pointer_back((void *)&sample_name);
	 
		 l.add_option("plant_lane", 'p', "plant lane", false, "");
		 l.set_arg_pointer_back((void *)&pl);
	 
		 l.add_option("company_name", 'c', "company name", false, "");
		 l.set_arg_pointer_back((void *)&cn);
	 
		 l.add_option("readset", 'R', "readset", false, "");
		 l.set_arg_pointer_back((void *)&readset);
	 
		 /*
			 WebP compression options
	 
			 -L 1 : lossless
			 -L 0 : lossy
			 -Q   : lossy quality, 0-100
			 -M   : WebP method, 0-6
		 */
		 webp_lossless = 1;
		 webp_quality = 80;
		 webp_method = 6;
		 webp_width = 4096;
		 webp_height = 4096;
		 zstd_level = 19;
	 
		 l.add_option("webp_lossless", 'W', "[0/1] Quality mode: 1=zstd lossless quality scores, 0=WebP lossy", true, 1);
		 l.add_alias_back("quality_lossless");
		 l.set_arg_pointer_back((void *)&webp_lossless);
	 
		 l.add_option("webp_quality", 'Q', "WebP lossy quality, range 0-100, only used when webp_lossless=0", true, 80);
		 l.set_arg_pointer_back((void *)&webp_quality);
	 
		 l.add_option("webp_method", 'M', "WebP method, range 0-6. Higher means slower but better compression; only used when webp_lossless=0", true, 6);
		 l.set_arg_pointer_back((void *)&webp_method);

		 l.add_option("webp_width", 'x', "WebP tile width, range 1-16383. Rounded down to a read_length multiple", true, 4096);
		 l.add_alias_back("webp-width");
		 l.set_arg_pointer_back((void *)&webp_width);

		 l.add_option("webp_height", 'y', "WebP tile height, range 1-16383", true, 4096);
		 l.add_alias_back("webp-height");
		 l.set_arg_pointer_back((void *)&webp_height);

		 l.add_option("zstd_level", 'Z', "zstd compression level, range 1-19. Higher means slower but better compression", true, 19);
		 l.add_alias_back("zstd-level");
		 l.set_arg_pointer_back((void *)&zstd_level);
	 
		 /*
			 Parse options
		 */
		 if (l.default_option_handler(argc, argv)) {
			 l.show_c_value(stderr);
			 return 1;
		 }
	 
		 l.show_c_value(stderr);
	 
		 /*
			 New positional arguments:
	 
			 argv[optind + 0] = idx_dir
			 argv[optind + 1] = read_bam
			 argv[optind + 2] = work_dir
			 argv[optind + 3] = fasta_path
			 argv[optind + 4] = rg, optional
		 */
		 int positional_arg = optind;
		 const char *idx_dir_arg = opt_idx_dir;
		 const char *read_bam_arg = opt_read_bam;
		 const char *work_dir_arg = opt_work_dir;
		 const char *fasta_path_arg = opt_fasta_path;
		 const char *rg_arg = opt_rg;
 
		 if (idx_dir_arg == NULL && positional_arg < argc) {
			 idx_dir_arg = argv[positional_arg++];
		 }
		 if (read_bam_arg == NULL && positional_arg < argc) {
			 read_bam_arg = argv[positional_arg++];
		 }
		 if (work_dir_arg == NULL && positional_arg < argc) {
			 work_dir_arg = argv[positional_arg++];
		 }
		 if (fasta_path_arg == NULL && positional_arg < argc) {
			 fasta_path_arg = argv[positional_arg++];
		 }
		 if (rg_arg == NULL && positional_arg < argc) {
			 rg_arg = argv[positional_arg++];
		 }
 
		 if (idx_dir_arg == NULL || read_bam_arg == NULL ||
			 work_dir_arg == NULL || fasta_path_arg == NULL) {
			 return l.output_usage();
		 }
	 
		 /*
			 Check normal options
		 */
		 xassert((thread_n >= 1) && (thread_n <= 256),
				 "Input error: thread_n cannot be less than 1 or more than 256\n");
 
		 xassert((paired_end == 0 || paired_end == 1),
				 "Input error: paired_end must be 0 or 1\n");
	 
		 if (min_base_qual == 0) {
			 min_base_qual = 15;
		 }
	 
		 /*
			 Check WebP options
		 */
		 xassert((webp_lossless == 0 || webp_lossless == 1),
				 "Input error: webp_lossless must be 0 or 1\n");
	 
		 xassert((webp_quality >= 0 && webp_quality <= 100),
				 "Input error: webp_quality must be between 0 and 100\n");
	 
		 xassert((webp_method >= 0 && webp_method <= 6),
				 "Input error: webp_method must be between 0 and 6\n");

		 xassert((webp_width >= 1 && webp_width <= 16383),
				 "Input error: webp_width must be between 1 and 16383\n");

		 xassert((webp_height >= 1 && webp_height <= 16383),
				 "Input error: webp_height must be between 1 and 16383\n");

		 xassert((webp_width >= read_length),
				 "Input error: webp_width must be at least read_length\n");

		 xassert((zstd_level >= 1 && zstd_level <= 19),
				 "Input error: zstd_level must be between 1 and 19\n");
	 
		 /*
			 Read positional arguments
		 */
		 std::string idx_dir_str = ensure_trailing_slash(idx_dir_arg);
		 std::string read_bam_str = read_bam_arg;
		 std::string work_dir_str = remove_trailing_slash(work_dir_arg);
		 std::string fasta_path_str = fasta_path_arg;
	 
		 /*
			 Auto-generated subdirectories
		 */
		 std::string tmp_dir_str           = join_path(work_dir_str, "tmp");
		 std::string pos_dir_str           = join_path(work_dir_str, "final_store_with_name");
		 std::string quality_score_dir_str = join_path(work_dir_str, "quality_score");
		 std::string webp_dir_str          = join_path(work_dir_str, "store_2_zero_lossless");
	 
		 /*
			 Other internal directories.
			 Keep these if later code writes into them.
		 */
		 std::string diff_base_dir_str     = join_path(work_dir_str, "diff_base");
		 std::string quality_diff_dir_str  = join_path(work_dir_str, "quality_diff_base");
		 std::string byte_flags_dir_str    = join_path(work_dir_str, "byte_flags");
	 
		 /*
			 Create directories.
			 Your create_directory() only creates one level, so create work_dir first.
		 */
		 if (!create_directory(work_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(tmp_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(pos_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(quality_score_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(webp_dir_str)) {
			 return 1;
		 }
	 
		 /*
			 Optional internal directories.
			 If one of these fails, I suggest returning 1, because later code may write into them.
		 */
		 if (!create_directory(diff_base_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(quality_diff_dir_str)) {
			 return 1;
		 }
	 
		 if (!create_directory(byte_flags_dir_str)) {
			 return 1;
		 }
	 
		 /*
			 Assign back to original char* members
		 */
		 idx_dir = strdup(idx_dir_str.c_str());
		 read_bam = strdup(read_bam_str.c_str());
		 work_dir = strdup(work_dir_str.c_str());
	 
		 pos_dir = strdup(pos_dir_str.c_str());
		 quality_score_dir = strdup(quality_score_dir_str.c_str());
		 webp_dir = strdup(webp_dir_str.c_str());
	 
		 fasta_path = strdup(fasta_path_str.c_str());
	 
		 /*
			 Optional rg
		 */
		 if (rg_arg != NULL && rg_arg[0] != '\0') {
			 rg = strdup(rg_arg);
		 } else {
			 rg = NULL;
		 }
	 
		 /*
			 Print compression mode for sanity check
		 */
		 fprintf(stderr,
				 "[Quality] mode=%s, webp_quality=%d, webp_method=%d, webp_width=%d, webp_height=%d, zstd_level=%d\n",
				 webp_lossless ? "zstd-lossless" : "webp-lossy",
				 webp_quality,
				 webp_method,
				 webp_width,
				 webp_height,
				 zstd_level);
		 fprintf(stderr, "[Input] paired_end=%d (%s)\n",
				 paired_end,
				 paired_end ? "paired-end" : "single-end");
	 
		 return 0;
	 }
 
	 int get_decompress_option(int argc, char *argv[])
	 {
		 options_list l;
		 l.add_title_string("\n");
		 l.add_title_string("  Usage:     ");
		 l.add_title_string(PACKAGE_NAME);
		 l.add_title_string(" decompress [Options] <idx_dir> <pos_dir> <webp_dir> <fasta_path>\n");
		 l.add_title_string("             or use --idx-dir --pos-dir --webp-dir --fasta-path\n");
		 l.add_title_string("    -W / --webp_lossless / --quality_lossless\n");
		 l.add_title_string("                         1=zstd lossless quality scores, 0=WebP lossy. Must match compression. Default: 1\n");
		 l.add_title_string("\n");
		 l.add_option("thread", 't', "Number of threads", true, 4); l.set_arg_pointer_back((void *)&thread_n);
		 l.add_option("read_length", 'n', "read length", true, 150); 
		 l.add_alias_back("rlen");
		 l.set_arg_pointer_back((void *)&read_length);
		 l.add_option("paired_end", 'e', "[0/1] Input layout used during compression: 1=paired-end, 0=single-end", true, 1);
		 l.set_arg_pointer_back((void *)&paired_end);
		 webp_lossless = 1;
		 l.add_option("webp_lossless", 'W', "[0/1] Quality mode used during compression: 1=zstd lossless, 0=WebP lossy", true, 1);
		 l.add_alias_back("quality_lossless");
		 l.set_arg_pointer_back((void *)&webp_lossless);
 
		 const char *opt_idx_dir = NULL;
		 const char *opt_pos_dir = NULL;
		 const char *opt_webp_dir = NULL;
		 const char *opt_fasta_path = NULL;
 
		 l.add_option("idx_dir", 'i', "Index directory", false, "");
		 l.set_arg_pointer_back((void *)&opt_idx_dir);
 
		 l.add_option("pos_dir", 'p', "Position directory", false, "");
		 l.set_arg_pointer_back((void *)&opt_pos_dir);
 
		 l.add_option("webp_dir", 'w', "WebP directory; in zstd-lossless mode it is only used to locate the work directory", false, "");
		 l.set_arg_pointer_back((void *)&opt_webp_dir);
 
		 l.add_option("fasta_path", 'f', "Reference FASTA file", false, "");
		 l.set_arg_pointer_back((void *)&opt_fasta_path);
 
		 if(l.default_option_handler(argc, argv)) {
			 l.show_c_value(stderr);
			 return 1;
		 }
		 l.show_c_value(stderr);
 
		 int positional_arg = optind;
		 const char *idx_dir_arg = opt_idx_dir;
		 const char *pos_dir_arg = opt_pos_dir;
		 const char *webp_dir_arg = opt_webp_dir;
		 const char *fasta_path_arg = opt_fasta_path;
 
		 if (idx_dir_arg == NULL && positional_arg < argc) {
			 idx_dir_arg = argv[positional_arg++];
		 }
		 if (pos_dir_arg == NULL && positional_arg < argc) {
			 pos_dir_arg = argv[positional_arg++];
		 }
		 if (webp_dir_arg == NULL && positional_arg < argc) {
			 webp_dir_arg = argv[positional_arg++];
		 }
		 if (fasta_path_arg == NULL && positional_arg < argc) {
			 fasta_path_arg = argv[positional_arg++];
		 }
 
		 if (idx_dir_arg == NULL || pos_dir_arg == NULL ||
			 webp_dir_arg == NULL || fasta_path_arg == NULL) {
			 return l.output_usage();
		 }
 
		 xassert((thread_n >= 0) && (thread_n <= 256), "Input error: thread_n cannot be less than 1 or more than 256\n");
		 xassert((paired_end == 0 || paired_end == 1), "Input error: paired_end must be 0 or 1\n");
		 xassert((webp_lossless == 0 || webp_lossless == 1), "Input error: webp_lossless must be 0 or 1\n");
 
 
		 idx_dir = strdup(idx_dir_arg);
		 if (idx_dir[strlen(idx_dir) - 1] != '/')
			 strcat(idx_dir, "/");
		 pos_dir = strdup(pos_dir_arg);
		 //read_fastq2 = strdup(argv[optind + 2]);
		 webp_dir = strdup(webp_dir_arg);
		 fasta_path = strdup(fasta_path_arg);
		 fprintf(stderr,
				 "[Quality] mode=%s\n",
				 webp_lossless ? "zstd-lossless" : "webp-lossy");
		 return 0;
	 }
 };
 
 
 
 #endif /* SRC_BWT_ONLINE_ALN_PARA_HPP_ */
 
