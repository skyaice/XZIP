/*
 * minimizer.cpp
 *
 *  Created on: 2022年4月20日
 *      Author: fenghe
 */

#include"minimizer.hpp"
	// emulate 128-bit integers and arrays
	unsigned char seq_nt4_table[256] = {
		0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};
	void Minimizer_generater::mm_sketch(const char *seq, int seq_len, std::vector<mm128_t> &p)
	{
		if(seq_len == 0) return;
		const char *str = seq;
		int len = seq_len;
		uint64_t shift1 = 2 * (mini_kmer_len - 1), mask = (1ULL<<2*mini_kmer_len) - 1, kmer[2] = {0,0};
		int i, j, l, buf_pos, min_pos, kmer_span = 0;
		mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
		tiny_queue_t tq;

		xassert(len > 0 && (window_size > 0 && window_size < 256) && (mini_kmer_len > 0 && mini_kmer_len <= 28), ""); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
		memset(buf, 0xff, window_size * 16);
		memset(&tq, 0, sizeof(tiny_queue_t));

		for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
			int c = seq_nt4_table[(uint8_t)str[i]];
			mm128_t info = { UINT64_MAX, UINT64_MAX };
			if (c < 4) { // not an ambiguous base
				int z;
				if (is_hpc) {
					int skip_len = 1;
					if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
						for (skip_len = 2; i + skip_len < len; ++skip_len)
							if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
								break;
						i += skip_len - 1; // put $i at the end of the current homopolymer run
					}
					tq_push(&tq, skip_len);
					kmer_span += skip_len;
					if (tq.count > mini_kmer_len) kmer_span -= tq_shift(&tq);
				} else kmer_span = l + 1 < mini_kmer_len? l + 1 : mini_kmer_len;
				kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
				kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
				if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1]? 0 : 1; // strand
				++l;
				if (l >= mini_kmer_len && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					info.y = (uint64_t)i<<1 | z;
				}
			} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
			buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
			if (l == window_size + mini_kmer_len - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
				for (j = buf_pos + 1; j < window_size; ++j)
					if (min.x == buf[j].x && buf[j].y != min.y){
						p.emplace_back(buf[j]);
					}
				for (j = 0; j < buf_pos; ++j)
					if (min.x == buf[j].x && buf[j].y != min.y) {
						p.emplace_back(buf[j]);
					}
			}
			if (info.x <= min.x) { // a new minimum; then write the old min
				if (l >= window_size + mini_kmer_len && min.x != UINT64_MAX){
					p.emplace_back(min);
				}
				min = info, min_pos = buf_pos;
			} else if (buf_pos == min_pos) { // old min has moved outside the window
				if (l >= window_size + mini_kmer_len - 1 && min.x != UINT64_MAX) {
					p.emplace_back(min);
				}
				for (j = buf_pos + 1, min.x = UINT64_MAX; j < window_size; ++j) // the two loops are necessary when there are identical k-mers
					if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
				for (j = 0; j <= buf_pos; ++j)
					if (min.x >= buf[j].x) min = buf[j], min_pos = j;
				if (l >= window_size + mini_kmer_len - 1 && min.x != UINT64_MAX) { // write identical k-mers
					for (j = buf_pos + 1; j < window_size; ++j) // these two loops make sure the output is sorted
						if (min.x == buf[j].x && min.y != buf[j].y){
							p.emplace_back(buf[j]);
						}
					for (j = 0; j <= buf_pos; ++j)
						if (min.x == buf[j].x && min.y != buf[j].y){
							p.emplace_back(buf[j]);
						}
				}
			}
			if (++buf_pos == window_size) buf_pos = 0;
		}
		if (min.x != UINT64_MAX)
			p.emplace_back(min);
	}


