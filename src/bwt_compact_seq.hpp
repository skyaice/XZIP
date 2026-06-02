#include <stdint.h>
#include <string.h>
typedef struct {
	uint64_t bwt_k;
    uint8_t length;
	char* left_seq;
    uint8_t left_length;
}bwt_compact_seq;