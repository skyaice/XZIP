/*
 * bwt.cpp
 *
 *  Created on: 2022年5月26日
 *      Author: fenghe
 */


#include "bwt.hpp"


char bwa_rg_id[256];
char *bwa_pg;

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[8192];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		if (scanres != 3) goto badread;
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			if (scanres != 2) goto badread;
			p->name = strdup(str);
			// read fasta comments
			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);
			if (c == EOF) {
				scanres = EOF;
				goto badread;
			}
			*q = 0;
			if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			if (scanres != 3) goto badread;
			p->offset = xx;
		}
		err_fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		if (scanres != 3) goto badread;
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			if (scanres != 3) goto badread;
			p->offset = xx;
			p->amb = str[0];
		}
		err_fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;

 badread:
	if (EOF == scanres) {
		//err_fatal(__func__, "Error reading %s : %s\n", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %s\n", fname);
}
bntseq_t *bns_restore(const char *prefix)
{
	char ann_filename[1024], amb_filename[1024], pac_filename[1024];
	bntseq_t *bns;
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	bns = bns_restore_core(ann_filename, amb_filename, pac_filename);
	if (bns == 0) return 0;
	return bns;
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) err_fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}



/***********************
 * SAM header routines *
 ***********************/

void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line, bam_hdr_t *header,  char **sn)
{
	int i, n_SQ=0, j = 0;
	extern char *bwa_pg;
	if (hdr_line) {
		const char *smPtr = strstr(hdr_line, "SM:");
		if (smPtr != NULL) {
			const char *tabPtr = strchr(smPtr, '\t');
			size_t length = tabPtr - smPtr-3;
			*sn = (char *)xcalloc(length + 1, 1);
			strncpy(*sn, smPtr+3, length);
			*sn[length] = '\0';
			
		}
	}

	header->n_targets = 24;
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, sizeof(uint32_t));
	if (hdr_line) {
		const char *p = hdr_line;
		while ((p = strstr(p, "@SQ\t")) != 0) {
			if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
			p += 4;
		}
	}
	if (n_SQ == 0) {
		for (i = bns->n_seqs-24; i < bns->n_seqs; ++i) {
			//err_printf("@SQ\tSN:%s\tLN:%d", bns->anns[i].name, bns->anns[i].len);
			header->target_name[j] = strdup(bns->anns[i].name);
			header->target_len[j] =  bns->anns[i].len;
			++j;

			if (bns->anns[i].is_alt) err_printf("\tAH:*\n");
			else err_fputc('\n', stdout);
		}
	} else if (n_SQ != bns->n_seqs)
		fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
	if (hdr_line) {
		//err_printf("%s\n", hdr_line);
		header->l_text = strlen(hdr_line)+1;
		header->text = (char *)xcalloc(header->l_text,1);
		memcpy(header->text, hdr_line, strlen(hdr_line));
		//memcpy(header->text+strlen(hdr_line), "\n", 1);
		header->text[strlen(hdr_line)]='\n';

	} 
	if (bwa_pg) err_printf("%s\n", bwa_pg);
}

static char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}


//char *bwa_set_rg(const char *sample_name, const char *readset, const char *pl,const char *cn)
char *bwa_set_rg(const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(bwa_rg_id, 0, 256);
	if (strstr(s, "@RG") != s) {
		fprintf(stderr, "[E::%s] the read group line is not started with @RG\n", __func__);
		goto err_set_rg;
	}

	// if (strstr(s, "\t") != NULL) {
	// 	fprintf(stderr, "[E::%s] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n", __func__);
	// 	goto err_set_rg;
	// }
	rg_line = strdup(s);
	bwa_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		fprintf(stderr, "[E::%s] no ID within the read group line\n", __func__);
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		fprintf(stderr, "[E::%s] @RG:ID is longer than 255 characters\n", __func__);
		goto err_set_rg;
	}
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return rg_line;

err_set_rg:
	free(rg_line);
	return 0;
}


char *bwa_insert_header(const char *s, char *hdr)
{
	int len = 0;
	if (s == 0 || s[0] != '@') return hdr;
	if (hdr) {
		len = strlen(hdr);
		//hdr = realloc(hdr, len + strlen(s) + 2);
		hdr[len++] = '\n';
		strcpy(hdr + len, s);
	} else hdr = strdup(s);
	bwa_escape(hdr + len);
	return hdr;
}



