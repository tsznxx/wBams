#include <stdio.h>
#include "sam.h"
#include <string.h>


/*****************************
Modified from calDepth.c
gcc -g -Wall -O2 -I.. wBamToWigBaseResolution.c -o bamToWigbr -lm -lz -L.. -lbam
*****************************/


char      chrom[10];
uint32_t  lastpos = 0;
float     ratio = 1;

typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	int r;
	tmpstruct_t *tmp = (tmpstruct_t*)data;
	if ((int)pos >= tmp->beg && (int)pos < tmp->end)
	{
		r=strcmp(tmp->in->header->target_name[tid],chrom);
		if (r || pos-lastpos !=1) //chrom change or incontinuous regions. 
		{
			printf("fixedStep chrom=%s start=%d step=1\n",tmp->in->header->target_name[tid],pos+1);
			if (r) strcpy(chrom,tmp->in->header->target_name[tid]);
		}
		printf("%-5.1f\n", n*ratio);
		lastpos=pos;
	}
	return 0;
}

int main(int argc, char *argv[])
{
	tmpstruct_t tmp;
	if (argc == 1) {
		fprintf(stderr, "wBamToWig version 1.1\n");
		fprintf(stderr, "Usage: wBamToWigbr   <in.bam>  [normratio]\n\n");
		fprintf(stderr, "\tBAM should be sorted.\n");
		fprintf(stderr, "\tNormalization ratio: nratio= librarysize/1000000. [def=1.0].\n");
		return 1;
	}

	if (argc ==3 ) ratio=1/atof(argv[2]); // Normalization ratio
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.in = samopen(argv[1], "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
	printf("track type=wiggle_0\n");
	sampileup(tmp.in, -1, pileup_func, &tmp);
	samclose(tmp.in);
	return 0;
}
