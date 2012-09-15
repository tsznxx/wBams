#include <stdio.h>
#include "sam.h"
#include <string.h>


/*****************************
Modified from calDepth.c
cd samtools-0.1.8/examples
gcc -g -Wall -O2 -I.. wbamToWig2.c -o wBamToWig -lm -lz -L.. -lbam
*****************************/

char chrom[128];
int  curtid=-1;
uint32_t  curbin=-1;
uint32_t  spancnt=0;
uint32_t  binsize = 50;
float     ratio = 1;
float     binsum=0;


typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

// print bin
void print_bin(int end)
{
	if (spancnt++ ==0 || end !=binsize) // first span or last span
		printf("variableStep chrom=%s span=%d\n", chrom ,end);
	if (binsum && end ) printf("%d\t%-5.2f\n",curbin*binsize+1,binsum/end*ratio);
	binsum=0;
}
// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
	if ((int)pos >= tmp->beg && (int)pos < tmp->end)
	{
		//
		if (tid != curtid)
		{
			if (curtid !=-1)
			{
				if (tmp->in->header->target_len[curtid]/binsize == curbin) // the end of the chrom
					print_bin(tmp->in->header->target_len[curtid]%binsize);
				else
					print_bin(binsize);
			}

			curtid=tid;
			strcpy(chrom,tmp->in->header->target_name[tid]);
			spancnt=0;
		}
		if( pos/binsize != curbin )
		{
			print_bin(binsize); //clear the bin at the same time
			curbin = pos/binsize;
		}
		binsum+=n;
	}
	return 0;
}

int main(int argc, char *argv[])
{
	tmpstruct_t tmp;
	if (argc == 1) {
		fprintf(stderr, "wBamToWig version 2.1\n");
		fprintf(stderr, "Usage: wBamToWig   <in.bam>  [binsize=50] [normratio=1 ]\n\n");
		fprintf(stderr, "\tBAM should be sorted.\n");
		fprintf(stderr, "\tNormalization ratio: nratio= librarysize/1000000. [def=1.0].\n");
		return 1;
	}

	if (argc ==3 ) binsize=atoi(argv[2]);
	if (argc ==4 ) ratio=1/atof(argv[3]); // Normalization ratio
	
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
