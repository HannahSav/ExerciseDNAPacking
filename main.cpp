#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hicio.h"
#include "counter.h"

int main(int argc, char* argv[]){

	if(argc < 2){
		printf("Usage %s <filename>\n", argv[0]);
	}
	printf("Reading %s\n", argv[1]);
	
	// Reading the chromosome data using hicio utility
	HiCData hic;
	readHiC(argv[1], &hic);

	// for (int i = 0; i < hic.pairCount; i++)
	// {
	// 	printf("%d-%d: %d\n", hic.pairs[i].i, hic.pairs[i].j, hic.pairs[i].n);
	// }

	counting(hic);
}