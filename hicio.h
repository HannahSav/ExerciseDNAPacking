#pragma once

typedef struct {
	int i; // Index of the first particle in the pair
	int j; // Index of the second particle
	int n; // Number of times particles i and j appear crosslinked in experiment
} HiCPair;
 

typedef struct {
	int atomCount; // Number of DNA fragments
	int pairCount; // Number of pairs that were crosslinked at least once
	HiCPair* pairs; // Array containing the list of all crosslinked pairs
} HiCData;

void readHiC(const char* filename, HiCData* hic); // Read HiC data from the file
void writeHiC(const char* filename, HiCData* hic);
