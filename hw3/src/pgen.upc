#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>
#include "packingDNAseq.h"
#include "new_kmer_hash.h"


int main(int argc, char *argv[]){
	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	
	
	char * input_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_name);
	int64_t lines = nKmers/THREADS;
	int64_t skip = lines * MYTHREAD;
	if (MYTHREAD < (nKmers % THREADS)) {
		lines++;
		skip += MYTHREAD; 
	} else {
		skip += (nKmers % THREADS);
	}
	int64_t read_char = lines * LINE_SIZE;
	int64_t skip_char = skip * LINE_SIZE;
	unsigned char* buffer =
		(unsigned char*) malloc((read_char) * sizeof(unsigned char)); // maybe +5
	
	upc_file_t *input;
	input = upc_all_fopen(input_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL);
	upc_all_fseek(input, skip_char*sizeof(unsigned char), UPC_SEEK_SET);
	int64_t cur_chars_read = upc_all_fread_local(input, buffer, sizeof(unsigned char), read_char,
																							 UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
	
	upc_barrier;
	upc_all_fclose(input);
	
	
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	
	
	 init_LookupTable();
	
	shared int64_t* next = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
	shared kmer_t* memory_heap = (shared kmer_t*)upc_all_alloc(nKmers, sizeof(kmer_t)); 
	
	shared int64_t* hash_table = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
	int64_t i;
	upc_forall(i=0; i<nKmers; i++; &hash_table[i])
		hash_table[i] = -1;
	upc_forall(i=0; i<nKmers; i++; &next[i]) 
		next[i] = -1;
	
	upc_barrier; 
	int64_t k = skip; 
	int64_t point = 0;
	
	
	char left_ext, right_ext;
	
	while (point < cur_chars_read) {
	
		left_ext = (char) buffer[point+KMER_LENGTH+1];
		right_ext = (char) buffer[point+KMER_LENGTH+2];
		add_kmer(next, skip, nKmers, hash_table, memory_heap, &buffer[point], left_ext, right_ext);
		
		point += LINE_SIZE;
		skip ++; 
	}	
	upc_barrier;
	constrTime += gettime();


	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	
	
	
	char cur_contig[MAXIMUM_CONTIG_SIZE+1], unpackedKmer[KMER_LENGTH+1]; 
	
	char output_file_name[50];
	sprintf(output_file_name, "pgen%d.out", MYTHREAD);
	FILE *output = fopen(output_file_name, "w");
	
	point = 0;
	for (;  point < lines * LINE_SIZE; point += LINE_SIZE) {
		char left_ext = (char) buffer[point+KMER_LENGTH+1];
		if (left_ext != 'F') continue;

		memcpy(cur_contig, &buffer[point], KMER_LENGTH * sizeof(char));
		int64_t posInContig = KMER_LENGTH;
		char right_ext = (char) buffer[point+KMER_LENGTH+2];
		
		while(right_ext != 'F' && right_ext != 0) {
			cur_contig[posInContig] = right_ext;
			posInContig++;
			right_ext = lookup_kmer_rext(memory_heap, next, hash_table, nKmers,
																	 (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
		}
		
		
		cur_contig[posInContig] = '\0';
		fprintf(output, "%s\n", cur_contig);
	}
	
	fclose(output);
	
	
	upc_barrier;
	
	
	if(MYTHREAD == 0) {
		upc_free(memory_heap);
		upc_free(next);
		upc_free(hash_table);
	}
	
	
	upc_barrier;
	traversalTime += gettime();

	if(MYTHREAD == 0) {
		output = fopen("pgen.out", "w");
		
		for(int i = 0; i < THREADS; ++ i) {
			char str[100];
			sprintf(str, "pgen%d.out", i);
			FILE*input = fopen(str, "r");
			while(fscanf(input, "%s", cur_contig) == 1) {
				fprintf(output, "%s\n", cur_contig); 
			}
			fclose(input);
		}
		fclose(output);
	}


	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
