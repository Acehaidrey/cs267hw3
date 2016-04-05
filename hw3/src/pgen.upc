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


//#define SINGLE_OUTPUT_FILE


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
	// All the data are stored in  * buffer *
	
	// close the open file
	upc_barrier;
	upc_all_fclose(input);
	
	
	//printf("Reading Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);
	
	
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	
	//int ok=1;
	
	/* Initialize lookup table that will be used for the DNA packing routines */
	 init_LookupTable();
	
	// generate hash table and heap
	shared int64_t* next = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
	shared kmer_t* memory_heap = (shared kmer_t*)upc_all_alloc(nKmers, sizeof(kmer_t)); 
	
	// initialize hash_table
	shared int64_t* hash_table = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
	int64_t i;
	upc_forall(i=0; i<nKmers; i++; &hash_table[i]) // initial hash table
		hash_table[i] = -1;
	upc_forall(i=0; i<nKmers; i++; &next[i]) // initial next index
		next[i] = -1;
	
	upc_barrier; // need to synchronize before we actually start!
	
	int64_t k = skip; // global kmer index
	int64_t point = 0;
	
	
//  i = 0; // test purpose
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

	//printf("Building Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	
	
	//ok=1;
	
	
	
	char cur_contig[MAXIMUM_CONTIG_SIZE+1], unpackedKmer[KMER_LENGTH+1]; 
	
	char output_file_name[50];
	sprintf(output_file_name, "pgen%d.out", MYTHREAD);
	FILE *output = fopen(output_file_name, "w"); // asynchronized & independent output
	
	point = 0;
	for (;  point < lines * LINE_SIZE; point += LINE_SIZE) {
		char left_ext = (char) buffer[point+KMER_LENGTH+1];
		if (left_ext != 'F') continue;

		memcpy(cur_contig, &buffer[point], KMER_LENGTH * sizeof(char));
		int64_t posInContig = KMER_LENGTH;
		char right_ext = (char) buffer[point+KMER_LENGTH+2];
		
		// Keep traversing the graph
		while(right_ext != 'F' && right_ext != 0) {
			cur_contig[posInContig] = right_ext;
			posInContig++;
			right_ext = lookup_kmer_rext(memory_heap, next, hash_table, nKmers,
																	 (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
		}
		
		
		// print the contig
		cur_contig[posInContig] = '\0';
		fprintf(output, "%s\n", cur_contig);
	}
	
	// close the output file
	fclose(output);
	
	
	// clean the allocated memory
	upc_barrier;
	
	
	//printf("Traversal Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);
	
	
	if(MYTHREAD == 0) {
		upc_free(memory_heap);
		upc_free(next);
		upc_free(hash_table);
	}
	
	
	upc_barrier;
	traversalTime += gettime();

// hack implementation to produce a single output file for correctness checking
#ifdef SINGLE_OUTPUT_FILE

	if(MYTHREAD == 0) {
		output = fopen("pgen.out", "w");
		
		for(int t = 0; t < THREADS; ++ t) {
			char str[50];
			sprintf(str, "pgen%d.out", t);
			FILE*input = fopen(str, "r");
			
			while(fscanf(input, "%s", cur_contig) == 1)
				fprintf(output, "%s\n", cur_contig);
			 
			fclose(input);
		}
		
		fclose(output);
	}

#endif


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
