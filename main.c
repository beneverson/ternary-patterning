#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include </Users/bernardeverson/Desktop/Research/Koder_Lab/4 Helix Bundle/TernaryPatterning/residue_scores.h>

#define kSequenceInputFilePath "/Users/bernardeverson/Desktop/Research/Koder_Lab/4 Helix Bundle/TernaryPatterning/sequence_input.data"
#define kSequenceOutputFilePath "/Users/bernardeverson/Desktop/Research/Koder_Lab/4 Helix Bundle/TernaryPatterning/sequence_output.data"
#define kNeutralPositionCode 'X'
#define kSurfacePositionCode 'U'
#define kBuriedPositionCode 'Z'
#define kNumberOfProteins 5

// return the surface residue corresponding to the score we input
//
char surfaceResidueForOrder (int order) {
	
	char res; 
	
	switch (order) {
		case 0:
			res = 'E';
			break;
		case 1:
			res = 'A';
			break;
		case 2:
			res = 'Q';
			break;
		case 3:
			res = 'K';
			break;
		case 4:
			res = 'R';
			break;
	}
	
	return res;
}


/*
 
 input: order	an int representing the order, from highest(=0)
				to lowest, of the score to return. 
 
 output: the score of the corresponding surface residue
 
 */

float surfaceScoreForOrder (int order) {
	
	float score; 
	
	switch (order) {
		case 0:
			score = (float)kSurfGluScore;
			break;
		case 1:
			score = (float)kSurfAlaScore;
			break;
		case 2:
			score = (float)kSurfGlnScore;
			break;
		case 3:
			score = (float)kSurfLysScore;
			break;
		case 4:
			score = (float)kSurfArgScore;
			break;
	}
	
	return score;
}

/*
	The inverse cumulative distribution function (CDF) of the 
	distribution of surface residues. 
	
 input: input_num -- the randomly generated float (between 0 and 1)
 
 output: residue character
 
 */
char surfaceInverseCDF (float input_num) {
	
	float integral = (float)0; 
	
	char res = 'X';
	
	// go through each residue, in order
	for (int order = 0; order < 5; order++) {
		
		// add the score of this residue to integral
		integral += surfaceScoreForOrder(order);
	
		// when the incremented number is higher than the 
		// random input number
		if (integral >= input_num) {
			
			// return the last residue we added
			// to 'integral'
			return surfaceResidueForOrder(order);
		}
	
	}
	
	return res;
}

/*
	Generate a surface residue according to the scoring fucntion 
	as established in the <residue_scores.h> header file. 
 
	If this function retuns 'X', then something has gone wrong. 
 */
char sampleSurfaceResidues (void) {
	
	float random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	
	while (random_number == 0.0) {
		random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	}
	
	return surfaceInverseCDF(random_number);
}

// return the core residue corresponding to the score we input
//
char coreResidueForOrder (int order) {
	
	char res; 
	
	switch (order) {
		case 0:
			res = 'M';
			break;
		case 1:
			res = 'L';
			break;
		case 2:
			res = 'I';
			break;
		case 3:
			res = 'A';
			break;
		case 4:
			res = 'F';
			break;
		case 5:
			res = 'V';
			break;
	}
	
	return res;
}


/*
 
 input: order	an int representing the order, from highest(=0)
 to lowest, of the score to return. 
 
 output: the score of the corresponding core residue
 
 */

float coreScoreForOrder (int order) {
	
	float score; 
	
	switch (order) {
		case 0:
			score = (float)kCoreMetScore;
			break;
		case 1:
			score = (float)kCoreLeuScore;
			break;
		case 2:
			score = (float)kCoreIleScore;
			break;
		case 3:
			score = (float)kCoreAlaScore;
			break;
		case 4:
			score = (float)kCorePheScore;
			break;
		case 5:
			score = (float)kCoreValScore;
			break;
	}
	
	return score;
}

/*
 The inverse cumulative distribution function (CDF) of the 
 distribution of core residues. 
 
 input: input_num -- the randomly generated float (between 0 and 1)
 
 output: residue character
 
 */
char coreInverseCDF (float input_num) {
	
	float integral = (float)0; 
	
	char res = 'Z';
	
	// go through each residue, in order
	for (int order = 0; order < 6; order++) {
		
		// add the score of this residue to integral
		integral += coreScoreForOrder(order);
		
		// when the incremented number is higher than the 
		// random input number
		if (integral >= input_num) {
			
			// return the last residue we added
			// to 'integral'
			return coreResidueForOrder(order);
		}
		
	}
	
	return res;
}

/*
 Generate a core residue according to the scoring fucntion 
 as established in the <residue_scores.h> header file. 
 
 If this function retuns 'X', then something has gone wrong. 
 */
char sampleCoreResidues (void) {
	
	float random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	
	while (random_number == 0.0) {
		random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	}
	
	return coreInverseCDF(random_number);
}

// return the neutral residue corresponding to the score we input
//
char neutralResidueForOrder (int order) {
	
	char res; 
	
	switch (order) {
		case 0:
			res = 'A';
			break;
		case 1:
			res = 'Q';
			break;
		case 2:
			res = 'K';
			break;
		case 3:
			res = 'E';
			break;
	}
	
	return res;
}


/*
 
 input: order	an int representing the order, from highest(=0)
 to lowest, of the score to return. 
 
 output: the score of the corresponding neutral residue
 
 */

float neutralScoreForOrder (int order) {
	
	float score; 
	
	switch (order) {
		case 0:
			score = (float)kNeuAlaScore;
			break;
		case 1:
			score = (float)kNeuGlnScore;
			break;
		case 2:
			score = (float)kNeuLysScore;
			break;
		case 3:
			score = (float)kNeuGluScore;
			break;
	}
	
	return score;
}

/*
 The inverse cumulative distribution function (CDF) of the 
 distribution of neutral residues. 
 
 input: input_num -- the randomly generated float (between 0 and 1)
 
 output: residue character
 
 */
char neutralInverseCDF (float input_num) {
	
	float integral = (float)0; 
	
	char res = 'U';
	
	// go through each residue, in order
	for (int order = 0; order < 4; order++) {
		
		// add the score of this residue to integral
		integral += neutralScoreForOrder(order);
		
		// when the incremented number is higher than the 
		// random input number
		if (integral >= input_num) {
			
			// return the last residue we added
			// to 'integral'
			return neutralResidueForOrder(order);
		}
		
	}
	
	return res;
}

/*
 Generate a neutral residue according to the scoring function 
 as established in the <residue_scores.h> header file. 
 
 If this function retuns 'X', then something has gone wrong. 
 */
char sampleNeutralResidues (void) {
	
	float random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	
	while (random_number == 0.0) {
		random_number = (float)rand()/((float)RAND_MAX + (float)1); 
	}
	
	return neutralInverseCDF(random_number);
}


char residueForPosition(char position) {
	
	char residue;

		switch (position) {
			case kSurfacePositionCode:
				residue = sampleSurfaceResidues();
				break;
			case kNeutralPositionCode:
				residue = sampleNeutralResidues();
				break;
			case kBuriedPositionCode:
				residue = sampleCoreResidues();
				break;
		}
	
	return residue;
}

int incrementalChargeForResidue(char residue) {
	int charge = 0; 
	
	// positively charged residues
	if (residue == 'R' || residue == 'K') {
		charge = 1;
	}
	else if (residue == 'E') { // negative residues
		charge = -1;
	}
	
	// all others are neutral
	return charge;
	
}

int main (int argc, const char * argv[]) {
	
	/*
		Get a pointer to the files at 'kSequenceInputFilePath' and 'kSequenceOutputFilePath'. 
	 */
	
	FILE * input_data_file;
	FILE * output_data_file;
	char * input_file_path = kSequenceInputFilePath;
	char * output_file_path = kSequenceOutputFilePath;
	
	// seed the random number generator
	struct timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
	
	// show an error if file at path can't be opened
	if ((input_data_file = fopen(input_file_path, "r")) == NULL) {
		printf("\nError opeining file at path: %s\n", input_file_path);
	}
	
	// show an error if file at path can't be opened
	if ((output_data_file = fopen(output_file_path, "a")) == NULL) {
		printf("\nError opeining file at path: %s\n", output_file_path);
	}

	// initialize our charge counter
	int charge = 0;
	
	// initialize a progress counter
	int progress_counter = 0;
	
	// find the length of the input file
	fseek(input_data_file, 0, SEEK_END);
	int input_data_file_size = ftell(input_data_file);
	fseek(input_data_file, 0, SEEK_SET);
	
	// create a string with the size of the input file 
	char final_sequence[input_data_file_size]; 
	char current_char; 
	char associated_residue;
	
	// go through entire process a specified number of times (kNumberOfProteins)
	for (int proteinsCount = 0; proteinsCount < (int)kNumberOfProteins ; proteinsCount++) {
		
		// the overall charge of the sequence should be less than -10
		while (charge > -15.0) {
			
			/*
				For each position in the input file, calculate a corresponding residue
				and print the result to the string final_sequence. 
			 */
			
			// first seek to the beginning of the file
			fseek(input_data_file, 0, SEEK_SET);
			
			// init a counter to help keep track of where we are
			int file_position_counter = 0;
			
			// re-initialize the charge of the protein
			charge = 0;
			
			// go through the entire file
			while (!feof(input_data_file)) {
				
				current_char = fgetc(input_data_file);
				
				if (current_char == kSurfacePositionCode || current_char == kNeutralPositionCode || 
					current_char == kBuriedPositionCode) {
					
					associated_residue = residueForPosition(current_char);
					
					charge += incrementalChargeForResidue(associated_residue);
					
					// add the associated residue to the output string
					final_sequence[file_position_counter] = associated_residue;
				}
				else {
					
					// add the current char to the output string
					final_sequence[file_position_counter] = current_char;
					
				}
				
				file_position_counter++;
			}
			
			progress_counter++;
			
			if (progress_counter % 1000000 == 0) printf("Progress: %d \n", (int)progress_counter/1000000);
			
		}
		
		printf("%s\n", final_sequence);
		fprintf(output_data_file, "charge = %d \n", charge);
		fputs(final_sequence, output_data_file);
		fputs("\n\n", output_data_file);
		
		// reset the charge again
		charge = 0;

	}
	
	fclose(output_data_file);
	fclose(input_data_file);
	
    return 0;
}
