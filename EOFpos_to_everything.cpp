#include <iostream>
#include <assert.h>
#include <sdsl/bit_vectors.hpp>
#include "Parameters.h"

using namespace std;
using namespace sdsl;

dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char

dataTypeNChar buildFreq(string fileName, dataTypeNChar empty_number);
dataTypeNChar remove_empty_symbols(string inputName);
int EOFpos_to_everything(string filename, dataTypeNChar BWT_length, dataTypeNChar empty_number);



int main(int argc, char *argv[]){
	if(argc != 2 ) {
		std::cerr << "usage: " << argv[0] << " input" << std::endl;
		exit(1);
	}
	
	if(strcmp(argv[1], "--help")==0 || strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-h")==0) {
		std::cerr << "usage: " << argv[0] << " input" << std::endl;
		exit(1);
	}
	
	string inputName = argv[1];	
	dataTypeNChar empty_number = remove_empty_symbols(inputName);
	dataTypeNChar BWT_length = buildFreq(inputName,empty_number);
	EOFpos_to_everything(inputName, BWT_length, empty_number);
	return 1;
}



int EOFpos_to_everything(string filename, dataTypeNChar BWT_length, dataTypeNChar empty_number){

	string ext = ".aux";

	tableOcc = new dataTypeNChar*[sizeAlpha];
	//Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		tableOcc[j] = new dataTypeNChar[sizeAlpha];
	}

	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++){
		for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++){
			tableOcc[j][h]=0;
		}
	}
	
	//Open EOFpos file
	dataTypeNSeq nText;
	string eof_s = filename + ".EOFpos";
	FILE* eof = fopen(eof_s.c_str(),"rb");
	if (eof == NULL) {
		std::cerr << "Error opening \"" << eof_s.c_str() << "\" file"<< std::endl;
		exit (1);
	}
	dataTypeNChar numchar = fread (&nText, sizeof(dataTypeNSeq), 1 ,eof);
	if(empty_number>0){
		dataTypeNChar n_bytes=(sizeof(dataTypeNSeq)+sizeof(dataTypeNChar)+sizeof(dataTypedimAlpha))*empty_number;
		fseek(eof,-n_bytes, SEEK_END);
	}
	
	//Open ebwt 
	string bwt_s = filename + ".ebwt";
	FILE* bwt = fopen(bwt_s.c_str(),"r");
	if (bwt == NULL) {
		std::cerr << "Error opening \"" << bwt_s << "\" file"<< std::endl;
		exit (1);		
	}

	FILE* OutFileBWT;
	string InfoFile_s = filename + "_info" + ext;
	FILE* InfoFile = fopen(InfoFile_s.c_str(),"wb");
	if (InfoFile == NULL){
		std::cerr << "Error opening: " << InfoFile_s << std::endl;
		exit (EXIT_FAILURE);
	}
	//Write BWT length, num. sequences and sizeAlpha
	fwrite(&BWT_length, sizeof(dataTypeNChar), 1, InfoFile);
	fwrite(&(freq[alphaInverse[0]]), sizeof(dataTypeNSeq), 1, InfoFile);
	fwrite(&sizeAlpha, sizeof(dataTypedimAlpha), 1, InfoFile);
	
	//Write alphabet symbols
	for(dataTypedimAlpha c = 0; c < sizeAlpha; c++){
		fwrite(&(alphaInverse[c]), sizeof(dataTypedimAlpha), 1, InfoFile);
	}
	
	size_t currentPile=0;

	//currentPile0
	assert(freq[alphaInverse[currentPile]] > 0);
	bit_vector b;
	b = bit_vector(freq[alphaInverse[currentPile]],0);

	//Open output file partial bwt
    string OutFileBWT_s = filename + "_bwt_" + to_string(currentPile) + ext;
    OutFileBWT = fopen(OutFileBWT_s.c_str(), "wb");
    if (OutFileBWT==NULL){
        std::cerr << "Error opening: " << OutFileBWT_s << std::endl;
        exit (EXIT_FAILURE);
    }
		
	dataTypeNChar i,pos = 1;
	dataTypeNSeq seqN=0;
	dataTypeNChar posN=0;
	dataTypedimAlpha pileN=0;
	dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[SIZEBUFFER];
		
	dataTypeNChar numcharBWT=1;
	dataTypeNChar numWrite=1;

	int toRead = SIZEBUFFER;
	while(numcharBWT>0){
		if(freq[alphaInverse[currentPile]] - (pos - 1) < SIZEBUFFER){
			toRead = freq[alphaInverse[currentPile]] - (pos - 1);
		}
		numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),toRead,bwt);
		numWrite = fwrite(BWTbuffer, sizeof(uchar), toRead, OutFileBWT);
		assert(numcharBWT == numWrite);

		for (i=0; i<numcharBWT; i++){
			if(BWTbuffer[i] == TERMINATE_CHAR){
				empty_number--;
				b[pos-1] = 1;
				numchar = fread (&seqN, sizeof(dataTypeNSeq), 1 , eof);
				numchar = fread (&posN, sizeof(dataTypeNChar), 1 , eof);
				numchar = fread (&pileN, sizeof(dataTypedimAlpha), 1 , eof);
				assert(numchar==1);
				fwrite(&seqN, sizeof(dataTypeNSeq), 1, InfoFile);					
			}
			
			//counting the number of occurrences in BWT of the currentPile
			tableOcc[(unsigned int)currentPile][alpha[(unsigned int)BWTbuffer[i]]]++;

			pos++;
		}
	}
	fclose(OutFileBWT);

	rrr_vector<> rrrb(b);
	string bvName = filename + "_bv_" + to_string(currentPile) + ext;
	store_to_file(rrrb,bvName);

	assert(empty_number==0);
	fseek(eof,sizeof(dataTypeNSeq), SEEK_SET);

    for (currentPile = 1 ; currentPile < sizeAlpha; ++currentPile) {
            
        assert(freq[alphaInverse[currentPile]] > 0);
		bit_vector b;
		b = bit_vector(freq[alphaInverse[currentPile]],0);

		//Open output file partial bwt
        string OutFileBWT_s = filename + "_bwt_" + to_string(currentPile) + ext;
        OutFileBWT = fopen(OutFileBWT_s.c_str(), "wb");
        if (OutFileBWT==NULL){
            std::cerr << "Error opening: " << OutFileBWT_s << std::endl;
            exit (EXIT_FAILURE);
        }
		pos = 1;
		dataTypeNSeq seqN=0;
		numcharBWT=1;
		toRead = SIZEBUFFER;
		while(numcharBWT>0){

			if(freq[alphaInverse[currentPile]] - (pos - 1) < SIZEBUFFER){
				toRead = freq[alphaInverse[currentPile]] - (pos - 1);
			}

			numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),toRead,bwt);
			numWrite = fwrite(BWTbuffer, sizeof(uchar), toRead, OutFileBWT);
			assert(numcharBWT == numWrite);

			for (i=0; i<numcharBWT; i++){

				if(BWTbuffer[i] == TERMINATE_CHAR){
					b[pos-1] = 1;
					numchar = fread (&seqN, sizeof(dataTypeNSeq), 1 , eof);
					numchar = fread (&posN, sizeof(dataTypeNChar), 1 , eof);
					numchar = fread (&pileN, sizeof(dataTypedimAlpha), 1 , eof);
					assert(numchar==1);
					fwrite(&seqN, sizeof(dataTypeNSeq), 1, InfoFile);					
				}
			
				//counting the number of occurrences in BWT of the currentPile
				tableOcc[(unsigned int)currentPile][alpha[(unsigned int)BWTbuffer[i]]]++;

				pos++;
			}
		}

		fclose(OutFileBWT);

		rrr_vector<> rrrb(b);
		string bvName = filename + "_bv_" + to_string(currentPile) + ext;
		store_to_file(rrrb,bvName);
	}


	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++){
		for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++){
			fwrite(&(tableOcc[j][h]), sizeof(dataTypeNChar), 1, InfoFile);
		}
	}
	fclose(InfoFile);
	fclose(bwt);
	fclose(eof);


	std::cout << "\nPrint TableOcc:\n";
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
            std::cout << tableOcc[j][h] << "\t";
		std::cout << std::endl;
    }
	
	int status;
	#if KEEPEBWT == 0
		//Remove the auxiliary file .ebwt
		string ebwtFileName = filename + ".ebwt";
		status = remove(ebwtFileName);
		if (status != 0) {
			cerr << "Error deleting file" << ebwtFileName << endl;
		}
	# endif
	//Remove EOFpos file
	status = remove(eof_s);
    if (status != 0) {
        cerr << "Error deleting file" << eof_s << endl;
    }

	return 1;
}



dataTypeNChar buildFreq(string fileName, dataTypeNChar empty_number) {
   
    //Open BWT file
    string fnBWT = string(fileName) + ".ebwt";
   
    FILE *InBWT = fopen(fnBWT.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << fnBWT << "!" << std::endl;
        exit (EXIT_FAILURE);
    }
    fseek(InBWT, 0, SEEK_SET);
    
	for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
		freq[z]=0;
	freq[SIZE_ALPHA-1]=0;
	
    //Reading to find the alphabet
    std::cerr << "Find the alphabet by reading " << fnBWT << " file" << std::endl;

	dataTypeNChar i;
	dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[SIZEBUFFER];
	dataTypeNChar numcharBWT=1;
	dataTypeNChar totalCharRead=0;

	while(numcharBWT>0){

		numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),SIZEBUFFER,InBWT);
		
		for (i=0; i<numcharBWT; i++){
			freq[(unsigned int)(BWTbuffer[i])]+=1;
			totalCharRead++;
		}
	}

	if(freq[TERMINATE_CHAR]==empty_number){
		std::cerr << "ERROR: The end-marker must be " << TERMINATE_CHAR << endl;
		std::cerr << "If you want to use a different end-marker, set the parameter TERMINATE_CHAR in Parameters.h" << endl;
		exit(1);
	}
	fclose(InBWT);
    
	//set alpha and alphaInverse
	alphaInverse = new dataTypedimAlpha[SIZE_ALPHA];
	sizeAlpha=0;
	for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
		if (freq[i] > 0) {
			alpha[i] = sizeAlpha;
			alphaInverse[sizeAlpha]=i;
			sizeAlpha++;
		}
	if (freq[SIZE_ALPHA-1] > 0) {
		alpha[SIZE_ALPHA-1] = sizeAlpha;
		alphaInverse[sizeAlpha]=SIZE_ALPHA-1;
		sizeAlpha++;
	}

	std::cerr << "\nFrom .ebwt file:\n";
	std::cerr << "\tNumber of symbols in the input file: " << totalCharRead << "\n";
	std::cerr << "\tSize alpha: " << (int)sizeAlpha << "\n";

	return totalCharRead;
	
}


dataTypeNChar remove_empty_symbols(string fileName){

	dataTypeNChar bwt_length;
	dataTypeNChar empty_number;

	//Open the file containing the number of symbols in the bwt and the number of empty-word symbols contained within
	string emptyName = fileName+".empty.info";
	FILE* emptyFile = fopen(emptyName.c_str(),"rb");
	if (emptyFile == NULL) {
		std::cerr << "Error opening \"" << emptyName << "\" file"<< std::endl;
		exit (1);
	}
	dataTypeNChar* empty_info = new dataTypeNChar[2];
	if(fread(empty_info,sizeof(dataTypeNChar),2,emptyFile) == 2){
		bwt_length = empty_info[0];
		empty_number = empty_info[1];
	}
	else{
		std::cerr << "Error reading \"" << emptyName << "\" file. File is smaller than expected"<< std::endl;
		exit (1);
	}

	fclose(emptyFile);

	//Compute the number of characters of the bwt to be kept
	dataTypeNChar toKeep = bwt_length - empty_number;


    string bwtFileName = string(fileName) + ".ebwt";
    FILE *bwtFile = fopen(bwtFileName.c_str(), "r");
	if (bwtFile == NULL) {
		std::cerr << "Error opening \"" << bwtFileName << "\" file"<< std::endl;
		exit (1);
	}

    string ebwtFileName = string(fileName) + ".tmp.ebwt";
    FILE *ebwtFile = fopen(ebwtFileName.c_str(), "w");
	if (ebwtFile == NULL) {
		std::cerr << "Error opening \"" << ebwtFileName << "\" file"<< std::endl;
		exit (1);
	}


	dataTypeNChar i;
	dataTypeNChar totalCharRead=0;
	dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[SIZEBUFFER];
	dataTypeNChar numcharBWT=1;
	char terminate_char=TERMINATE_CHAR;

	while(numcharBWT>0){

		//read the bwt in chunks of SIZEBUFFER
		numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),SIZEBUFFER,bwtFile);

		for (i=0; i<numcharBWT; i++){

			//stop reading if reaching the symbols corresponding to the (artificially added) EMPTY_CHARs
			if (totalCharRead < toKeep){

				//if the current character is the EMPTY_CHAR, change it to the end-of-string character, otherwise leave it unchanged
				if (BWTbuffer[i] == EMPTY_CHAR){
					fwrite(&(terminate_char), sizeof(dataTypedimAlpha), 1, ebwtFile);
				}
				else{
					fwrite(&(BWTbuffer[i]), sizeof(dataTypedimAlpha), 1, ebwtFile);
				}
				totalCharRead++;
			}
			else{
				numcharBWT=0;
			}
		}
	}
	
	fclose(bwtFile);
	fclose(ebwtFile);

	//Remove the auxiliary file .empty.info
    int status = remove(emptyName);
    if (status != 0) {
        cerr << "Error deleting file" << emptyName << endl;
    }

	//Move the auxiliary file .tmp.ebwt into .ebwt
    status = rename(ebwtFileName,bwtFileName);
    if (status != 0) {
        cerr << "Error renaming file" << ebwtFileName << endl;
    }
	

	return empty_number;
}
