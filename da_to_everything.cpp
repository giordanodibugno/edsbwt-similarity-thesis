#include <iostream>
#include <assert.h>
#include <sdsl/bit_vectors.hpp>
#include "Parameters.h"

#define END_MARKER '#' //end-marker used for the string collection 

using namespace std;
using namespace sdsl;

dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char

dataTypeNChar buildFreq(string fileName, dataTypeNChar BWT_seqs);
int da_to_everything(string filename, dataTypeNChar BWT_length);
dataTypeNChar remove_empty_symbols(string inputName);



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
	dataTypeNChar BWT_seqs =remove_empty_symbols(inputName);
	dataTypeNChar BWT_length = buildFreq(inputName,BWT_seqs);
	da_to_everything(inputName, BWT_length);
	return 1;
}



int da_to_everything(string filename, dataTypeNChar BWT_length){

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
	

	//Open document array
	string da_s = filename + ".4.da";
	FILE* da = fopen(da_s.c_str(),"rb");
	if (da == NULL) {
		std::cerr << "Error opening \"" << da_s.c_str() << "\" file"<< std::endl;
		exit (1);
	}


	//Open ebwt
	//n.b. The file does not contain the tail, so the length is different from document array.
	string bwt_s = filename + ".ebwt";
	FILE* bwt = fopen(bwt_s.c_str(),"r");
	if (bwt == NULL) {
		std::cerr << "Error opening \"" << bwt_s << "\" file"<< std::endl;
		exit (1);		
	}

	FILE* OutFileBWT;
	//FILE* dollarIndexOut;
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
	
	size_t currentPile;

    for (currentPile = 0 ; currentPile < sizeAlpha; ++currentPile) {
            
        assert(freq[alphaInverse[currentPile]] > 0);
		bit_vector b;
		b = bit_vector(freq[alphaInverse[currentPile]],0);
		//vector<dataTypeNSeq> dollarIndex;



        /*Open output file dollarIndex
        string dollarIndex_s = filename + "_EOFpos_" + to_string(currentPile) + ext;
        dollarIndexOut = fopen(dollarIndex_s.c_str(), "wb");
        if (dollarIndexOut == NULL){
            std::cerr << "Error opening: " << dollarIndex_s << std::endl;
            exit (EXIT_FAILURE);
        }*/

		//Open output file partial bwt
        string OutFileBWT_s = filename + "_bwt_" + to_string(currentPile) + ext;
        OutFileBWT = fopen(OutFileBWT_s.c_str(), "wb");
        if (OutFileBWT==NULL){
            std::cerr << "Error opening: " << OutFileBWT_s << std::endl;
            exit (EXIT_FAILURE);
        }
		


		dataTypeNChar i;
		dataTypeNChar pos = 1;
		dataTypeNSeq* DAbuffer = new dataTypeNSeq[SIZEBUFFER];
		dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[SIZEBUFFER];
		dataTypeNChar numcharBWT=1;
		dataTypeNChar numcharDA=1;
		dataTypeNChar numWrite=1;

		int toRead = SIZEBUFFER;
		while(numcharBWT>0 and numcharDA>0){

			if(freq[alphaInverse[currentPile]] - (pos - 1) < SIZEBUFFER){
				toRead = freq[alphaInverse[currentPile]] - (pos - 1);
			}

			numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),toRead,bwt);
			numcharDA = fread(DAbuffer,sizeof(dataTypeNSeq),toRead,da);
			assert(numcharBWT == numcharDA);

			numWrite = fwrite(BWTbuffer, sizeof(uchar), toRead, OutFileBWT);
			assert(numcharBWT == numWrite);

			for (i=0; i<numcharBWT; i++){

				if(BWTbuffer[i] == TERMINATE_CHAR){
					b[pos-1] = 1;
					//fwrite(&(DAbuffer[i]), sizeof(dataTypeNSeq), 1, dollarIndexOut);					
					fwrite(&(DAbuffer[i]), sizeof(dataTypeNSeq), 1, InfoFile);					
				}
			
				//counting the number of occurrences in BWT of the currentPile
				tableOcc[(unsigned int)currentPile][alpha[(unsigned int)BWTbuffer[i]]]++;

				pos++;
			}
		}

		//fclose(dollarIndexOut);
		fclose(OutFileBWT);

		rrr_vector<> rrrb(b);
		string bvName = filename + "_bv_" + to_string(currentPile) + ext;
		store_to_file(rrrb,bvName);
	}

	
	/*
	string tableOccOut_s = filename + "_tableOcc" + ext;
	FILE* tableOccOut = fopen(tableOccOut_s.c_str(),"wb");
	if (tableOccOut == NULL){
		std::cerr << "Error opening: " << tableOccOut_s << std::endl;
		exit (EXIT_FAILURE);
	}
	*/

	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++){
		for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++){
			//fwrite(&(tableOcc[j][h]), sizeof(dataTypeNChar), 1, tableOccOut);
			fwrite(&(tableOcc[j][h]), sizeof(dataTypeNChar), 1, InfoFile);
		}
	}
	fclose(InfoFile);


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

	//Remove the auxiliary file .4.da
    status = remove(da_s);
    if (status != 0) {
        cerr << "Error deleting file" << da_s << endl;
    }


	return 1;
}










dataTypeNChar buildFreq(string fileName, dataTypeNChar BWT_seqs) {
	
    string fnBWT = string(fileName) + ".ebwt";
   
    FILE *InBWT = fopen(fnBWT.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << fnBWT << "!" << std::endl;
        exit (EXIT_FAILURE);
    }
    fseek(InBWT, 0, SEEK_SET);

	dataTypeNChar numEle=0;
    
	for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
		freq[z]=0;
	freq[SIZE_ALPHA-1]=0;


    //First reading in order to find the alphabet
    std::cerr << "Find the alphabet by reading " << fnBWT << " file" << std::endl;



	dataTypeNChar i;
	dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[SIZEBUFFER];
	dataTypeNChar numcharBWT=1;

	while(numcharBWT>0){

		numcharBWT = fread(BWTbuffer,sizeof(dataTypedimAlpha),SIZEBUFFER,InBWT);

			for (i=0; i<numcharBWT; i++){
				freq[(unsigned int)(BWTbuffer[i])]+=1;
				numEle++;
			}
	}

	if(freq[TERMINATE_CHAR]!= BWT_seqs){
		std::cerr << "ERROR: The end-marker must be " << TERMINATE_CHAR << endl;
		std::cerr << "If you want to use a different end-marker, set the parameter TERMINATE_CHAR in Parameters.h" << endl;
		exit(1);
	}

    
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
	std::cerr << "\tNumber of symbols in the input file: " << numEle << "\n";
	std::cerr << "\tSize alpha: " << (int)sizeAlpha << "\n";

	return numEle;
}




dataTypeNChar remove_empty_symbols(string fileName){

	dataTypeNChar bwt_length;
	dataTypeNChar empty_number;
	dataTypeNChar numSeq=0;

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


    string bwtFileName = string(fileName) + ".bwt";
    FILE *bwtFile = fopen(bwtFileName.c_str(), "r");
	if (bwtFile == NULL) {
		std::cerr << "Error opening \"" << bwtFileName << "\" file"<< std::endl;
		exit (1);
	}

    string ebwtFileName = string(fileName) + ".ebwt";
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
				else if (BWTbuffer[i] == END_MARKER){
					fwrite(&(terminate_char), sizeof(dataTypedimAlpha), 1, ebwtFile);
					numSeq++;
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

	//Remove the auxiliary file .bwt
    status = remove(bwtFileName);
    if (status != 0) {
        cerr << "Error deleting file" << bwtFileName << endl;
    }


	return numSeq+empty_number;
}
