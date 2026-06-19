
#ifndef EDSBWT_hpp
#define EDSBWT_hpp

/**
 ** Citations:
 **
 ** TITLE: The Burrows-Wheeler transform of an elastic-degenerate string
 ** AUTHORS: Lapo Cioni, Veronica Guerrini and Giovanna Rosone
 ** ICTCS 2024
 ** 
 **
 ** Copyright (c) 2024, Lapo Cioni, Veronica Guerrini and Giovanna Rosone.  All rights reserved.
 ** Use of this source code is governed
 ** by a MIT license that can be found in the LICENSE file. 
 **
 **/


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
//#include <stdlib.h>

#include "Parameters.h" // Defines ulong and uchar.
#include "Sorting.h" // Defines ulong and uchar.
#include <string.h>

#include <deque>

////////2024///////
#include <map>
#include <sdsl/bit_vectors.hpp>
////////////////


#define BUFFERSIZE 1024

#define SIZE_ALPHA 256

#define DIMBLOCK  1024     // 1048576

using std::string;
using std::vector;
using std::map;
using namespace sdsl;

class EDSBWT
{
public:
    EDSBWT(string, int, int);
    EDSBWT(string, string, int, int);
    ~EDSBWT();

    string  fileOutBwt;
    string  fileOutCyc;
    string  ext;
    string  fileInputBase;
    bool resourcesLoaded;
    
	std::vector<rangeElement> vectRangeDollarPile;
	std::vector<rangeElement> vectRangeOtherPile;
	dataTypedimAlpha symbPile;
    
    dataTypeNChar SIZEBUFFERcycFiles;
    
    dataTypelenSeq lengthRead;    //Length of each text
    dataTypeNChar lengthTot;   //Total length of all texts without $-symbols
    dataTypeNChar lengthTot_plus_eof;   //Total length of all texts with $-symbols

    dataTypeNSeq nText;   //number total of texts in filename1
    dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
    
    int numthreads;
    
    dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
    dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
    dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
    dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char
    
	//Added/modified by Lapo 06-06-2024
	bit_vector rrrb;
	dataTypeNSeq first_symbol_index;

	//Added by Vero 15_05_2024
	std::vector< rrr_vector<> > EOFpos; 
	std::vector< rrr_vector<>::rank_1_type > EOFpos_rank; 
	dataTypeNChar** EOF_ID;

    std::vector <bool> vectInsTexts;
    dataTypeNSeq textToInsert;
    vector <dataTypeNChar> vectSizeCurrentPile;

    vector< vector< vector<dataTypeNChar> > > vectorOcc;
    vector <dataTypeNChar> numBlocksInPartialBWT;
    
    int recoverInfo(string);

    int splitIntoPartial(string, int);
    
    int buildFreq(string);
  
    int decodeBCRmultipleReverse(string);
    
    //dataTypeNChar rankManySymbolsByVector(FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *, FILE &);
    
    int computeVectorOcc(string filename);
    
    int convertFromCycFileToFastaOrFastq( string );
	
	////////2024///////
	float computeSimilarityFromEDS(string fileEDS);
	int backwardSearch(string fileInput, string fileOutDecode, dataTypeNSeq n_kmer, string kmers, dataTypelenSeq lenKmer, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
    int computeCountersInSingleInterval(FILE *InFileBWT, dataTypeNChar toRead, dataTypedimAlpha currentPile,dataTypeNChar * counters);

    int updateSingleInterval(std::vector<rangeElement> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, uchar symbol, dataTypeNChar * counters, dataTypeNChar *numBlockCounter, dataTypeNChar * contInCurrentBlock, dataTypeNChar toRead, uchar *bufferBlock);
	int updateIntervals(std::vector<rangeElement> &vectRange, string fileInput, string fileOutDecode, uchar symbol, dataTypedimAlpha currentPile);
	//int updateSingleInterval(std::vector<rangeElement> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, uchar symbol);
	dataTypeNChar rankManySymbols(FILE & InFileBWT, dataTypeNChar *counters, dataTypeNChar toRead, uchar *foundSymbol, uchar *bufferBlock);
	int link(rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
	void dollars_in_interval(dataTypelenSeq symbPile, std::deque<dataTypeNSeq> &output,dataTypeNChar i,dataTypeNChar j, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
	rangeElement preceding_dollars_finder(dataTypeNSeq i, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
	
	dataTypeNChar a;
	

	
	#if RECOVERBW==1
		int findMultipleDollarsBackward(std::vector<rangeElementBW> &vectRange, string fileInput,string fileOutDecode,dataTypeNSeq n_kmer,rank_support_v<1> &rb_1,bit_vector::select_1_type &bsel_1);
		int updateSingleIntervalBW(std::vector<rangeElementBW> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, dataTypeNChar * counters, dataTypeNChar *numBlockCounter, dataTypeNChar * contInCurrentBlock, dataTypeNChar toRead, uchar *bufferBlock);
		int findBlockToReadBWT(dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock);
		void printElementBW(std::vector<rangeElementBW> &vectRange);
	#endif
	void indices_of_dollars_in_interval(std::ostream& searchOutput,dataTypelenSeq symbPile, dataTypeNChar i, dataTypeNChar j, dataTypeNSeq n_kmer, dataTypeNChar occPos, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
	
	int rankInverseManyByVector (string fileInput, dataTypeNSeq numKmersInput, uchar *toFindSymbols, uchar *bufferBlock);
	int update_Pos_Pile_Blocks(dataTypeNChar *toRead, dataTypeNChar *numBlock, dataTypedimAlpha currentPile, uchar toFindSymbol);
	int update_Pos_Pile(rangeElement *tripla);

	////////////////////
    
private:

    int findBlockToRead(dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar *);
    
    FILE * openFilePartialIn(string filename, dataTypedimAlpha currentPile);
    dataTypeNChar readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT);
    dataTypeNChar writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT);
    int closeFilePartial(FILE * pFile);
	
	void print (std::vector<rangeElement> &vectRange);
	void print_interval_number ();

    
};




#endif /* EDSBWT_hpp */
