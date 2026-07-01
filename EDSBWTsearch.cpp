//
//  EDSBWT.cpp
//  EDSBWT
//


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
 
#include "EDSBWTsearch.hpp"
#include "Sorting.h"
#include "Parameters.h"

#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <string.h>     // std::string, std::to_string
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>      // std::stringstream

#include <math.h>

#include <omp.h>

#include <sys/resource.h>

////////2024///////
#include <map>
#include <sdsl/bit_vectors.hpp>
#include "malloc_count/malloc_count.h"
////////////////

#ifndef RECOVER_INFO
	#define RECOVER_INFO 1
#else
	#define RECOVER_INFO 0
#endif

using namespace std;
using namespace sdsl;

static vector<vector<string> > readEDSByLocus(string filename) // Parser di file .eds
{
	vector<vector<string> > loci;
	ifstream edsFile(filename);
	if (!edsFile.is_open()) {
		cerr << "ERROR opening EDS file " << filename << endl;
		exit(1);
	}

	vector<string> currentLocus;
	string currentString;
	bool insideLocus = false; // Variabile booleana che indica se in un istante stiamo analizzando dentro un locus
	char current; // Carattere corrente

	while (edsFile.get(current)) {
		if (current == '\n' || current == '\r' || current == '\t' || current == ' ') {
			continue; // Skip
		}

		if (current == '{') {
			if (insideLocus) {
				cerr << "ERROR parsing EDS file " << filename << ": nested '{' found" << endl;
				exit(1);
			}
			insideLocus = true;
			currentLocus.clear();
			currentString.clear();
		}
		else if (current == '}') {
			if (!insideLocus) {
				cerr << "ERROR parsing EDS file " << filename << ": unmatched '}' found" << endl;
				exit(1);
			}
			if (currentString == string(1, EMPTY_CHAR_EDS)) {
				currentString.clear();
			}
			currentLocus.push_back(currentString);
			loci.push_back(currentLocus);
			currentString.clear();
			currentLocus.clear();
			insideLocus = false;
		}
		else if (current == ',') {
			if (!insideLocus) {
				cerr << "ERROR parsing EDS file " << filename << ": comma outside a locus" << endl;
				exit(1);
			}
			if (currentString == string(1, EMPTY_CHAR_EDS)) {
				currentString.clear();
			}
			currentLocus.push_back(currentString);
			currentString.clear();
		}
		else {
			if (!insideLocus) {
				cerr << "ERROR parsing EDS file " << filename << ": symbol outside a locus" << endl;
				exit(1);
			}
			currentString.push_back(current);
		}
	}

	if (insideLocus) {
		cerr << "ERROR parsing EDS file " << filename << ": missing closing '}'" << endl;
		exit(1);
	}

	return loci;
}

EDSBWT::EDSBWT (string fileInput, string fileEDS, int mode, int num_threads)
	: EDSBWT(fileInput, mode, num_threads)
{
	computeSimilarityFromEDS(fileEDS);
}

EDSBWT::EDSBWT (string fileInput, int mode, int num_threads)
{
    
	std::cerr << "Backward Search\n";
	ext = ".aux";
	fileInputBase = fileInput;
	resourcesLoaded = false;
	
	cout << "DEBUG: " << DEBUG << endl;
	
	#if RECOVER_INFO
		#if DEBUG==1
			fprintf(stderr, "##BEFORE recoverInfo\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
		#endif
		assert(recoverInfo(fileInput)==1);
		//cout << "Lettura corretta?" << endl;
		//exit(1);
	#else
	if (mode == 1){
		time_t startS,endS;
		double difS;
		//run buildFreq
		assert (buildFreq(fileInput) == 1);

		#if DEBUG == 1
		time (&startS);
		//run splitIntoPartial to compute BWT-partial
		std::cerr << "Start splitIntoPartial " << startS << " seconds\n";
		#endif
		
		assert (splitIntoPartial(fileInput,1) == 1);
		
		#if DEBUG == 1
		time (&endS);
		difS = difftime (endS,startS);
		std::cerr << "End splitIntoPartial " << endS << " seconds\n";
		std::cerr << "splitIntoPartial tooks " << difS << " seconds\n";
		#endif

		#if (DEBUG == 1) //|| (verboseDecode == 1)
			std::cerr << "i" << "\t" << "freq" << "\t" << "Code" << "\t" << "ASCII" << "\n";
			dataTypedimAlpha mmm=0;
			for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; i++)
				if (freq[i] > 0) {
					std::cerr << i << "\t" << freq[i] << "\t" << (unsigned int)alpha[i] << "\t" << (unsigned int)alphaInverse[mmm] << "\n";
					mmm++;
				}

			std::cerr << "TableOcc: "  << "\n";
			for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
				std::cerr << int(g)  << ":\t";
				for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
					std::cerr << tableOcc[g][j]  << "\t";
				std::cerr << "\n";
			}
		#endif
	}
	#endif

	#if (DEBUG == 1)
	fprintf(stderr, "##BEFORE computeVectorOcc\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	std::cerr << "For the computation of the new positon useful for unbuildBCR we use the vector of the occurrences of " << DIMBLOCK << " size" << std::endl;
	#endif

	//time_t startC,endC;
	//double difC;
	//time (&startC);
	//std::cerr << "\nStart computeVectorOcc " << startC << " seconds\n";

	//run computeVectorOcc to set vectorOcc    

	assert ( computeVectorOcc(fileInput) == 1);
	
	
	//time (&endC);
	//difC = difftime (endC,startC);
	//std::cerr << "End    computeVectorOcc " << endC << " seconds\n";
	//std::cerr << "computeVectorOcc tooks " << difC << " seconds\n";

	#if DEBUG
		for (int x = 0 ; x < sizeAlpha; x++)  {
			for(int z = 0; z < sizeAlpha ; z++)   {
				cout << "vectorOcc[" << x << "][" << z << "]=" << vectorOcc[x][z][0] << "\t";
			}
			cout << endl;
		}
	#endif

    
	#if RECOVER_INFO==0
		char *fileEndPos = new char[256];
		sprintf (fileEndPos,"%s%s",fileInput.c_str(),".EOFpos");
		fprintf(stderr, "##BEFORE MAP\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
		m = map_create(fileEndPos);	//n.b. non usiamo più map_create, quindi credo che le cose in RECOVER_INFO debbano essere rimosse tutte.
	#endif
	
	char *fileBitVector = new char[256];
	sprintf (fileBitVector,"%s%s",fileInput.c_str(),".bitvector");
	
	#if DEBUG == 1
	cout<<"The bitvector is in: "<<fileBitVector<<endl;
	fprintf(stderr, "##BEFORE rank_1_type\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	
	if (load_from_file(rrrb,fileBitVector) != 1) {
		std::cerr << "Error opening \"" << fileBitVector << "\" file"<< std::endl;
		exit (1);
	}
	rank_support_v<1> rb_1(&rrrb);
	
	#if DEBUG == 1
		fprintf(stderr, "##BEFORE select_1_type\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	bit_vector::select_1_type bsel_1(&rrrb);
	dataTypeNSeq numLoci = rb_1(rrrb.size());
	if (numLoci == 0) {
		std::cerr << "Error: the EDS bitvector contains no loci" << std::endl;
		exit(1);
	}
	first_symbol_index = numLoci > 1 ? bsel_1(2)-1 : rrrb.size()-1;

	#if DEBUG == 1
	fprintf(stderr, "##AFTER select_1_type\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
		cout<<"This is the bitvector"<<endl;
		for (dataTypeNSeq ttt=0; ttt < rrrb.size() ; ttt++ )
			cout<< rrrb[ttt] << "\t";
		cout << endl;
	#endif

	delete [] fileBitVector;
	resourcesLoaded = true;
	
}

float EDSBWT::computeSimilarityFromEDS(string fileEDS)
{
	rank_support_v<1> rb_1(&rrrb);
	bit_vector::select_1_type bsel_1(&rrrb);

	#if RECOVERBW==1
		//string searchOutput_s = fileInputBase + ".output.csv";
		string searchOutput_s = fileEDS + "output.csv";
		std::ofstream searchOutput;
		searchOutput.open(searchOutput_s,ios::out);
		if(searchOutput.is_open()){
//			searchOutput << "Pattern_#" << "\t" << "word_index" << "\t" << "segment_index" << "\t" << "WordInSeg_index" << "\t" << "position_in_string\n";
			searchOutput << "#Pat" << "\t" << "$_i" << "\t" << "D[i]" << "\t" << "S_j" << "\t" << "S_j[r]\n";
			searchOutput.close();
		}
		else{
			cerr << "ERROR opening file " << searchOutput_s << " to write output\n";
			exit(1);
		}
	#endif

	#if DEBUG == 1
		fprintf(stderr, "##BEFORE BACKWARD\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
		std::cerr << "backwardSearch\n";
	#endif

	vector<vector<string> > loci = readEDSByLocus(fileEDS);
	vector<float> locusSimilarities;
	dataTypeNChar len_match;
	std::string kmer;
	dataTypeNChar globalStringId=0;
	dataTypeNChar lenKmer=0;
	float totalNum = 0.0; // Numeratore frazione finale
	float totalDen = 0.0; // Denominatore frazione finale

	cout << "Loci read from EDS file: " << loci.size() << endl;

	for (dataTypeNSeq locusId=0; locusId<loci.size(); locusId++) {
		float locusNum = 0.0;
		float locusDen = 0.0;
		cout << "Locus " << locusId + 1 << " contains " << loci[locusId].size() << " strings" << endl;

		for (dataTypeNSeq stringId=0; stringId<loci[locusId].size(); stringId++) {
			kmer = loci[locusId][stringId];
			lenKmer = kmer.length();

// 			cout << "Pattern: " << kmer << " of length " << lenKmer << endl;

			if (lenKmer == 0) {
				cout << "len_match=0 (empty string skipped)" << endl << endl;
				continue;
			}

			locusDen += lenKmer;
			totalDen += lenKmer;

			#if DEBUG==1
			time_t startI,endI;
			double difI;
				time (&startI);
			#endif
			
			len_match=backwardSearch(fileInputBase.c_str(), fileEDS.c_str(), globalStringId+1, kmer, lenKmer, rb_1, bsel_1);
			
			locusNum += len_match;
			totalNum += len_match;
			std::cout << "len_match=" << len_match << endl;
			std::cout << endl;

			#if DEBUG==1
				time (&endI);
				difI = difftime (endI,startI);
				std::cerr << "End backwardSearch " << endI << " seconds\n";
				std::cerr << "backwardSearch tooks " << difI << " seconds\n";
			#endif
			globalStringId++;
		}

		float locusSim = 0.0;
		if (locusDen > 0.0) {
			locusSim = locusNum/locusDen;
		}
		locusSimilarities.push_back(locusSim);
		std::cout << "Similarita locus " << locusId + 1 << " (" << locusNum << "/" << locusDen << "):" << locusSim << endl;
		std::cout << endl;
	}

	float sim = 0.0;
	if (totalDen > 0.0) {
		sim = totalNum/totalDen;
	}
	std::cout << "Array similarita loci: [";
	for (dataTypeNSeq locusId=0; locusId<locusSimilarities.size(); locusId++) {
		if (locusId > 0) {
			std::cout << ", ";
		}
		std::cout << locusSimilarities[locusId];
	}
	std::cout << "]" << endl;
	std::cout << "Similarita totale EDS ("<<totalNum<<"/"<<totalDen<<"):" << sim << endl;

	return sim;
}




#if RECOVERBW==1
int EDSBWT::findMultipleDollarsBackward(std::vector<rangeElementBW> &vectRange, string fileInput,string fileOutDecode,dataTypeNSeq n_kmer,rank_support_v<1> &rb_1,bit_vector::select_1_type &bsel_1)	{
	#if DEBUG==1
	time_t start,end;
	time (&start);
	#endif

	//string searchOutput_s = fileInput + ".output.csv";
	string searchOutput_s = fileOutDecode + "output.csv";

	std::ofstream searchOutput;
	searchOutput.open(searchOutput_s, std::ios_base::app);

	
	//searchOutput << "There is at least one occurrence of the pattern for each line, and the pattern begins in some position of the corresponding word. " 
	//	<< "word index is the index of the word, starting from 0. degSymb index is the index of the degenerate symbol that contains that word, starting from 1. "
	//	<< "inSeg index is the the index of the word with respect to the degenerate symbol in which it is contained, starting from 0.\n";
	//searchOutput << "Pattern #" << n_kmer << "\t" << "word index" << "\t" << "segment index" << "\t" << "WordInSeg index" << "\t" << "position in word" << endl;
	//searchOutput << "\nPattern # " << n_kmer << endl;
	//	char *filenameIn = new char[128];
	//	char *filename = new char[9];
		string filenameIn;
		const char *ext = ".aux";
		FILE *InFileBWT;
				
		#if DEBUG == 1
		dataTypeNSeq numTotKmers = vectRange.size();
		std::cerr << "findMultipleDollarsBackward: We want to compute the seqID of " << numTotKmers  << " sequences." << std::endl;
		#endif
		
		std::vector< rangeElementBW > vectRangeCopy;
		rangeElementBW rangeEle;
		
  
	dataTypeNChar * countersDiff = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT 
	dataTypeNChar * countersStart = new dataTypeNChar[sizeAlpha];
	dataTypeNChar * countersEnd = new dataTypeNChar[sizeAlpha];
	
	dataTypeNChar numBlockCounterStart = 0, numBlockCounterEnd= 0, numBlock=0;   //number of the blocks read
	dataTypeNChar contInCurrentBlockStart = 0, contInCurrentBlockEnd= 0;    //number of the symbols read
		
	uchar *bufferBlock = new uchar[DIMBLOCK];
	uchar foundSymbol = '\0';  //here, it is not useful

		
	dataTypeNSeq j = 0;
	dataTypeNChar occPos = 0;
	
	
	while (!vectRange.empty()) {
		quickSort(vectRange);
		
		#if DEBUG == 1
			std::cerr << "\n findMultipleDollarsBackward - Their (sorted) intervals are:" << std::endl;
			printElementBW(vectRange);
		#endif

			
		j=0;
		while (j < vectRange.size()) {   //If there is at least one TERMINATE_CHAR symbol
			
			//For each symbol in the kmer we have to update vectRange (First and Last vectors)
			dataTypedimAlpha currentPile = vectRange[j].pileN;
			filenameIn = fileInput + "_bwt_" + to_string((int)(currentPile))+ext;


			InFileBWT = fopen(filenameIn.c_str(), "rb");
			if (InFileBWT==NULL) {
				std::cerr << "findMultipleDollarsBackward: BWT file " << filenameIn << " j= " << (unsigned int)j << ": Error opening " << std::endl;
				exit (EXIT_FAILURE);
			}

			dataTypeNSeq k=j;
			for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
				countersDiff[i]=0;
				countersStart[i]=0;
				countersEnd[i]=0;
			}
			numBlockCounterStart = 0, numBlock=0;   //number of the blocks read
			contInCurrentBlockStart = 0;    //number of the symbols read
			numBlockCounterEnd = 0, numBlock=0;   //number of the blocks read
			contInCurrentBlockEnd = 0;    //number of the symbols read
			
			if (vectRange[k].startPosN <= vectRange[k].endPosN) {
				vectRange[k].startPosN --;   //So we compute rank until position First - 1
			
				dataTypeNChar toRead = vectRange[k].startPosN;
				
				if (toRead > 0) {
					//Before, we find the block where toRead position is.
					assert ( findBlockToReadBWT(currentPile, &toRead, &numBlock) == 1);
					
					//First time! Update counters con the symbols up to the previous block
					if (numBlock > 0) {
						for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
							countersStart[r] =  vectorOcc[currentPile][r][(numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
					}					
				}
				//else (toRead == 0) --> numBlock=0
			
				#if DEBUG == 1
					//if (numBlock >= numBlocksInPartialBWT[currentPile])
						std::cerr << "\t First interval: toRead " << toRead << " numBlockCounterStart " << numBlockCounterStart << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";				    
				#endif
				
				//Move file pointer to the block numBlock and update counter at position toRead
				fseek (InFileBWT, numBlock*DIMBLOCK, 0);
				//Count the symbols in the selected block
					
				dataTypeNChar numberRead = rankManySymbols(*InFileBWT, countersStart, toRead, &foundSymbol, bufferBlock);
							#if DEBUG == 1
								std::cerr << "\t After first rank********foundSymbol " << (unsigned int)foundSymbol <<  "\n";
							#endif
				assert (toRead <= numberRead);
				
				numBlockCounterStart = numBlock;
				contInCurrentBlockStart = toRead;
				#if DEBUG == 1				
					std::cerr << "\t UPDATE: numBlockCounterStart " << numBlockCounterStart << " contInCurrentBlockStart " << contInCurrentBlockStart << "\n";
					std::cerr << "countersStart:\t";
					for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
								   std::cerr << " " << countersStart[i];
								std::cerr << "\n";
				#endif
				
				//END
				toRead = vectRange[k].endPosN;
				//Commented since countersEnd already computes its value independently from counterStart
				//if ( (dataTypeNChar)floor((long double)((toRead-1)/DIMBLOCK)) == numBlockCounterStart) {
				//	for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
				//		countersEnd[i] += countersStart[i];  //+1???
				//	}
				//	contInCurrentBlockEnd = contInCurrentBlockStart;
				//}
				
				assert (updateSingleIntervalBW(vectRange, InFileBWT, k, currentPile, countersEnd, &numBlockCounterEnd, &contInCurrentBlockEnd, toRead, bufferBlock) == 1);
				#if DEBUG == 1				
					std::cerr << "\t UPDATED: numBlockCounterEnd " << numBlockCounterEnd << " contInCurrentBlockEnd " << contInCurrentBlockEnd << "\n";
					std::cerr << "countersEnd:\t";
					for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
						std::cerr << " " << countersEnd[i];
					std::cerr << "\n";
				#endif
				
				#if DEBUG == 1
						std::cerr << "findMultipleDollarsBackward: counters (simboli presenti) countersDiff: \t";
				#endif
				
				//Compute new positions
				for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
					countersDiff[i] =  countersEnd[i] - countersStart[i];  //+1???
					
							#if DEBUG == 1
								std::cerr << " " << countersDiff[i];
							#endif				
				}

				
				#if DEBUG == 1
					std::cerr << "\n";
				#endif
				//For end-marker
				if (countersDiff[0] > 0) { //Since at least a dollar belongs to this interval, we call LINK
						rangeEle.startPosN = vectRange[k].startPosN + 1;   //First = c[c] + rank (c, First - 1) + 1
						rangeEle.endPosN = vectRange[k].endPosN;
								//rangeEle.seqN = vectRange[k].seqN;
						rangeEle.pileN = vectRange[k].pileN;

						// this is the LINK CALL.
						indices_of_dollars_in_interval(searchOutput, rangeEle.pileN, rangeEle.startPosN, rangeEle.endPosN, n_kmer, occPos, rb_1, bsel_1);

						//seqID.push_back(rangeEle);
						countersStart[0] = countersEnd[0];  //For the next iteration
						countersEnd[0] = 0;
				}
				//For the other symbols
				for (dataTypedimAlpha i = 1 ; i < sizeAlpha; i++) {
							if (countersDiff[i] > 0) { //Simbolo alpha[i] presente
								rangeEle.pileN=i;
								rangeEle.startPosN=countersStart[i];
								rangeEle.endPosN=countersEnd[i];
								for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {									
									rangeEle.startPosN = rangeEle.startPosN + tableOcc[g][i];
									rangeEle.endPosN = rangeEle.endPosN + tableOcc[g][i];
								}
								//First = c[c] + rank (c, First - 1) + 1
								rangeEle.startPosN ++;  //We must to sum 1 to first
								//rangeEle.seqN = vectRange[k].seqN;
								#if DEBUG == 1
								std::cerr << "--> Insert pileN " << (int)rangeEle.pileN << " , pos (" << rangeEle.startPosN << ", " << rangeEle.endPosN << ")\n";
								#endif
								vectRangeCopy.insert(std::end(vectRangeCopy), rangeEle);
							}
							countersStart[i] = countersEnd[i];  //For the next iteration
							countersEnd[i] = 0;
				}
				numBlockCounterStart = numBlockCounterEnd;
				contInCurrentBlockStart = contInCurrentBlockEnd;
			}
			/////////////////////
			
			#if DEBUG == 1
			std::cerr << "-->Prima del WHILE, numBlockCounterStart " << numBlockCounterStart << " contInCurrentBlockStart: " << contInCurrentBlockStart << " countersStart:\n";
			for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
				std::cerr << " " << countersStart[i];
			}
			#endif	
			
			k++;
			
			//BOTH
			while ((k< vectRange.size()) && (vectRange[k].pileN == currentPile)) {
				
					if (vectRange[k].startPosN <= vectRange[k].endPosN) {
						
						//START
						vectRange[k].startPosN --;   //So we compute rank until position First - 1
						dataTypeNChar toRead = vectRange[k].startPosN;	

						assert (updateSingleIntervalBW(vectRange, InFileBWT, k, currentPile, countersStart, &numBlockCounterStart, &contInCurrentBlockStart, toRead, bufferBlock) == 1);
						
						//END
						toRead = vectRange[k].endPosN;
						//Commented since countersEnd already computes its value independently from counterStart
						//if ( (dataTypeNChar)floor((long double)((toRead-1)/DIMBLOCK)) == numBlockCounterStart) {
						//	for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
						//		countersEnd[i] += countersStart[i];  //+1???
						//	}
						//	contInCurrentBlockEnd = contInCurrentBlockStart;
						//}
						assert (updateSingleIntervalBW(vectRange, InFileBWT, k, currentPile, countersEnd, &numBlockCounterEnd, &contInCurrentBlockEnd, toRead, bufferBlock) == 1);
						
						for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++) {
							countersDiff[i] =  countersEnd[i] - countersStart[i];  //+1???
						}
						numBlockCounterStart = numBlockCounterEnd;
						contInCurrentBlockStart = contInCurrentBlockEnd;
						
						//I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
						//Si potrebbe unire al precedente
						if (countersDiff[0] > 0) { //Ci sono dei dollari, va chiamato link
							rangeEle.startPosN = vectRange[k].startPosN + 1; //Avevamo sottratto 1
							rangeEle.endPosN = vectRange[k].endPosN;
							//rangeEle.seqN = vectRange[k].seqN;
							rangeEle.pileN = vectRange[k].pileN;

							// this is the LINK CALL.
							indices_of_dollars_in_interval(searchOutput, rangeEle.pileN, rangeEle.startPosN, rangeEle.endPosN, n_kmer, occPos, rb_1, bsel_1);

							//seqID.push_back(rangeEle);
							countersStart[0] = countersEnd[0];  //For the next iteration
							countersEnd[0] = 0;
							
						}
						//For each other symbol in the range
						for (dataTypedimAlpha i = 1 ; i < sizeAlpha; i++) {
							if (countersDiff[i] > 0) { //Simbolo alpha[i] presente
								rangeEle.pileN=i;
								rangeEle.startPosN=countersStart[i];
								rangeEle.endPosN=countersEnd[i];
								for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {
									rangeEle.startPosN = rangeEle.startPosN + tableOcc[g][i];
									rangeEle.endPosN = rangeEle.endPosN + tableOcc[g][i];
								}
								//First = c[c] + rank (c, First - 1) + 1
								rangeEle.startPosN ++;  //We must to sum 1 to first
								//rangeEle.endPosN ++;  //We must to sum 1 to last
								//rangeEle.seqN = vectRange[k].seqN;
								vectRangeCopy.insert(std::end(vectRangeCopy), rangeEle);
							}
							countersStart[i] = countersEnd[i];  //For the next iteration
							countersEnd[i] = 0;
						}
						numBlockCounterStart = numBlockCounterEnd;
						contInCurrentBlockStart = contInCurrentBlockEnd;				
					}
					
				k++;
			}
			fclose(InFileBWT);

			j=k;
	
		}
		
		occPos++;  		
	
		vectRangeCopy.shrink_to_fit();
		vectRange.swap(vectRangeCopy);
		vector<rangeElementBW>().swap(vectRangeCopy);   // clear vectRangeCopy reallocating 
	}

	searchOutput.close();
	#if DEBUG==1
	time (&end);
    std::cerr << "findMultipleDollarsBackward: occPos = " << occPos << " and took " << difftime (end,start) << " seconds" << "\n";
	#endif


	//delete[] toFindSymbols;
	//return seqID;
	return 1;
}
#endif


int EDSBWT::backwardSearch(string fileInput, string fileOutDecode, dataTypeNSeq n_kmer, string kmer, dataTypelenSeq lenKmer, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1)
{
	//std::cerr << "Vector size of the occurrences: " << DIMBLOCK  << std::endl;

//#if DEBUG==1
// 	std::cerr << "Pattern: " << kmer << "\n";
//#endif

	//Initialization
	uchar symbol = kmer[lenKmer-1];
	bool symbolInAlphabet = false;
	for (dataTypedimAlpha mm = 0; mm < sizeAlpha; mm++) {
		if (alphaInverse[mm] == symbol) {
			symbolInAlphabet = true;
			break;
		}
	}
	if (!symbolInAlphabet) {
		vectRangeOtherPile.clear();
		vectRangeDollarPile.clear();
		return 0;
	}
	vectRangeOtherPile.resize(1);
	vectRangeOtherPile[0].startPosN = 1;  //The first occurrence of symbol in F is in the first position in the pile Symbol
		//The last occurrence of the symbol prevSymbol in F is in the last position in the pile prevSymbol
		//It also corresponds to C[(unsigned int)(symbol) + 1]
		
	vectRangeOtherPile[0].endPosN = 0;
	
	for (dataTypedimAlpha mm = 0 ; mm < sizeAlpha; mm++){
		vectRangeOtherPile[0].endPosN += tableOcc[alpha[(unsigned int)(symbol)]][mm];
	}

	//symbPile stores the current pile
	symbPile = alpha[(unsigned int)(symbol)]; 
	
	#if DEBUG==1
		std::cerr << "backwardSearch - Initialization: Symbols of patterns in positions " << lenKmer << "\n";
		std::cerr << symbol  << "\t";		
		print(vectRangeOtherPile);
	#endif

	#if DEBUG==1
		dataTypeNChar maxInterval=0;
	#endif

	//std::cerr << "sizeof(rangeElement) " << sizeof(rangeElement) << "\n";
	#if DEBUG==1
	time_t start,end;
	#endif

	//START backwardsearch (prima iterazione con ultimo carattere fatta sopra)
	
	//vectRangeOtherPile è un vettore che memorizza gli intervalli della backwardsearch
	
	dataTypeNChar posSymb=lenKmer-1;

	for (; posSymb>0 && !vectRangeOtherPile.empty(); posSymb--) {   //For each symbol of the kmer

//		time(&start);

		#if DEBUG == 1
		std::cerr << "\nIteration: posSymb in Pattern = "  << (int) posSymb << "\n";
		#endif
		
		//Preceding symbol in the pattern
		symbol= kmer[posSymb-1];
		symbolInAlphabet = false;
		for (dataTypedimAlpha mm = 0; mm < sizeAlpha; mm++) {
			if (alphaInverse[mm] == symbol) {
				symbolInAlphabet = true;
				break;
			}
		}
		if (!symbolInAlphabet) {
			vectRangeOtherPile.clear();
			vectRangeDollarPile.clear();
			return lenKmer-posSymb;
		}

		#if DEBUG == 1
		std::cerr << "\n Iteration: posSymb in Pattern = "  << (int) posSymb << " symbol "<< symbol  << "\t";
		std::cerr << "symbPile= "<< (char)alphaInverse[(int)symbPile]  << "\n";
		#endif
		
		assert(link(rb_1,bsel_1)==1);

		#if DEBUG == 1
		cerr<<"backwardSearch - dopo link (including merge)"<<"\n";
		print_interval_number();
		cerr << "vectRangeDollarPile (computed): ";
		print(vectRangeDollarPile);		
		cerr << "-symbPile: " << (int)symbPile << "\n";
		cerr << "vectRangeOtherPile (original): ";
		print(vectRangeOtherPile);		
		#endif

		//For each symbol in the kmer we have to update both vectRangeDollarPile (if not empty) and vectRangeOtherPile 
		
		if (!vectRangeDollarPile.empty())
			assert ( updateIntervals (vectRangeDollarPile, fileInput, fileOutDecode, symbol, 0) == 1);

		assert ( updateIntervals (vectRangeOtherPile, fileInput, fileOutDecode, symbol, symbPile) == 1);
		
		//Append elements of vectRangeOtherPile to vectRangeDollarPile
		vectRangeDollarPile.insert(std::end(vectRangeDollarPile), std::begin(vectRangeOtherPile), std::end(vectRangeOtherPile));
		
		vectRangeDollarPile.shrink_to_fit();  //Requests the container to reduce its capacity to fit its size.
		vectRangeOtherPile.swap(vectRangeDollarPile);
		vector<rangeElement>().swap(vectRangeDollarPile);   // clear output reallocating 

		
		#if DEBUG==1
			if (maxInterval < vectRangeOtherPile.size())
				maxInterval = vectRangeOtherPile.size();
		#endif
		
		symbPile = alpha[(unsigned int)(symbol)]; 
		
		#if DEBUG == 1
			cerr << "After updateIntervals symbPile: " << (int)symbPile << "\n";
		#endif
		
		//MODIFICATO
		//#if RECOVERBW==0
		//	if (vectRangeOtherPile.empty()) {
		//		std::cerr << "Pattern not found, interrupted in position " << posSymb << "\n";			
		//	}
		//#endif

//		std::cerr << "--------------------------------------backward search: the cycle for " << "symbol in position " << (int)posSymb << " took " << difftime(end,start) << " seconds\n\n";
	
	}
	//END backwardsearch

	#if DEBUG==1
	cout<<"++++backwardSearch - maxInterval= "<< maxInterval << "\n";
	#endif
	
	
	#if DEBUG==1	
	print(vectRangeOtherPile);
	std::cerr << "++++++++++++++++++\n";
	fprintf(stderr, "##\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	
	//MODIFICATO -> tolto il for
	/*
		dataTypeNSeq intervalLength;
		for (dataTypeNSeq tt=0; tt<vectRangeOtherPile.size(); tt++) {
			intervalLength = vectRangeOtherPile[tt].endPosN - vectRangeOtherPile[tt].startPosN +1;
			if (intervalLength > 0){
				res = res + intervalLength;
			}
		}
		cout << "Found " << res << " starting positions of occurrences of pattern number " << n_kmer << endl;
	*/

	#if RECOVERBW == 1
		#if DEBUG == 1
			std::cerr << "\nbackwardSearch: compute starting point (findMultipleDollarsBackward)\n";
		#endif
		std::vector<rangeElementBW> vectRange;
		while(! vectRangeOtherPile.empty()){
			rangeElement eleR = vectRangeOtherPile.back();
			vectRangeOtherPile.pop_back();
			rangeElementBW eleC;
			eleC.pileN = symbPile;
			eleC.startPosN = eleR.startPosN;
			eleC.endPosN = eleR.endPosN;
			vectRange.insert(std::end(vectRange), eleC);
		}
		assert(findMultipleDollarsBackward(vectRange,fileInput,fileOutDecode,n_kmer,rb_1,bsel_1));
	//#else //COMMENTATO
	//	vector<rangeElement>().swap(vectRangeOtherPile);   // clear vectRange reallocating 
	#endif

	//MODIFICATO -> valore di ritorno lunghezza del match
	//return res;
	//posSymb -> contiene l'indice del mismatch, se vectRangeOtherPile è vuoto
	if (vectRangeOtherPile.empty())
		return lenKmer-posSymb-1;
	else
		return lenKmer;

}

#if RECOVERBW==1
int EDSBWT::updateSingleIntervalBW(std::vector<rangeElementBW> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, dataTypeNChar * counters, dataTypeNChar *numBlockCounter, dataTypeNChar * contInCurrentBlock, dataTypeNChar toRead, uchar *bufferBlock) {
			
				dataTypeNChar numBlock=0;
				uchar foundSymbol = '\0';  //here, it is not useful
				
				if (toRead > 0) {
					//we need to know how many occurrences of each symbol there are up to the position toRead.
					//if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
					//Before, we find the block where toRead position is.
					
					assert ( findBlockToRead(counters, currentPile, &toRead, &numBlock) == 1);
				}
				//else toTead==0 --> numBlock=0
		
				
				{
					//We are in another block and we need to read a new block in the BWT
					if(numBlock>0){
					for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
						counters[r] =  vectorOcc[currentPile][r][(numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
					}
					fseek (InFileBWT, numBlock*DIMBLOCK, 0);
			
					
					dataTypeNChar numberRead = rankManySymbols(*InFileBWT, counters, toRead, &foundSymbol, bufferBlock);
					assert (toRead <= numberRead); 
				}
				
				*contInCurrentBlock=toRead;
				*numBlockCounter = numBlock;
										
				#if DEBUG == 1				
					std::cerr << "\n updateSingleIntervalBW UPDATE: numBlockCounter " << *numBlockCounter << " contInCurrentBlock(toRead) " << *contInCurrentBlock << "\n";				
					std::cerr << "counters at the end of updateSingleIntervalBW:\t";
					for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
					   std::cerr << " " << counters[i];
					std::cerr << "\n";
				#endif
		
				return 1;
}
#endif

  int EDSBWT::computeCountersInSingleInterval(FILE *InFileBWT, dataTypeNChar toRead, dataTypedimAlpha currentPile,dataTypeNChar * counters) {
	//dataTypeNChar * counters = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
				
				dataTypeNChar numBlock = 0;
				#if DEBUG == 1
				uchar foundSymbol = '\0';  //here, it is not useful	
				#endif
				dataTypeNChar numberRead=0;
				
				//cont is the number of symbols already read!
				//toRead = vectTriple[k].posN - cont;				
				//for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
				//	counters[i]=0;
				
				#if DEBUG == 1
				std::cerr << "toRead of the vectRange[k].startPosN " << toRead << "\n";
				#endif
				if (toRead > 0) {
					//we need to know how many occurrences of each symbol there are up to the position toRead.
					//if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
					//Before, we find the block where toRead position is.
					assert ( findBlockToRead(counters, currentPile, &toRead, &numBlock) == 1);
				}

				if (toRead <= DIMBLOCK) {   //If toRead == DIMBLOCK, because I may need to known foundSymbol character
					//std::cerr << "Move file to the position " << numBlock*DIMBLOCK <<  "\n";
					fseek (InFileBWT, numBlock*DIMBLOCK, 0);
					#if DEBUG == 1
					std::cerr << "\t********foundSymbol " << (unsigned int)foundSymbol <<  "\n";
					#endif
					assert (toRead == numberRead);
					//cont += numberRead;
				}	
				
				return 1;
}

int EDSBWT::updateSingleInterval(std::vector<rangeElement> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, uchar symbol, dataTypeNChar * counters, dataTypeNChar *numBlockCounter, dataTypeNChar * contInCurrentBlock, dataTypeNChar toRead, uchar *bufferBlock) {
				dataTypeNChar numBlock=0;
				uchar foundSymbol = '\0';  //here, it is not useful
				
				#if DEBUG == 1
				std::cerr << "toRead of the vectRange[k].PosN " << toRead << "\n";
				#endif
				if (toRead > 0) {
					//we need to know how many occurrences of each symbol there are up to the position toRead.
					//if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks precede the block where the position toRead is.
					//Before, we find the block where toRead position is.
					
					numBlock = (dataTypeNChar)floor((long double)((toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
					
					#if DEBUG == 1
					 if (numBlock >= numBlocksInPartialBWT[currentPile])
						std::cerr << "     updateSingleInterval: numBlock " << numBlock << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";
					#endif
					assert(numBlock < numBlocksInPartialBWT[currentPile]);     //CONTROLLO MA SI PUO' TOGLIERE
					toRead = toRead - (numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK

				}
				//else toTead==0 --> numBlock=0
		
				if (numBlock == *numBlockCounter) {
					#if DEBUG == 1
						std::cerr << "Same numBlock: " << numBlock << " *contInCurrentBlock " << *contInCurrentBlock << " *numBlockCounter " << *numBlockCounter << "\n";
					#endif
					//We are in the same interval where counter is computed
					//counter is computed at the position contInCurrentBlock
					//In bufferBlock we already have the symbols
					//For each symbol in the buffer, it updates the number of occurrences into counters
					for (dataTypeNChar r=(*contInCurrentBlock); r<toRead; r++)        
						counters[alpha[(unsigned int)bufferBlock[r]]]++;    //increment the number of letter symbol into counters
				}
				else if (numBlock > *numBlockCounter) {
					//We are in another interval and we need to read a new block in the BWT
					for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
						counters[r] =  vectorOcc[currentPile][r][(numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
					
					//std::cerr << "Move file to the position " << numBlock*DIMBLOCK <<  "\n";
					fseek (InFileBWT, numBlock*DIMBLOCK, 0);
					
					//if (toRead > DIMBLOCK) {
					//	std::cerr << "(toRead > DIMBLOCK): toRead=" << toRead << " numBlock= "<<  numBlock << " *contInCurrentBlock " << *contInCurrentBlock << " *numBlockCounter " << *numBlockCounter << "\n";
					//}
					
					dataTypeNChar numberRead = rankManySymbols(*InFileBWT, counters, toRead, &foundSymbol, bufferBlock);
					#if DEBUG == 1
						std::cerr << "\t********foundSymbol " << (unsigned int)foundSymbol <<  "\n";
					#endif
					assert (toRead <= numberRead);  
				}
				
				*contInCurrentBlock=toRead;
				*numBlockCounter = numBlock;
		
								
				#if DEBUG == 1
					std::cerr << "counters  after vectRange[k]:\t";
					for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
					   std::cerr << " " << counters[i];
					std::cerr << "\n";
				#endif
		
		

				return 1;
}

void MergeAndRemove (std::vector< rangeElement > &vectRange, dataTypeNSeq &k,dataTypeNSeq &k_tmp){
		
	if (vectRange[k_tmp].startPosN > vectRange[k_tmp].endPosN)//k_tmp not a valid interval
	{
		#if DEBUG==1
		cerr << "\t REPLACE (" << vectRange[k_tmp].startPosN << " " <<  vectRange[k_tmp].endPosN << ") with"; 	
		cerr << " (" << vectRange[k].startPosN << " " <<  vectRange[k].endPosN << ")\n";
		#endif
		
		vectRange[k_tmp].startPosN = vectRange[k].startPosN;
		vectRange[k_tmp].endPosN = vectRange[k].endPosN;
	}
	else{ //k_tmp is a valid interval
		if ( (vectRange[k].startPosN == vectRange[k_tmp].endPosN +1) && (vectRange[k].startPosN <= vectRange[k].endPosN) ){
			#if DEBUG==1
			cerr << "\tMERGE (" << vectRange[k].startPosN << " " <<  vectRange[k].endPosN << ") in"; 	
			cerr << " (" << vectRange[k_tmp].startPosN << " " <<  vectRange[k_tmp].endPosN << ") -->"; 	
			#endif

			vectRange[k_tmp].endPosN = vectRange[k].endPosN;

			#if DEBUG==1
			cerr << "(" << vectRange[k_tmp].startPosN << " " <<  vectRange[k_tmp].endPosN << ")\n";
			#endif
		}		
		else if ( vectRange[k].startPosN > vectRange[k].endPosN)   
		{
			#if DEBUG==1
				cerr << "\tSKIP (" << vectRange[k].startPosN << " " <<  vectRange[k].endPosN << ") "; 	
				cerr << " k_tmp points to (" << vectRange[k_tmp].startPosN << " " <<  vectRange[k_tmp].endPosN << ")\n"; 
			#endif
		}		
		else {
			k_tmp++;
			#if DEBUG==1
			cerr << "\t COPY (" << vectRange[k].startPosN << " " <<  vectRange[k].endPosN << ") in"; 	
			cerr << " (" << vectRange[k_tmp].startPosN << " " <<  vectRange[k_tmp].endPosN << ")\n";
			#endif
			vectRange[k_tmp].startPosN = vectRange[k].startPosN;
			vectRange[k_tmp].endPosN = vectRange[k].endPosN;
		}
	}
}


int EDSBWT::updateIntervals(std::vector<rangeElement> &vectRange, string fileInput, string fileOutDecode, uchar symbol, dataTypedimAlpha currentPile){
	
	dataTypeNSeq nIntervals = vectRange.size();
	
	char *filenameIn = new char[128];
	char *filename = new char[9];
	const char *ext = ".aux";
	FILE *InFileBWT;

	#if DEBUG
	cout << "updateIntervals -- nIntervals = " << nIntervals << " currentPile= " << (int)currentPile << endl;
	#endif
	//Last = C[c] + rank (c, Last)				--> (vectTriple[1].pileN, vectTriple[1].posN)
	//First = C[c] + rank (c, First - 1) + 1    --> (vectTriple[0].pileN, vectTriple[0].posN)
	//So we write:
	
	for (dataTypeNSeq i=0; i<nIntervals; i++)    //For each kmer
		if (vectRange[i].endPosN >= vectRange[i].startPosN)   //if not, the kmer is not in the collection
			vectRange[i].startPosN --;   //So we compute rank until position First - 1


//2024-06-06
		dataTypeNChar * counters = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
		for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
			counters[i]=0;
		dataTypeNChar numBlockCounter = 0;   //number of the blocks read
		dataTypeNChar contInCurrentBlock = 0;    //number of the symbols read
		
		uchar *bufferBlock = new uchar[DIMBLOCK];
		uchar foundSymbol = '\0';  //here, it is not useful
		/////////////////////////

	dataTypeNSeq j = 0;
	dataTypeNSeq k_tmp=0;
	while (j < nIntervals) {   //nIntervals=vectRange.size()
		//The symbol for the sequences seqN in F[posN]  is the symbol
		//Now, I have to find the new symbol, it is in B[pileN] in position posN and so I can update pileN and posN
		
		
		sprintf (filename, "%s%d", "_bwt_",currentPile);
		sprintf (filenameIn,"%s%s%s",fileInput.c_str(), filename,ext);
		#if DEBUG==1
		std::cerr << "\nupdateIntervals: Current BWT-partial= " << (unsigned int)currentPile << "\n";
		#endif


		InFileBWT = fopen(filenameIn, "rb");
		if (InFileBWT==NULL) {
			std::cerr << "\nupdateIntervals: BWT file " << filenameIn << " j= " << (unsigned int)j << ": Error opening " << std::endl;
			exit (EXIT_FAILURE);
		}

	dataTypeNChar toRead = 0;
	dataTypeNChar numBlock=0;
		dataTypeNSeq k=j;
		#if DEBUG == 1
		cerr << "\t::::::Ora aggiornamento del PRIMO singolo intervallo nella pila " << (unsigned int)currentPile << ". Valgono j= " << j << " k= " <<  k << "\t";
		cerr << " vectRange[k].startPosN=" << vectRange[k].startPosN << " vectRange[k].endPosN " << vectRange[k].endPosN << "\n";
		#endif
		if (k< nIntervals) {
			if (vectRange[k].startPosN <= vectRange[k].endPosN) {
				
				//START
				toRead = vectRange[k].startPosN;   //startPosN has been decreased by 1, it can be 0
				if (toRead > 0) {
					numBlock = (dataTypeNChar)floor((long double)((toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
					//Update counters con the symbols up to the previous block
					if (numBlock > 0) {
						for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
							counters[r] =  vectorOcc[currentPile][r][(numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
					}					
				}
				//else (toRead == 0) --> numBlock=0
					
				#if DEBUG == 1
				//if (numBlock >= numBlocksInPartialBWT[currentPile])
					std::cerr << "\t First interval: toRead " << toRead << " numBlock " << numBlock << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";				    
				#endif
				assert(numBlock < numBlocksInPartialBWT[currentPile]);    
					
				toRead = toRead - (numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK
				
				//Move file pointer to the block numBlock and update counter at posizion toRead
				fseek (InFileBWT, numBlock*DIMBLOCK, 0);
				//Count the symbols in the selected block
				
				//if (toRead > DIMBLOCK) {
				//		std::cerr << "(toRead > DIMBLOCK): toRead=" << toRead << " numBlock= "<<  numBlock << " *contInCurrentBlock " << contInCurrentBlock << " *numBlockCounter " << numBlockCounter << "\n";
				//	}
				dataTypeNChar numberRead = rankManySymbols(*InFileBWT, counters, toRead, &foundSymbol, bufferBlock);
						#if DEBUG == 1
							std::cerr << "\t After first rank********foundSymbol " << (unsigned int)foundSymbol <<  "\n";
						#endif
				assert (toRead <= numberRead); //2024-06-04 	//== --> >=
				
				contInCurrentBlock=toRead;
				numBlockCounter = numBlock;		
				#if DEBUG == 1				
					std::cerr << "\t UPDATE: numBlockCounter " << numBlockCounter << " contInCurrentBlock " << contInCurrentBlock << "\n";
				#endif
				//assert (updateSingleInterval(vectRange, InFileBWT, k, currentPile, symbol, counters, &numBlockCounter, &contInCurrentBlock, vectRange[k].startPosN, bufferBlock) == 1);
				//Update the value in vectTriple[k].posN, it must contain the position of the symbol in F
				//Symbol is
				//newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
				//PosN is
				vectRange[k].startPosN = counters[alpha[(unsigned int)symbol]];
				for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
					vectRange[k].startPosN = vectRange[k].startPosN + tableOcc[g][alpha[(unsigned int)symbol]];
				}
				//First = c[c] + rank (c, First - 1) + 1
				vectRange[k].startPosN ++;  //We must to sum 1 to first
				
				//END
				assert (updateSingleInterval(vectRange, InFileBWT, k, currentPile, symbol, counters, &numBlockCounter, &contInCurrentBlock, vectRange[k].endPosN, bufferBlock) == 1);
				//I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
				//Symbol is
				//newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
				//PosN is
				vectRange[k].endPosN = counters[alpha[(unsigned int)symbol]];
				for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
					vectRange[k].endPosN = vectRange[k].endPosN + tableOcc[g][alpha[(unsigned int)symbol]];
				}
								
			}
		}
		
		//consecutive intervals and empty intervals
		if (k>0){
			MergeAndRemove(vectRange,k,k_tmp);
		}

		/////////////////


		k++;
		#if DEBUG == 1
			cerr << "\t::::ADESSO aggiornamento progressivo degli altri intervalli nella pila " << (unsigned int)currentPile << ". Valgono k_tmp= " << k_tmp << " k= " <<  k << "\n";
		#endif

		while (k< nIntervals) {
			if (vectRange[k].startPosN <= vectRange[k].endPosN) {
				//START
				assert (updateSingleInterval(vectRange, InFileBWT, k, currentPile, symbol, counters, &numBlockCounter, &contInCurrentBlock, vectRange[k].startPosN, bufferBlock) == 1);
				//I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
				//Symbol is
				//newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
				//PosN is
				vectRange[k].startPosN = counters[alpha[(unsigned int)symbol]];
				for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
					vectRange[k].startPosN = vectRange[k].startPosN + tableOcc[g][alpha[(unsigned int)symbol]];
				}
				//First = c[c] + rank (c, First - 1) + 1
				vectRange[k].startPosN ++;  //We must to sum 1 to first
				
				//END
				assert (updateSingleInterval(vectRange, InFileBWT, k, currentPile, symbol, counters, &numBlockCounter, &contInCurrentBlock, vectRange[k].endPosN, bufferBlock) == 1);
				//I have to update the value in vectTriple[k].posN, it must contain the position of the symbol in F
				//Symbol is
				//newSymb[vectTriple[k].seqN] = symbol;   //it is not useful here
				//PosN is
				vectRange[k].endPosN = counters[alpha[(unsigned int)symbol]];
				for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
					vectRange[k].endPosN = vectRange[k].endPosN + tableOcc[g][alpha[(unsigned int)symbol]];
				}
			}
			
//			if ((vectRange[k].seqN == vectRange[k_tmp].seqN) ) {//&& (vectRange[k].pileN == vectRange[k_tmp].pileN) )  {					
			MergeAndRemove(vectRange,k,k_tmp);
			//////////////////

			
			k++;
		}
		fclose(InFileBWT);

		j=k;
	}
	
	if ( (k_tmp==0) && (vectRange[k_tmp].startPosN > vectRange[k_tmp].endPosN)) {//in k_tmp not a valid interval
		#if DEBUG == 1
		cerr << "\tEND update (pattern non trovato): k_tmp= " << k_tmp << " start= " << vectRange[k_tmp].startPosN << " END: " <<  vectRange[k_tmp].endPosN  <<endl;
		#endif
		vectRange.resize(0);
	}
	
	if ( (vectRange[k_tmp].startPosN <= vectRange[k_tmp].endPosN)) {//in k_tmp is a valid interval
		vectRange.resize(k_tmp+1);
	}
	
	#if DEBUG == 1
	print_interval_number();
	print(vectRange);
	#endif
	
	delete [] bufferBlock;
	delete [] counters;
	
	delete [] filenameIn;
	delete [] filename;

	return 1;
}


dataTypeNChar EDSBWT::rankManySymbols(FILE & InFileBWT, dataTypeNChar *counters, dataTypeNChar toRead, uchar *foundSymbol, uchar *bufferBlock)
{
	dataTypeNChar numchar=0;
	//it reads toRead symbols from the fp file (Partial BWT)
	//while (toRead > 0) {            //((numchar!=0) && (toRead > 0)) {
		if (toRead <= DIMBLOCK) {    //Read toRead characters
			
			#if DEBUG==1
			std::cerr << "\t\t In rankManySymbols: toRead " << toRead << " is less than DIMBLOCK " << DIMBLOCK << "\n";
			#endif
			
			//fewer than x symbols could be read, but not less than toRead
			numchar = fread(bufferBlock,sizeof(uchar),DIMBLOCK,&InFileBWT);
			assert(numchar >= toRead);   //2024-06-04 	//== --> >=
			//The symbol of the sequence k.  It is the symbol in the last position in the partial BWT that we have read.
			if ((toRead > 0)) //2024-06-04
				*foundSymbol = bufferBlock[toRead-1];     
		}
		else {   //Read sizebuffer characters
			std::cerr << "rankManySymbols: Error to read is " << toRead << std::endl;
			exit (EXIT_FAILURE);
		}

		#if DEBUG==1
		std::cerr << "\t\t In rankManySymbols - counters before:\t";
		for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
			std::cerr << " " << counters[i];
		std::cerr << "\n";
		#endif

		//For each symbol in the buffer, it updates the number of occurrences into counters
		for (dataTypeNChar r=0; r<toRead; r++)
			counters[alpha[(unsigned int)bufferBlock[r]]]++;    //increment the number of letter symbol into counters

		#if DEBUG==1
		std::cerr << "\t\t In rankManySymbols - counters after: \t";
		for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
			std::cerr << " " << counters[i];
		std::cerr << "\n";
		std::cerr << "\t\t In rankManySymbols - bufferBlock after: \t";
		for (dataTypedimAlpha i = 0 ; i < numchar; i++)
			std::cerr << " " << bufferBlock[i];
		std::cerr << "\n";
		#endif

	//}
	

	return numchar;
}


int EDSBWT::link(rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1){

	// find the indexes of each dollar contained in each interval [startPosN,endPosN] of vectRangeOtherPile and store it in d.
	deque<dataTypeNSeq> d;
	dataTypeNSeq sizeVectRange = vectRangeOtherPile.size();

	dataTypeNChar start,end;

	for (dataTypeNSeq itVectRange=0; itVectRange < sizeVectRange; itVectRange++) {

		start = vectRangeOtherPile[itVectRange].startPosN;	
		end = vectRangeOtherPile[itVectRange].endPosN;

		if(start <= end){
			dollars_in_interval(symbPile,d,start,end,rb_1,bsel_1);
		}
	}

	//after obtaining the vector of the indexes (note that there are no repetitions), sort it.
	std::sort(d.begin(),d.end());


	//TO DO: limitation
	vectRangeDollarPile.reserve(d.size()*2/3);

	rangeElement current;


	while(!(d.empty())){

		current = preceding_dollars_finder(d.back(),rb_1,bsel_1);
		d.pop_back();

		if (!(vectRangeDollarPile.empty()) and (current.endPosN + 1 >= (vectRangeDollarPile.back()).startPosN)){
			if(current.startPosN < (vectRangeDollarPile.back()).startPosN){	
			
				(vectRangeDollarPile.back()).startPosN = current.startPosN;
			}
		}
		else {
			vectRangeDollarPile.insert(vectRangeDollarPile.end(),current);
		}

		dollars_in_interval(0,d,current.startPosN,current.endPosN,rb_1,bsel_1);
		//d.shrink_to_fit();

	}

	reverse(vectRangeDollarPile.begin(),vectRangeDollarPile.end());
	return 1;
}


void EDSBWT::dollars_in_interval(dataTypelenSeq symb, deque<dataTypeNSeq> &d_out,dataTypeNChar i,dataTypeNChar j,rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1){


	dataTypeNSeq l = EOFpos_rank[symb](i-1);	
	dataTypeNSeq u = EOFpos_rank[symb](j);	
	dataTypeNSeq index;

	for (dataTypeNSeq k = l; k < u; k++){
		index = EOF_ID[symb][k];

		if(index > first_symbol_index){	
			d_out.push_back(index);
			
		}
	}
	
}

// This function takes the position i of a dollar and the associated rank and select support (for 1) rb_1 and bsel_1
// of a compressed bitvector, computes the indexes of the dollars belonging to the starting and ending words
// of the segment on its left, and lastly returns the POSITIONS in the BWT of the characters
// which precede those dollars in the words. The positions are given as the beginning and end of the interval.
rangeElement EDSBWT::preceding_dollars_finder(dataTypeNSeq i, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1){
	dataTypeNChar a;


	a = rb_1(i+1); //a is the index of the segment containing the i-th dollar

	dataTypeNChar start = bsel_1(a-1); //bsel_1(k) gives the index of the k-th 1. That is, the index of the first word of the k-th segment.
	
	
	dataTypeNChar end = bsel_1(a)-1; //same as before. We subtract 1 because we want to find the index of the last word of the (a-1)-th segment.

	#if DEBUG == 1
	std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";
	#endif

	rangeElement output;
	output.startPosN = start+1;
	output.endPosN = end+1;
										
	return output;
}


void EDSBWT::indices_of_dollars_in_interval(std::ostream& searchOutput,dataTypelenSeq symbPile, dataTypeNChar i,dataTypeNChar j, dataTypeNSeq n_kmer, dataTypeNChar occPos, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1){


	dataTypeNSeq l = EOFpos_rank[symbPile](i-1);	
	dataTypeNSeq u = EOFpos_rank[symbPile](j);	

	dataTypeNSeq dollar_index_in_word;
	dataTypeNSeq symb_index;
	dataTypeNSeq dollar_index_in_symb;

	for (dataTypeNSeq k = l; k < u; k++){
		dollar_index_in_word = EOF_ID[symbPile][k];;
		symb_index = rb_1(dollar_index_in_word + 1);	
		dollar_index_in_symb = dollar_index_in_word - bsel_1(symb_index);

		//cout << "There is at least one occurrence of pattern number " << n_kmer << " starting at position " << occPos << " of word number " << dollar_index_in_word << ", which is word number " << dollar_index_in_symb << " inside of degenerate symbol number " << symb_index << endl;	

		searchOutput << n_kmer << "\t" << dollar_index_in_word << "\t" << symb_index << "\t" << dollar_index_in_symb << "\t" << occPos << endl;
	}



}


///////////////////////////////////////////////////////

    
int EDSBWT::computeVectorOcc(string filename)
{
	numBlocksInPartialBWT.resize(sizeAlpha);
	//Set number of blocks for each BWT-partial
	for (dataTypedimAlpha x = 0 ; x <= sizeAlpha-1; x++) {
		dataTypeNChar tot_size = EOFpos_rank[x].size();
		numBlocksInPartialBWT[x] = (dataTypeNChar)ceil((long double)tot_size/DIMBLOCK);
        //#if DEBUG==1
			std::cerr << "computeVectorUnbuildBCR: " << "freq[" << (unsigned int)x << "]= " << freq[alphaInverse[x]] << " and numBlocksInPartialBWT[" << (unsigned int)x << "]= " << numBlocksInPartialBWT[x] << "\n";
		//#endif
	}
	
	// Start by allocating an array for array of arrays
	vectorOcc.resize(sizeAlpha);    //For each BWT-partial
	time_t start,end;
	
	// alphaInverse[x] is the symbol to which correspond bwt_x
	
	//For each BWT-partial
	//Read BWT-partials in parallel
	static FILE *InFileBWT;
	dataTypedimAlpha x=0;
	
	for ( x=0; x < sizeAlpha; x++) {
		InFileBWT = openFilePartialIn(filename, x);
		fseek(InFileBWT, 0, SEEK_SET);
		
		// Allocate an array for each block of BWT-partial
		vectorOcc[x].resize(sizeAlpha);
		
		// Allocate an array of integers for each element of vectorOcc[x]
		for (dataTypedimAlpha y = 0 ; y < sizeAlpha; y++)   {      //For each block
			vectorOcc[x][y].resize(numBlocksInPartialBWT[x],0);
		}
		
		uchar *bufBlock = new uchar[DIMBLOCK];
		
		dataTypeNChar numBlock = 0;
		dataTypeNChar num_read = 0;
		
		time (&start);
		#if DEBUG==1
            std::cerr << "\ncomputeVectorUnbuildBCR: " << " start= " << start << "= " << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
        #endif

        //Read DIMBLOCK symbols in BWT-partial
		while( ( (num_read =  readOnFilePartial(bufBlock, DIMBLOCK, InFileBWT) ) && (num_read > 0) )  &&  (numBlock < numBlocksInPartialBWT[x]))   //Added check on numBlocks
        {
            for (dataTypeNChar i=0; i<num_read; i++) {
                vectorOcc[x][alpha[(unsigned int)(bufBlock[i])]][numBlock]++;
				//std::cerr << "---x = " << (unsigned int)x << " alpha " << (unsigned int)(bufBlock[i]) << " -numBlock " << numBlock << " vectorOcc is " << vectorOcc[x][alpha[(unsigned int)(bufBlock[i])]][numBlock] << ".\n";
            }
			numBlock++;
		}//end-while
			
		if ( !feof(InFileBWT) && (numBlock > numBlocksInPartialBWT[x])) {
			std::cerr << "computeVectorUnbuildBCR: Error - The file contains more blocks than allocates.\n" ;
            exit(1);
		}
                
        //For FIXED x
		//Compute the sum cumulative for each BWT-partial
		for (dataTypedimAlpha z = 0 ; z < sizeAlpha; z++)  {      //For each symbol z
			for(dataTypeNChar y = 1; y < numBlocksInPartialBWT[x] ; y++)   {      //For each block y>1 of partial-BWT x
				vectorOcc[x][z][y]=vectorOcc[x][z][y-1] + vectorOcc[x][z][y];   //Sum the previous one: ie Blcok y and block y-1
			}
		}
			
		//VectorOcc[x] stores info similar to row tableOcc[x] but keeping track of symbol occurrences into blocks
		//N.B. VectorOcc[x][z][y] stores the total number of z-occurrences in the BWT-partial corresponding to alpha[x] symbol up to the y-th block
		
		time (&end);

		#if DEBUG == 1
			double dif = difftime (end,start);
			std::cerr << "computeVectorUnbuildBCR: the cycle for " << "symbol = " << (int)x << " tooks " << dif << " seconds" << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
		#endif
			
		#if DEBUG==1
			std::cerr << "computeVectorUnbuildBCR: " << " end= " << end << "]= " << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
		#endif
			
		delete [] bufBlock;
		
		closeFilePartial(InFileBWT);
                
	}//end-for
        
        #if DEBUG==1
            for (dataTypedimAlpha x = 0 ; x < sizeAlpha; x++) {
                std::cerr << "x = " << (unsigned int)x << " For the " << alphaInverse[x] << "-BWT-partial: the #symbols is " << freq[alphaInverse[x]] << ".\n";
                std::cerr << "Number of blocks of the symbol " << numBlocksInPartialBWT[x] << "\n";
                for(dataTypedimAlpha z = 0; z < sizeAlpha; ++z) {
                    std::cerr << "Symbol " << (unsigned int)z << ":\t";
                    for(dataTypeNChar y = 0; y < numBlocksInPartialBWT[x]; ++y) {
                        std::cerr << (int)vectorOcc[x][z][y];
                    }
                    std::cerr << "\n";
                }
                std::cerr << "\n";
            }
        #endif
        
        return 1;
    }

int EDSBWT::findBlockToRead(dataTypeNChar *counters, dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock) {
        //Find the block numblock, where the position toRead is
    
        *numBlock = (dataTypeNChar)floor((long double)((*toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
        //if (*numBlock >= numBlocksInPartialBWT[currentPile])
        //std::cerr << "     findBlockToRead: numBlock " << *numBlock << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";
        assert(*numBlock < numBlocksInPartialBWT[currentPile]);
    
        if (*numBlock > 0) {
           for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
              counters[r] =  vectorOcc[currentPile][r][(*numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
           *toRead = *toRead - (*numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK
        }
        return 1;
    }
    
#if RECOVERBW==1
int EDSBWT::findBlockToReadBWT(dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock) {
        //Find the block numblock, where the position toRead is
    
        *numBlock = (dataTypeNChar)floor((long double)((*toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
        //if (*numBlock >= numBlocksInPartialBWT[currentPile])
        //std::cerr << "     findBlockToRead: numBlock " << *numBlock << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";
        assert(*numBlock < numBlocksInPartialBWT[currentPile]);       //2024-06-07 SUPERFLUO????
    
        if (*numBlock > 0) {           
           *toRead = *toRead - (*numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK
        }
        return 1;
    }
#endif
    
int EDSBWT::buildFreq(string fileWithExt) {
   
    //Open LCP and DA and BWT files
    string fnBWT = string(fileWithExt) + ".ebwt";
   
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
    uchar bwt;
    //uchar c;
    dataTypeNChar numcharBWT;
    numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT);
    
    freq[(unsigned int)(bwt)]=1;
    numEle++;
    
    //set freq
    dataTypeNSeq numSeq=0;
    numcharBWT=1;
    while ( ( numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT) ) && ( numcharBWT > 0 )  ) {
        
        if (bwt == TERMINATE_CHAR)
            numSeq++;

        freq[(unsigned int)(bwt)]++;
        
        numEle++;
        //cerr << numEle << " " << bwt << " " << (int)qs << "\n";
    }//end-while
    
	if(freq[TERMINATE_CHAR]==0){
		std::cerr << "ERROR: The end-marker must be " << TERMINATE_CHAR << endl;
		std::cerr << "If you want to use a different end-marker, set the parameter TERMINATE_CHAR in Parameters.h" << endl;
		exit(1);
	}
    //set nText
    nText = numSeq;
    
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
    std::cerr << "\tNumber of sequences: " << numSeq << "\n";
    std::cerr << "\tNumber of symbols in the input file: " << numEle << "\n";
    std::cerr << "\tSize alpha: " << (int)sizeAlpha << "\n";
    #if BUILD_LCP == 1
        std::cerr << "\tminLCP: " << minLCP << "\n";
        std::cerr << "\tmaxLCP: " << maxLCP << "\n";
    #endif

    fclose(InBWT);
    
    return 1;
}
    
	
int EDSBWT::recoverInfo(string filename) {         
    
	//Read FileInfo
    string fnInfoFile = filename + "_info" + ext ;
    FILE* InfoFile = fopen(fnInfoFile.c_str(), "rb");
    if (InfoFile==NULL) {
        std::cerr << "Error opening " << fnInfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
    }
	//Read BWT-length, num. sequences and sizeAlpha
    assert (fread(&lengthTot_plus_eof,sizeof(dataTypeNChar),1,InfoFile) ==1);
    assert (fread(&nText,sizeof(dataTypeNSeq),1,InfoFile) ==1);
    assert (fread(&sizeAlpha,sizeof(dataTypedimAlpha),1,InfoFile) ==1);
	//set alpha and alphaInverse
	alphaInverse = new dataTypedimAlpha[sizeAlpha];
    for (dataTypedimAlpha i = 0; i < sizeAlpha; ++i){
		assert (fread(&(alphaInverse[i]), sizeof(dataTypedimAlpha), 1, InfoFile)==1);
		alpha[alphaInverse[i]] = i;
	}
	
	std::cout << "\nFrom " << fnInfoFile << " file:\n";
    std::cout << "\tNumber of sequences: " << nText << "\n";
    std::cout << "\tTotal length (with $): " << lengthTot_plus_eof << "\n";
    std::cout << "\tSize alpha: " << (int)sizeAlpha << "\n";
	std::cout << "\tAlphabet: ";
	for (dataTypedimAlpha i = 0; i < sizeAlpha; ++i)
		std::cout << (char)alphaInverse[i] << "\t";
	std::cout << "\n";
	
	//Upload partial bitvectors with EOF positions
	EOFpos.resize(sizeAlpha);
	EOFpos_rank.resize(sizeAlpha);
	
	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
		char *fileBitVector = new char[256];
		sprintf (fileBitVector,"%s%s%u%s",filename.c_str(),"_bv_",j,ext.c_str());
		
		#if DEBUG==1
			cerr<<"Uploaded bitvector "<< fileBitVector << endl;
		#endif

		if (load_from_file(EOFpos[j],fileBitVector) != 1) {
			std::cerr << "Error loading bitvector from " << fileBitVector << std::endl;
			exit (1);
		}
		rrr_vector<>::rank_1_type tmp(&EOFpos[j]);
		EOFpos_rank.insert(EOFpos_rank.begin()+j,tmp);
	}
	#if DEBUG==1
		fprintf(stderr, "##AFTER uploading bitvectors with EOF pos\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	//for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
	//	cerr << EOFpos_rank[j].size() << endl;
	//}
	
	//Read dollars ID from FileInfo
	//set EOF_ID
    EOF_ID = new dataTypeNChar*[sizeAlpha];
	dataTypeNChar *numEOF = new dataTypeNChar[sizeAlpha];
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		numEOF[j] = EOFpos_rank[j](EOFpos_rank[j].size());
		//cout << numEOF[j] << endl;
        EOF_ID[j] = new dataTypeNChar[numEOF[j]];
    }
	
	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
		for (dataTypeNChar h = 0 ; h < numEOF[j]; h++) {
            assert (fread(&EOF_ID[j][h],sizeof(dataTypeNSeq),1,InfoFile) == 1);
        }
	}
	#if DEBUG==1
		fprintf(stderr, "##AFTER uploading EOF_ID\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	
    //Read tableOcc from FileInfo
	//set tableOcc
    tableOcc = new dataTypeNChar*[sizeAlpha];
    //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        //Update tableOcc with the symbols of the previous partial BWT files
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++) {
            assert(fread(&tableOcc[j][h],sizeof(dataTypeNChar),1,InfoFile) == 1);
        }
    }
    fclose(InfoFile);
	
    //#if DEBUG == 1
    std::cout << "\nFrom " << fnInfoFile << " file (TableOcc):\n";
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
            std::cout << tableOcc[j][h] << "\t";
		std::cout << std::endl;
    }
    //#endif

    return 1;
}


int EDSBWT::splitIntoPartial(string fileWithExt, int mode) {
    
    //Open BWT, LCP, DA, QS files
    string fnBWT = string(fileWithExt) + ".ebwt";
    vector < FILE * > InBWT;
    
    int t=0;
    
    #if OMP
    InBWT.resize(numthreads);
    for(t=0;t<numthreads; t++)
    #else
    InBWT.resize(1);
    #endif
    {
        InBWT[t] = fopen(fnBWT.c_str(), "rb");
        if (InBWT[t]==NULL) {
            std::cerr << "Error opening " << fnBWT << "!" << std::endl;
            exit (EXIT_FAILURE);
        }
    }
    
    if (mode == 1) {    //allocate tableOcc
        tableOcc = new dataTypeNChar*[sizeAlpha];
        //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
        for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
            tableOcc[j] = new dataTypeNChar[sizeAlpha];
        }
        for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
            for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
                tableOcc[j][h]=0;
    }
    

    vector < FILE *> OutFileBWT;
    #if OMP
        OutFileBWT.resize(numthreads);
    #else
        OutFileBWT.resize(1);
    #endif
        
    dataTypeNChar numEle=0;
    
    //sumFreq si to split ebwt
    vector<dataTypeNChar> sumFreq;
    sumFreq.resize(sizeAlpha+1);
    
    sumFreq[0]=0;
    for (dataTypedimAlpha j = 1 ; j <= sizeAlpha; j++)
        sumFreq[j]=sumFreq[j-1]+freq[alphaInverse[j-1]];
    
    size_t currentPile;
    #if OMP
#pragma omp parallel for default(shared) private(currentPile) firstprivate(t) num_threads(numthreads) schedule(dynamic, 1) reduction(+:numEle)
    #endif
        for (currentPile = 0 ; currentPile < sizeAlpha; ++currentPile) {
            
            assert(freq[alphaInverse[currentPile]] > 0);

            #if OMP
                t = omp_get_thread_num();//id_thread
                double start = omp_get_wtime();
            #endif
            
            //Open BCR partial files
            dataTypeNChar numcharBWT, numWrite;
            uchar *bufferBWT = new uchar[DIMBLOCK];
            
            string fnOutBWT = "bwt_" + to_string(currentPile)+ext;
            OutFileBWT[t] = fopen(fnOutBWT.c_str(), "wb");
            if (OutFileBWT[t]==NULL)
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening: " << fnOutBWT << std::endl;
                exit (EXIT_FAILURE);
            }
            fseek(InBWT[t], sumFreq[currentPile]*sizeof(uchar), SEEK_SET);
                      

            dataTypeNChar j = sumFreq[currentPile];
            
            dataTypeNChar toRead = DIMBLOCK;  //read DIMBLOCK symbols per time
            
            while (j < sumFreq[currentPile+1]) {
                
                if( (sumFreq[currentPile+1] - j) < toRead) //check remaining symbols to read
                    toRead = sumFreq[currentPile+1] - j;
                
                //read toRead symbols
                numcharBWT = fread(bufferBWT, sizeof(uchar),toRead, InBWT[t]);
                //write
                numWrite = fwrite(bufferBWT, sizeof(uchar), toRead, OutFileBWT[t]);
                assert(numcharBWT == numWrite);
                
                //set tableOcc
                if(mode == 1){
                    //counting the number of occurrences in BWT of the currentPile
                    uchar bwt;
                    for(dataTypeNChar i = 0; i < numcharBWT; i++ ){
                        bwt = bufferBWT[i];
                        tableOcc[(unsigned int)currentPile][alpha[(unsigned int)bwt]]++;
                    }
                }
    
                numEle++;
                
                j+=toRead;
                
            }  //end-for
            
            //Close partial files
            fclose(OutFileBWT[t]);
            
            #if OMP
            #pragma omp critical
            {
                std::cerr << "splitIntoPartial: THREAD = " << t << " tooks " << omp_get_wtime()-start << " seconds " << "\n\n";
            }
            #endif
            
        }  //end-for
    
    //Close BWT, LCP, DA, QS files
    #if OMP
    for(t=0;t<numthreads; t++)
    #endif
    {
        fclose(InBWT[t]);

    }
        
        
    std::cerr <<  "The total number of elements is " << numEle << "\n";
        
    return 1;
    
}







FILE * EDSBWT::openFilePartialIn(string filename, dataTypedimAlpha currentPile) {
    static FILE *inFile;
    string filenameIn = filename + "_bwt_" + to_string((int)(currentPile))+ext;
    inFile = fopen(filenameIn.c_str(), "rb");
    if (inFile==NULL) {
        #if OMP
        #pragma omp critical
        #endif
        {
            std::cerr << "openFilePartialIn: file currentPile=" << (unsigned int)currentPile << ": Error opening: " << filenameIn << std::endl;
            exit (EXIT_FAILURE);
        }
    }
    return inFile;
}

dataTypeNChar EDSBWT::readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT) {
   dataTypeNChar numchar;

   numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);

   return numchar;
}


int EDSBWT::closeFilePartial(FILE * pFile) {
    fclose(pFile);
    return 1;
}

void EDSBWT::print(std::vector<rangeElement> &vectRange){
			std::cerr << "size = " << vectRange.size() << std::endl;
			std::cerr << "startPosN and endPosN: ";
			for (dataTypeNSeq g = 0 ; g < vectRange.size(); g++) {
					std::cerr << vectRange[g].startPosN  << " " << vectRange[g].endPosN  << "\t";
			}
			std::cerr << std::endl;
}

#if RECOVERBW==1
void EDSBWT::printElementBW(std::vector<rangeElementBW> &vectRange){
			std::cerr << "size = " << vectRange.size() << std::endl;
			std::cerr << "pileN \t startPosN \t endPosN:\n";
			for (dataTypeNSeq g = 0 ; g < vectRange.size(); g++) {
					std::cerr << (int)vectRange[g].pileN << " \t " << vectRange[g].startPosN  << " \t " << vectRange[g].endPosN  << "\n";
			}
			std::cerr << std::endl;
}
#endif

void EDSBWT::print_interval_number (){
			std::cerr << "Number of intervals in vectRangeOtherPile: " << vectRangeOtherPile.size() << std::endl;			
}

EDSBWT::~EDSBWT()
{
	if (!resourcesLoaded) {
		return;
	}

	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		delete [] tableOcc[j];
		tableOcc[j] = NULL;
	}
	delete [] tableOcc;
	tableOcc = NULL;

	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		delete [] EOF_ID[j];
		EOF_ID[j] = NULL;
	}
	delete [] EOF_ID;
	EOF_ID = NULL;

	delete [] alphaInverse;
	alphaInverse = NULL;
	resourcesLoaded = false;
}
