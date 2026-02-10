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
 
 /*
 * Setting
 */

 
#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <fstream>

#define DEBUG 0

#define SIZEBUFFER 1024     //Size of the buffer for I/O partial ebwt/LCP/DA/SA

#define TERMINATE_CHAR '$'     //it is the symbol used as "end of strings", it must be lexicographically smaller than all the letters of the alphabet
#define TERMINATE_CHAR_LEN uchar(255)      //it is stored in cyc files, it is ignored by the algorithm, so it must not belong to the alphabet

#define EMPTY_CHAR_EDS 'E'     //it is the symbol (if any) used to represent the empty word in the input eds. It will be converted in NULL_CHAR during the preprocessing 
#define EMPTY_CHAR 'Z'     //it is the symbol intermally used as empty word. It must not appear in the input eds, and must be lexicographically greater than all the characters of the alphabet

#define KEEPEBWT 1     //set this to 0 to automatically delete the .ebwt file

#define SIZE_ALPHA 256  

#define sizeKmer 100
#define BACKWARD 1
#define DECODE 0

#define MODE 1
//MODE can be set equal to
// 1 --> for decoding eBWT by using the maximum read length only (maxLengthRead)
// 2 --> for decoding eBWT by using existing partial ebwt files, in addition to .info and .table files
// 3 --> for decoding eBWT by using only .info and .table files


typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))

/* USE: 0 - below 255 (unsigned char)
 *		1 - between 256 and 65.536 (unsigned short)
 *		2 - between 65.536 and 4.294.967.296 (unsigned int)
 *		3 - otherwise (unsigned long)
 */

// Type size for Sequences Length (in biologic case 100)
#define dataTypeLengthSequences 1		

// Type size for Number of Sequences in the input file
#define dataTypeNumSeq 2		

// Type size for Number of Character in the input file (length of the BWT)
#define dataTypeNumChar 3		


//Set the types
#if dataTypeLengthSequences == 0
	#define dataTypelenSeq uchar
#elif dataTypeLengthSequences == 1
	#define dataTypelenSeq ushort
#elif dataTypeLengthSequences == 2
	#define dataTypelenSeq uint
#elif dataTypeLengthSequences == 3
	#define dataTypelenSeq ulong	
#endif

#if dataTypeNumSeq == 0
	#define dataTypeNSeq uchar
#elif dataTypeNumSeq == 1
	#define dataTypeNSeq ushort
#elif dataTypeNumSeq == 2
	#define dataTypeNSeq uint
#elif dataTypeNumSeq == 3
	#define dataTypeNSeq ulong	
#endif


#if dataTypeNumChar == 0
	#define dataTypeNChar uchar
#elif dataTypeNumChar == 1
	#define dataTypeNChar ushort
#elif dataTypeNumChar == 2
	#define dataTypeNChar uint
#elif dataTypeNumChar == 3
	#define dataTypeNChar ulong	
#endif



//if you want to delete the partial files, please set it to 1
#define deletePreprocessingFile 0



//if KEEP_eBWT_IN_EXT_MEMORY==1, BCR uses files for partials ebwts
//if KEEP_eBWT_IN_EXT_MEMORY==0, BCR uses strings for partials ebwts
//In both cases, SA, DA, LCP are stored in files.
#define KEEP_eBWT_IN_EXT_MEMORY  1



#endif
