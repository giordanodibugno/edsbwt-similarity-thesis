//
//  EDS-BWT
//
//  Created by Giovanna on 03/01/2021.
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

#include <iostream>
#include <assert.h>
#include <string.h>     // std::string, std::to_string
#include <sstream>
#include <stdio.h>
#include <math.h>

#include <map>


using std::cout;
using std::endl;

#include "EDSBWTsearch.hpp"
#include "Parameters.h"

using namespace std;

int main(int argc, char *argv[])
{
	//Input 
	if( argc != 3 ) {
		std::cerr << "usage: " << argv[0] << " inputEBWTfile inputEDSfile" << std::endl;
		std::cerr << "where:" << std::endl;
		std::cerr << "  inputEBWTfile is the BWT filename without the extension .ebwt (and .ebwt.qs for the QS string)" << std::endl;
		std::cerr << "  inputEDSfile is the EDS file whose loci will be searched" << std::endl;
		exit(1);
    }

	std::cout << "BCR_eds: " << argv[0] << std::endl;
	std::cout << "BCR_eds: The input ebwt file is " << argv[1] << std::endl;
	std::cout << "BCR_eds: The EDS file is " << argv[2] << std::endl;

	#if DEBUG == 1
	if (MODE == 1)
		std::cout << "unBCR_QS: MODE is 1 --> unBCR " << std::endl;
	else if (MODE == 2)
		std::cout << "unBCR_QS: MODE is 2 --> unBCR by using existing partial files" << std::endl;
	else
		std::cout << "unBCR_QS: MODE is 3 --> unBCR by using .table file" << std::endl;
	#endif
	
	string fileInput=argv[1];
	string filePattern=argv[2];

	int num_threads = 1;
 
	EDSBWT *BCRdec;

	BCRdec = new EDSBWT(fileInput, MODE, num_threads);
	float similarity = BCRdec->computeSimilarityFromEDS(filePattern);
	std::cout << "BCR_eds: Returned similarity is " << similarity << std::endl;
    
	#if RECOVERBW==1
		std::cerr << "\nThe csv file is ready! \n";
	#else
		std::cerr << "\nThe search is finished! \n";
	#endif

	delete BCRdec;

	std::cerr << "The End!\n";

    return 0;
}
