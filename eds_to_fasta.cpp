#include <iostream>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include "Parameters.h"

using namespace std;
using namespace sdsl;


int main (int argc, char *argv[]){
	if(argc != 3 ) {
		std::cerr << "usage: " << argv[0] << " input output" << std::endl;
		std::cerr << "input is the full filename, while output is the filename without any extension.\n";
		exit(1);
	}

	char current;
	char next;
	char null_char_input = EMPTY_CHAR_EDS;
	char null_char=EMPTY_CHAR;	//must be lexicographically greater than every other letter in the input file.
	bit_vector b;
	dataTypeNChar empty_counter=0;


	string inputName = argv[1];
	string outputName = argv[2];

	string edsName = inputName;
	ifstream edsFileForLength (edsName);
	cerr<<"Checking number of characters in "<<edsName<<" to initialize the bitvector\n";
	if (edsFileForLength.is_open()){
		edsFileForLength.seekg(0,std::ios_base::end);
		std::ios_base::streampos end_pos = edsFileForLength.tellg();
		b = bit_vector(end_pos,0);
		cerr<<edsName<<" contains "<<end_pos<<" characters\n";
	}
	else{
		cerr << "Error: could not open file "<<edsName<<endl;
		exit(1);
	}
	edsFileForLength.close();


	ifstream edsFile (edsName);

	string newName = outputName + ".fasta";
	ofstream lista (newName);

//	int k=1;	
	long int i=0;	
	long int chars=0;	

	cerr<<"Start reading "<<edsName<<endl;
	
	if (edsFile.is_open() && lista.is_open()){

		edsFile.get(current);

		if(current==EMPTY_CHAR){
			cerr << "\n\nError: the input contains character " << EMPTY_CHAR << ", which is used by this tool to represent an empty string. Please change EMPTY" << "_CHAR in Parameters.h or change the character in your input.\n\n\n";
			edsFile.close();
			lista.close();
			remove(newName);
			exit(1);
		}

		next=edsFile.peek();
		if(current=='{'){
			b[i]=1;
			i=i+1;
//			lista<<">"<<k<<endl;
			lista<<">1"<<endl;
			if (next == ','){
				lista<<null_char;
				empty_counter++;
				chars++;
			}
		}
		else {
			cout<<"Error: the string does not start with {."<<endl;
			return 0;
		}


		while (edsFile.get(current)){

			if(current==EMPTY_CHAR){
				cerr << "\n\nError: the input contains character " << EMPTY_CHAR << ", which is used by this tool to represent an empty string. Please change EMPTY" << "_CHAR in Parameters.h or change the character in your input.\n\n\n";
				edsFile.close();
				lista.close();
				remove(newName);
				exit(1);
			}

			if(current=='{'){
//				lista<<endl<<">"<<k<<endl;
				lista<<endl<<">1"<<endl;
				b[i]=1;
				i=i+1;
				next=edsFile.peek();
				if (next == ','){
					lista<<null_char;
					empty_counter++;
					chars++;
				}
				if(i%1000000 == 0){
					cerr<<i/1000000<<" million words read"<<"\r";
				}
			}
			else if(current == '}'){
//				k=k+1;
//				if(max_l < l_tmp){
//					max_l = l_tmp;
//				}
//				l_tmp = 0;
			}
			else if(current == ','){
				b[i]=0;
				i=i+1;
				lista<<endl<<">2"<<endl;
				next=edsFile.peek();
				if (next == ',' or next == '}'){
					lista<<null_char;
					empty_counter++;
					chars++;
				}
			}
			else if(current == null_char_input){
				lista<<null_char;
				empty_counter++;
				chars++;
			}
			else {
				lista<<current;
//				l_tmp++;
				chars++;
			}
		}
		lista<<endl;
		edsFile.close();
		lista.close();
	}
	else {
		cerr << "Error: could not open file " << edsName << " or "<< newName << endl;
		return 0;
	}

	cerr<<"\nFinished reading"<<endl;


	b.resize(i);
	string bvName = outputName + ".bitvector";
	store_to_file(b,bvName);


	dataTypeNChar bwt_length = chars + i;
	string empty_name = outputName+".empty.info";
	FILE* empty_number_file = fopen(empty_name.c_str(),"wb");
	if (empty_number_file != NULL){
		fwrite(&bwt_length, sizeof(dataTypeNChar), 1, empty_number_file);
		fwrite(&empty_counter, sizeof(dataTypeNChar), 1, empty_number_file);
	fclose(empty_number_file);
	}
	else{
		cerr << "Error: could not open file "<<empty_name<<endl;
		return 0;
	}

	
	return 1;
}
