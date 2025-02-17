CC = g++ 

OMP = 0
FASTQ = 0

RECOVERBW = 1

DEFINES = -DRECOVERBW=$(RECOVERBW) -DFASTQ=$(FASTQ) -DOMP=$(OMP)

#OMP_LIB = -fopenmp

#SDSL
SDSL_INC = $(HOME)/local/include
SDSL_LIB = $(HOME)/local/lib

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES) 
#$(OMP_LIB)

#SDSL
CPPFLAGS += -I$(SDSL_INC) 
LDLIBS = -L$(SDSL_LIB) $(CPPFLAGS) -lsdsl -ldivsufsort -ldivsufsort64 -ldl

all: mainEDS-BWT converter-gsuf converter-BCR eds2fasta stringCheck

mainEDS-BWT_obs = mainEDS-BWT.o EDSBWTsearch.o Sorting.o malloc_count/malloc_count.o
mainEDS-BWT: $(mainEDS-BWT_obs)
	$(CC) -o EDSBWTsearch $(mainEDS-BWT_obs) $(LDLIBS)  

converter-gsuf_obs = da_to_everything.o
converter-gsuf: $(converter-gsuf_obs)
	$(CC) -o da_to_everything $(converter-gsuf_obs) $(LDLIBS)  
	
converter-BCR_obs = EOFpos_to_everything.o
converter-BCR: $(converter-BCR_obs)
	$(CC) -o EOFpos_to_everything $(converter-BCR_obs) $(LDLIBS)  

eds2fasta_obs = eds_to_fasta.o
eds2fasta: $(eds2fasta_obs)
	$(CC) -o eds_to_fasta $(eds2fasta_obs) $(LDLIBS)  

stringCheck_obs = stringCheck.o
stringCheck: $(stringCheck_obs)
	$(CC) -o stringCheck $(stringCheck_obs) $(LDLIBS)  

clean:
	rm -f core *.o *~ EDSBWTsearch da_to_everything EOFpos_to_everything eds_to_fasta stringCheck

depend:
	$(CC)  -MM *.cpp *.c > dependencies.mk

include dependencies.mk
