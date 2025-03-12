# EDS-BWT
This tool computes the Burrows-Wheeler transform of an elastic-degenerate string, called EDS-BWT, and can search patterns on it.

An *elastic degenerate string* (eds) is the concatenation of a collection of strings, as in the following example:

```sh
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

The collections of strings between curly brackets are called *degenerate symbols*, or *segments*.



## Install

```sh
git clone --recursive https://github.com/giovannarosone/EDS-BWT.git 
cd EDS-BWT
```

### Requirements

After dowloading EDS-BWT, install [BCR](https://github.com/giovannarosone/BCR_LCP_GSA)

```sh
git clone https://github.com/giovannarosone/BCR_LCP_GSA
cd BCR_LCP_GSA
make STORE_INDICES_DOLLARS=1
cd ..
```

or alternatively, [gsufsort](https://github.com/felipelouza/gsufsort)

```sh
git clone https://github.com/felipelouza/gsufsort.git
cd gsufsort
make TERMINATOR=0 DNA=1
cd ..
```

You can specify the path where BCR (or gsufsort) is installed in variable BCRPATH (or GSUFPATH) of file EDS-BWTransform.sh.
You also need to install the [sdsl library](https://github.com/simongog/sdsl-lite) and specify the path of SDSL/include and SDSL/lib in the parameters SDSL_INC and SDSL_LIB of file makefile.


### Compile

```sh
make
```

If you want the tool to just search for the number of occurrences of a pattern, without finding the positions, then compile with option RECOVERBW=0, as follows
```sh
make RECOVERBW=0
```





## Run

To compute the EDS-BWT of an elastic degenerate string written in file input.eds, please run:

```sh
./EDS-BWTransform.sh input output
```
where output is used as the base name of the output files. By default, the script `EDS-BWTransform.sh` computes the extended Burrow-Wheeler transform of the collection of string in the eds by using BCR. However, one can choose to select gsufsort by setting BCR=0 in `EDS-BWTransform.sh`.  

To search one or more patterns, contained in file patterns, please run:
```sh
./EDSBWTsearch output patterns
```


## Input files

### Elastic degenerate string
The input file must be a text file with extension .eds, where each degenerate symbol is contained in curly bracket and each string in a degenerate symbol is separated by comma.
Empty strings can be represented both by empty strings themselves, or by using the character for empty string EMPTY_CHAR_EDS (default is E, can be changed in Parameters.h).
Example:

```sh
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,}{CTGT}
```
and
```sh
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,E}{CTGT}
```
represent the same eds.

If your input file input.txt does not have curly brackets around all degenerate symbols, such as the following example
```sh
ATTGCT{CTA,TA,A}CTACGGACT{A,}CTGT
```
then please execute
```sh
./stringCheck input.txt output
```

This produces the file output.eds, which will have the brackets as required.


Finally, your input must not contain the character EMPTY_CHAR (default is Z), which is internally used to represent empty strings.
If your input does contain this character, please change the parameter EMPTY_CHAR in Parameters.h to a character that is not used in your input and is lexicographically greater than all the characters in your input. 

### Patterns

The patterns file is a text file where each pattern is in a different line, ending with UNIX end-of-line markers (LF), as opposed to Windows ones (CR+LF).
Empty lines are not allowed. Example:

```sh
TATT
ACT
TTAT
```


## Output

### Transform output

Output files are written in the same directory as the input eds file. The filenames are the same as the input, with different extensions, unless a filename is specified when running. The output is the following:

- input.bitvector: an sdsl bitvector corresponding to the bitvector associated to the eds;

- input.ebwt: the extended Burrow-Wheeler transform of the collection of string in the eds. This is not used for the search, and can be removed by setting the parameter KEEPEBWT to 0 in Parameters.h;

- for each character in the eds, a file input_bwt_i.aux is created, where i goes from 0 to the number of characters in the eds. It contains the part of the ebwt corresponding to the i-th character of the alphabet (the 0-th being the end-marker symbol, which is TERMINATE_CHAR in Parameters.h). Characters are encoded in binary and take 4 bytes each;

- for each character in the eds, a file input_bv_i.aux is created, where i goes from 0 to the number of characters in the eds. It contains an sdsl bitvector which have a 1 in the positions where the ebwt is the end-marker symbol;

- input_info.aux: a binary file containing auxiliary informations for the search, such as the size of the alphabet, the frequence table of the characters in the ebwt and an array which contains the indexes of each end-marker symbol


### Search output

EDS-BWT outputs the starting position for each occurrence of the pattern, giving the degenerate symbol, the string number and the position in the string. If more than one occurrence begins at the same position in the same string, they are counted as just one occurrence.

The ouput is a .csv file, written in the same directory as the pattern, with "output.csv" appended to the filename.
Each line is an occurrence. Columns are separated by tabulation, not commas, and give the following informations, in order:

- \#Pat: the index of the pattern found;

- $_i: the index i of the string (starting from 0) at which the occurrence begins;

- D\[i\]: the index i (starting from 1) of the degenerate symbol at which the occurrence begins, the eds being D\[1\]...D\[k\];

- S_j: the index j (starting from 0) of the string, contained in degenerate symbol D\[i\] = {S_0,...S_q}, at which the occurrence begins;

- S_j\[r\]: the position (starting from 0) at which the occurrence begins inside the string.


For example, the output for searching patterns
```sh
TATT
ACT
TTAT
```

in
```sh
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

is

```sh
\#Pat	$_i	D[i]	S_j	S_j[r]
1	3	2	0	0
1	1	1	1	0
1	4	3	0	1
1	7	4	0	1
3	4	3	0	0
3	7	4	0	0
3	0	1	0	1
```

It can be checked using the tool with the following commands:
```sh
./EDS-BWTransform.sh sample/test/test
./EDSBWTsearch sample/test/test sample/test/kmers.txt
```

## References

Cioni L., Guerrini V., Rosone G., [The Burrows-Wheeler Transform of an Elastic-Degenerate String](https://ceur-ws.org/Vol-3811/paper240.pdf). In *Proceedings of the 25th Italian Conference on Theoretical Computer Science* ICTCS 2024. CEUR-WS.org, vol 3811, pp 66--80.

--

<small> Supported by PNRR project “THE—Tuscany Health Ecosystem” — Spoke 6 “Precision medicine & personalized healthcare”, funded by the European Commission under the NextGeneration EU programme.

