# Modifiche a EDSBWTsearch per Similarita tra EDS

Questo documento riassume le modifiche fatte durante la sessione di lavoro su `EDSBWTsearch`, con l'obiettivo di avvicinare il programma al calcolo di una misura di similarita tra pangenomi rappresentati come Elastic Degenerate Strings (EDS).

## Obiettivo Iniziale

Il programma originale `EDSBWTsearch` prendeva in input:

```sh
./EDSBWTsearch indice file_pattern.txt
```

dove:

- `indice` e il basename di un indice EDS-BWT gia costruito;
- `file_pattern.txt` contiene una stringa/pattern per riga.

Il programma eseguiva una backward search per ogni pattern e, nella versione modificata inizialmente, calcolava un singolo valore globale:

```text
somma delle lunghezze matchate / somma delle lunghezze dei pattern
```

Per la tesi, invece, l'obiettivo e lavorare su EDS complete, organizzate in loci:

```text
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

e calcolare una similarita per ogni locus.

## Backup Creato

Prima della rifattorizzazione principale e stata salvata una copia della versione funzionante con similarita per locus:

```text
EDSBWTsearch.cpp.backup-locus-similarity
```

Questo file conserva la versione in cui:

- il programma leggeva gia un file `.eds`;
- calcolava le similarita per locus;
- ma il calcolo era ancora eseguito direttamente nel costruttore `EDSBWT::EDSBWT(...)`.

## Correzioni Preliminari

Nel ciclo originale sui k-mer c'erano due errori/parti incomplete.

Prima:

```cpp
den += kmer;
std::cout << "Similarita test con kmers ("<<num<<"/"<<den<<"):" << percentage << endl;
```

Dopo:

```cpp
den += lenKmer;
std::cout << "Similarita test con kmers ("<<num<<"/"<<den<<"):" << sim << endl;
```

Motivo:

- `den` e un numero, quindi non puo sommare una `std::string`;
- `percentage` non era dichiarata;
- `sim` era gia calcolata come `num / den`.

Questa e stata la prima versione della misura globale.

## Cambio Input: Da File TXT a File EDS

Il secondo input di `EDSBWTsearch` non e piu un file `.txt` con un pattern per riga.

Ora il programma si usa cosi:

```sh
./EDSBWTsearch indice query.eds
```

Esempio:

```sh
./EDSBWTsearch test sample/test/test.eds
```

Il primo argomento rimane il basename dell'indice EDS-BWT:

```text
test
```

Il programma continua quindi a usare file come:

```text
test_info.aux
test.bitvector
test_bwt_0.aux
test_bv_0.aux
...
```

Il secondo argomento e ora una EDS completa, da cui vengono estratti i loci.

## Parser EDS Aggiunto

In `EDSBWTsearch.cpp` e stata aggiunta la funzione:

```cpp
static vector<vector<string> > readEDSByLocus(string filename)
```

Questa funzione legge un file `.eds` e restituisce una struttura:

```cpp
vector<vector<string> > loci;
```

Esempio:

```text
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

diventa:

```cpp
loci[0] = {"TT", "TA", "C"};
loci[1] = {"T"};
loci[2] = {"TT", "A", "G"};
loci[3] = {"TT", ""};
loci[4] = {"ATT"};
```

Il parser:

- riconosce `{` come inizio locus;
- riconosce `}` come fine locus;
- riconosce `,` come separatore tra stringhe alternative;
- ignora spazi, tab e newline;
- gestisce `EMPTY_CHAR_EDS` come stringa vuota;
- segnala errori se trova parentesi annidate, virgole fuori da un locus, simboli fuori da un locus o una parentesi non chiusa.

## Similarita Per Locus

Il vecchio ciclo sui pattern e stato sostituito da un ciclo annidato:

```cpp
for each locus:
    for each string in locus:
        backwardSearch(...)
```

Per ogni stringa non vuota del locus viene calcolata:

```text
len_match / len_string
```

Nel codice, il valore del locus viene calcolato in forma equivalente:

```text
locusSim = somma len_match / somma len_string
```

Le stringhe vuote vengono saltate:

```cpp
if (lenKmer == 0) {
    cout << "len_match=0 (empty string skipped)" << endl << endl;
    continue;
}
```

Quindi una stringa vuota non contribuisce ne al numeratore ne al denominatore.

Per ogni locus il programma stampa:

```text
Similarita locus i (num/den):score
```

Alla fine stampa anche:

```text
Array similarita loci: [...]
Similarita totale EDS (num/den):score
```

## Modifica Di backwardSearch

La funzione:

```cpp
int EDSBWT::backwardSearch(...)
```

prima ritornava sostanzialmente un indicatore/numero legato alla presenza del pattern.

E stata modificata per ritornare la lunghezza matchata:

```cpp
if (vectRangeOtherPile.empty())
    return lenKmer-posSymb-1;
else
    return lenKmer;
```

Quindi:

- se il pattern viene matchato completamente, ritorna `lenKmer`;
- se la ricerca si interrompe, ritorna la lunghezza del match trovato prima del mismatch.

Questo valore e usato direttamente nel calcolo della similarita.

## Rifattorizzazione Del Costruttore

Prima, il costruttore:

```cpp
EDSBWT::EDSBWT(string fileInput, string filepatterns, int mode, int num_threads)
```

faceva troppe cose:

1. caricava l'indice;
2. caricava i bitvector;
3. leggeva l'input di query;
4. eseguiva tutte le backward search;
5. calcolava le similarita;
6. stampava il risultato;
7. liberava memoria.

Questo rendeva difficile riusare il codice per il futuro algoritmo:

```text
A(P1,P2)
A(P2,P1)
media finale
```

Per questo e stato aggiunto un nuovo costruttore:

```cpp
EDSBWT(string fileInput, int mode, int num_threads);
```

Questo costruttore carica solo l'indice EDS-BWT.

Il vecchio costruttore e stato mantenuto per compatibilita:

```cpp
EDSBWT(string fileInput, string fileEDS, int mode, int num_threads);
```

ma ora internamente chiama il nuovo costruttore e poi chiama la funzione di similarita:

```cpp
EDSBWT::EDSBWT(string fileInput, string fileEDS, int mode, int num_threads)
    : EDSBWT(fileInput, mode, num_threads)
{
    computeSimilarityFromEDS(fileEDS);
}
```

## Nuovo Metodo computeSimilarityFromEDS

In `EDSBWTsearch.hpp` e stato dichiarato:

```cpp
float computeSimilarityFromEDS(string fileEDS);
```

In `EDSBWTsearch.cpp` e stato implementato:

```cpp
float EDSBWT::computeSimilarityFromEDS(string fileEDS)
```

Questa funzione:

1. crea i supporti `rank` e `select` sul bitvector dell'indice;
2. se `RECOVERBW=1`, prepara il CSV delle posizioni;
3. legge la EDS con `readEDSByLocus`;
4. calcola la similarita di ogni locus;
5. stampa l'array delle similarita;
6. calcola la similarita totale della EDS rispetto all'indice;
7. restituisce il valore totale come `float`.

Questa e la modifica piu importante per il prossimo passo, perche ora il codice puo fare:

```cpp
EDSBWT indexP2("P2", MODE, 1);
float a12 = indexP2.computeSimilarityFromEDS("P1.eds");
```

e quindi potra calcolare:

```cpp
EDSBWT indexP1("P1", MODE, 1);
EDSBWT indexP2("P2", MODE, 1);

float a12 = indexP2.computeSimilarityFromEDS("P1.eds");
float a21 = indexP1.computeSimilarityFromEDS("P2.eds");

float sim = (a12 + a21) / 2.0;
```

## Campi Aggiunti Alla Classe EDSBWT

In `EDSBWTsearch.hpp` sono stati aggiunti:

```cpp
string fileInputBase;
bool resourcesLoaded;
```

`fileInputBase` conserva il basename dell'indice caricato, per esempio:

```text
test
```

Questo permette a `computeSimilarityFromEDS(...)` di chiamare:

```cpp
backwardSearch(fileInputBase.c_str(), ...)
```

senza dover ricevere di nuovo il basename dell'indice.

`resourcesLoaded` serve al distruttore per sapere se deve liberare le risorse allocate.

## Gestione Memoria Spostata Nel Distruttore

Prima, la memoria veniva liberata alla fine del costruttore.

Questo non era piu corretto dopo la rifattorizzazione, perche l'oggetto deve restare vivo e riusabile dopo il caricamento dell'indice.

La liberazione di:

```cpp
tableOcc
EOF_ID
alphaInverse
```

e stata spostata in:

```cpp
EDSBWT::~EDSBWT()
```

Il distruttore ora controlla:

```cpp
if (!resourcesLoaded) {
    return;
}
```

e poi libera le strutture allocate.

## Modifica Del Main

In `mainEDS-BWT.cpp` il messaggio di usage e stato aggiornato.

Prima:

```text
usage: ./EDSBWTsearch inputEBWTfile inputPATTERNfile
```

Ora:

```text
usage: ./EDSBWTsearch inputEBWTfile inputEDSfile
```

Il main ora fa:

```cpp
BCRdec = new EDSBWT(fileInput, MODE, num_threads);
float similarity = BCRdec->computeSimilarityFromEDS(filePattern);
std::cout << "BCR_eds: Returned similarity is " << similarity << std::endl;
```

Quindi il valore non e piu solo stampato dentro la classe, ma viene restituito al chiamante.

Inoltre il `main` ora termina con:

```cpp
return 0;
```

Prima tornava `1`, quindi la shell vedeva il programma come fallito anche quando l'esecuzione era corretta.

## README Aggiornato

Il `README.md` e stato aggiornato per descrivere il nuovo comportamento:

```sh
./EDSBWTsearch output input.eds
```

La sezione `Patterns` e stata sostituita da `Search EDS`.

La documentazione ora dice che:

- ogni degenerate symbol della EDS viene trattato come un locus;
- ogni stringa nel locus viene cercata come pattern;
- viene prodotto un valore di similarita per locus;
- le stringhe vuote vengono ignorate nel punteggio normalizzato;
- se `RECOVERBW=1`, continua a esistere anche il CSV delle posizioni.

## Comandi Di Compilazione Usati

La compilazione e stata fatta con `RECOVERBW=0`, usando la SDSL locale del repository:

```sh
make -B mainEDS-BWT RECOVERBW=0 SDSL_INC=/home/giordanodb/EDS-BWT/sdsl_lib/include SDSL_LIB=/home/giordanodb/EDS-BWT/sdsl_lib/lib
```

`-B` forza la ricompilazione completa degli oggetti con la macro corretta.

## Test Eseguiti

### Test Con EDS Identica All'Indice

Comando:

```sh
./EDSBWTsearch test sample/test/test.eds
```

Risultato atteso: similarita completa.

Output rilevante:

```text
Array similarita loci: [1, 1, 1, 1, 1]
Similarita totale EDS (15/15):1
BCR_eds: Returned similarity is 1
```

### Test Con EDS Di Query Non Perfetta

E stata creata una query temporanea:

```text
/tmp/edsbwt_noncomplete_query.eds
```

contenente:

```text
{TT,GG}{AAAA}{C}{TATA,}{ATT}
```

Comando:

```sh
./EDSBWTsearch test /tmp/edsbwt_noncomplete_query.eds
```

Output rilevante:

```text
Array similarita loci: [0.75, 0.5, 1, 1, 1]
Similarita totale EDS (13/16):0.8125
BCR_eds: Returned similarity is 0.8125
```

Questo conferma che non viene prodotto sempre un match completo e che il calcolo per locus risponde ai match parziali.

### Simulazione Dell'Algoritmo Completo Su Una EDS Con Se Stessa

E stata fatta una prova concettuale dell'algoritmo completo usando la stessa EDS sia come `P1` sia come `P2`.

```text
P1 = sample/test/test.eds
P2 = sample/test/test.eds
```

Poiche l'indice `test` e gia costruito sulla stessa EDS, il confronto completo e stato simulato eseguendo due volte la ricerca:

```sh
./EDSBWTsearch test sample/test/test.eds
./EDSBWTsearch test sample/test/test.eds
```

Questo corrisponde a:

```text
A(P1,P2) = cerca P1.eds nell'indice di P2
A(P2,P1) = cerca P2.eds nell'indice di P1
```

Nel caso particolare `P1 = P2`, entrambe le direzioni usano lo stesso indice e la stessa query.

Output rilevante in entrambe le esecuzioni:

```text
Array similarita loci: [1, 1, 1, 1, 1]
Similarita totale EDS (15/15):1
BCR_eds: Returned similarity is 1
```

Quindi:

```text
A(P1,P2) = 1
A(P2,P1) = 1
```

La media finale dell'algoritmo completo sarebbe:

```text
(1 + 1) / 2 = 1
```

Questo conferma che, nel caso di una EDS confrontata con se stessa, il comportamento attuale produce la massima similarita attesa.

### Esempi EDS Disponibili Nel Repository

Nel repository sono presenti due file `.eds` di esempio:

```text
sample/test/test.eds
sample/ExamplePaper/examplePaper.eds
```

Il contenuto di `sample/test/test.eds` e:

```text
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

Il contenuto di `sample/ExamplePaper/examplePaper.eds` e:

```text
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,}{CTGT}
```

Quindi esiste gia un secondo esempio per fare un confronto reale tra due EDS diverse. Per usarlo nell'algoritmo completo, pero, serve costruire anche il suo indice EDS-BWT, ad esempio con:

```sh
./EDS-BWTransform.sh sample/ExamplePaper/examplePaper sample/ExamplePaper/examplePaper
```

Poi si potranno calcolare le due direzioni:

```sh
./EDSBWTsearch sample/ExamplePaper/examplePaper sample/test/test.eds
./EDSBWTsearch test sample/ExamplePaper/examplePaper.eds
```

e infine mediare i due valori restituiti.

### Struttura Delle Due EDS Usate Nel Confronto

Le due EDS considerate per il confronto reale sono:

```text
P1 = sample/test/test.eds
P2 = sample/ExamplePaper/examplePaper.eds
```

#### P1: sample/test/test.eds

Contenuto:

```text
{TT,TA,C}{T}{TT,A,G}{TT,}{ATT}
```

Struttura per locus:

```text
1: {TT,TA,C}       -> 3 stringhe, lunghezze 2,2,1
2: {T}             -> 1 stringa, lunghezza 1
3: {TT,A,G}        -> 3 stringhe, lunghezze 2,1,1
4: {TT,}           -> 2 stringhe, una vuota, lunghezza non vuota 2
5: {ATT}           -> 1 stringa, lunghezza 3
```

Totali:

```text
numero loci: 5
numero stringhe/alternative: 10
stringhe non vuote: 9
stringhe vuote: 1
caratteri totali biologici non vuoti: 15
```

#### P2: sample/ExamplePaper/examplePaper.eds

Contenuto:

```text
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,}{CTGT}
```

Struttura per locus:

```text
1: {ATTGCT}        -> 1 stringa, lunghezza 6
2: {CTA,TA,A}      -> 3 stringhe, lunghezze 3,2,1
3: {CTACGGACT}     -> 1 stringa, lunghezza 9
4: {A,}            -> 2 stringhe, una vuota, lunghezza non vuota 1
5: {CTGT}          -> 1 stringa, lunghezza 4
```

Totali:

```text
numero loci: 5
numero stringhe/alternative: 8
stringhe non vuote: 7
stringhe vuote: 1
caratteri totali biologici non vuoti: 26
```

## Stato Attuale Del Programma

Attualmente `EDSBWTsearch` fa:

```text
indice EDS-BWT + query.eds
    -> array di similarita per locus
    -> similarita totale query.eds rispetto all'indice
```

Formalmente, se l'indice rappresenta `P2` e la query e `P1.eds`, il valore restituito e:

```text
A(P1, P2)
```

Il prossimo passo sara costruire un main/orchestratore che faccia:

```text
input: P1.eds, P2.eds

1. costruisci indice P1
2. costruisci indice P2
3. calcola A(P1,P2) cercando P1.eds nell'indice P2
4. calcola A(P2,P1) cercando P2.eds nell'indice P1
5. restituisci (A(P1,P2) + A(P2,P1)) / 2
```

## Note Aperte

Restano alcune scelte progettuali da decidere:

- se costruire il flusso finale con uno script shell oppure con un nuovo main C++;
- se salvare gli array per-locus in un file;
- se la media finale tra `A(P1,P2)` e `A(P2,P1)` deve essere aritmetica semplice o pesata;
- se le stringhe vuote debbano essere ignorate, come ora, oppure trattate con una regola esplicita diversa;
- se rendere silenzioso l'output di debug/stampa dei pattern per usare il programma in pipeline.


## Aggiornamento Finale: Programma Unico E Test Su Dataset Reale

In questa fase il progetto e stato portato dal prototipo basato su `EDSBWTsearch` a un primo programma unico per il calcolo della similarita simmetrica tra due EDS.

### Nuovo programma `EDSBWTsimilarity`

E stato aggiunto il file:

```text
EDSBWTsimilarity.cpp
```

Il nuovo main prende in input due file `.eds`:

```bash
./EDSBWTsimilarity P1.eds P2.eds
```

Il flusso implementato e:

```text
1. prende P1.eds e P2.eds;
2. costruisce l indice EDS-BWT di P1;
3. costruisce l indice EDS-BWT di P2;
4. carica l indice di P1 e l indice di P2;
5. calcola A(P1,P2), cercando i loci di P1 nell indice di P2;
6. calcola A(P2,P1), cercando i loci di P2 nell indice di P1;
7. restituisce la media aritmetica:

   Similarity(P1,P2) = (A(P1,P2) + A(P2,P1)) / 2
```

Il programma usa internamente i tool gia presenti nel progetto:

```text
eds_to_fasta
BCR_LCP_GSA/BCR_LCP_GSA
EOFpos_to_everything
EDSBWTsearch / classe EDSBWT
```

Quindi il flusso completo e diventato:

```text
P1.eds, P2.eds
   -> conversione EDS -> FASTA
   -> costruzione indici EDS-BWT
   -> ricerca direzionale sui loci
   -> A(P1,P2), A(P2,P1)
   -> similarita finale
```

### Modifica al `makefile` principale

Il `makefile` principale e stato esteso con un nuovo target:

```text
similarity
```

Questo target compila il binario:

```text
EDSBWTsimilarity
```

La compilazione usata nei test e stata:

```bash
make -B similarity RECOVERBW=0 \
  SDSL_INC=/home/giordanodb/Desktop/Programma/EDS-BWT/sdsl_lib/include \
  SDSL_LIB=/home/giordanodb/Desktop/Programma/EDS-BWT/sdsl_lib/lib
```

### Problema trovato con il dataset reale

Provando il confronto tra:

```text
../dataset/23A.eds
../dataset/23B.eds
```

la prima esecuzione si fermava durante la costruzione dell indice di `23A.eds`, dentro:

```text
BCR_LCP_GSA/BCR_LCP_GSA
```

Il crash non era nella parte di similarita, ma nella costruzione della BWT.

Dal FASTA generato da `23A.eds` e emerso:

```text
numero sequenze: 249
lunghezza minima: 1
lunghezza massima: 2187
caratteri totali: 35342
```

Il problema era che `BCR_LCP_GSA` era configurato con lunghezze di sequenza su 1 byte:

```cpp
#define dataTypeLengthSequences 0
```

cioe con tipo `unsigned char`, sufficiente solo per sequenze lunghe circa fino a 254 caratteri.

### Modifica a `BCR_LCP_GSA/Parameters.h`

Per permettere l uso del dataset reale, in:

```text
BCR_LCP_GSA/Parameters.h
```

e stato modificato:

```cpp
#define dataTypeLengthSequences 0
```

in:

```cpp
#define dataTypeLengthSequences 1
```

In questo modo `BCR_LCP_GSA` usa `unsigned short` per rappresentare le lunghezze delle sequenze, sufficiente per sequenze fino a circa 65535 caratteri.

Backup creato:

```text
BCR_LCP_GSA/Parameters.h.backup-lengthseq-20260603
```

### Modifica a `BCR_LCP_GSA/makefile`

Dopo la ricompilazione di `BCR_LCP_GSA`, si e visto che il programma non generava piu il file delle posizioni dei terminatori, necessario al flusso EDS-BWT.

Nel makefile di BCR:

```text
BCR_LCP_GSA/makefile
```

e stato quindi reso permanente:

```make
STORE_INDICES_DOLLARS = 1
```

Questo fa si che `BCR_LCP_GSA` venga compilato con:

```text
-DSTORE_ENDMARKER_POS=1
```

e produca le informazioni necessarie a `EOFpos_to_everything`.

Backup creato:

```text
BCR_LCP_GSA/makefile.backup-store-eofpos-20260603
```

### Ricompilazione di `BCR_LCP_GSA`

Dopo le modifiche, `BCR_LCP_GSA` e stato ricompilato con successo:

```bash
cd /home/giordanodb/Desktop/Programma/EDS-BWT/BCR_LCP_GSA
make clean
make BCR_BWTCollection STORE_INDICES_DOLLARS=1
```

Nel log del nuovo test si vede che ora BCR usa lunghezze a 2 byte:

```text
dataTypelenSeq: sizeof(type of seq length): 2 bytes
```

Per `23A.eds` la costruzione dell indice e arrivata correttamente a:

```text
The BWT et al. is ready!
The End!
```

Lo stesso e avvenuto per `23B.eds`.

### Test completo su dataset reale

Il test completo eseguito e stato:

```bash
./EDSBWTsimilarity ../dataset/23A.eds ../dataset/23B.eds
```

Il programma ha:

```text
1. costruito l indice EDS-BWT di 23A.eds;
2. costruito l indice EDS-BWT di 23B.eds;
3. calcolato A(23A,23B);
4. calcolato A(23B,23A);
5. calcolato la media finale.
```

Risultato finale ottenuto:

```text
A(P1,P2): 0.825671
A(P2,P1): 0.77625
Similarity: 0.80096
```

Quindi, con `P1 = 23A.eds` e `P2 = 23B.eds`, la similarita simmetrica calcolata dal programma e:

```text
0.80096
```

### File generati durante il test

Durante il test, nella cartella `../dataset` vengono generati file intermedi e file di indice con prefisso `23A` e `23B`, ad esempio:

```text
23A.ebwt
23A.bitvector
23A_info.aux
23A_bwt_*.aux
23A_bv_*.aux

23B.ebwt
23B.bitvector
23B_info.aux
23B_bwt_*.aux
23B_bv_*.aux
```

I file `.eds` originali non vengono modificati.

## Aggiunta funzione di pulizia temporanei






## Modifiche principali
- Cambio dell'input di EDSBWTSearch a Index - .eds
- Creazione parser EDS -> Vect Vect String
- Modifica di backwardSearch da match booleano a lunghezza
- computeSimilarityFromEDS (misura asym tra EDSBWTidx e EDS)


## Da aggiungere
Covid // Grafico
Batteri // Intra specie o Tra colori diversi
- 10 genomi per specie 