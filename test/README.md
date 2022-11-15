# br-index tests

## About

br-index test suite. 

Unit tests and integration tests are implemented.

## How to run

In order to build and run unit tests & integration tests, move to *build* directory and execute
```
cmake ..
make unit-tests
make integration-tests
```

## List of unit tests

- SparseSdVectorTest
- HuffmanStringTest
- RleStringTest
- PermutedLcpTest
- BrIndexInmemoryTest

## Detailed description of integration tests

### Pre-test

Both of in-memory construction and PFP construction are conducted on each text with *bl*=64. 
They are loaded with parameter *bl*=0 and *bl*=64.

In total, 4 different indexes are tested for each text:
- In-memory build, bl=0
- In-memory build, bl=64
- PFP build, bl=0
- PFP build, bl=64

### Test1: br-index functionality test

We divide the target text into substrings of the given length.
For each substring, we

- sequentially apply *left-extension* from an empty string to the entire substring. After that, we
    - sequentially apply *left-contraction* from the entire substring to an empty string, and
    - sequentially apply *right-contraction* from the entire substring to an empty string. Also we
- sequentially apply *right-extension* from an empty string to the entire substring. 

For each step of extension/contraction, we
- compare the returned value of *count* with Karp-Rabin hashtable, and
- confirm the returned list of *locate* correctly points to the starting positions of the currently searched pattern.

### Test2: computing Maximal Exact Matches

We mutate characters of the target text with the given probability (0.05 by default). We let this mutated text be the pattern text.

We compute Maximal Exact Matches on the pattern and the target.

### Datasets used

In addition to the dataset described below, you can add text files to *dataset* directory (.txt extension is necessary).
The br-index will be tested also on those texts.

Testing on large files (>1GB) is not recommended; It'll take too much time and run out of memory.

#### einstein.txt

100000 bytes prefix of einstein.en.txt taken from (http://pizzachili.dcc.uchile.cl/repcorpus/real/einstein.en.txt.gz).

#### influenza.txt

100000 bytes prefix of influenza DNA sequence taken from (http://pizzachili.dcc.uchile.cl/repcorpus/real/influenza.gz).

#### dna.txt

100000 bytes prefix of dna.50MB taken from (http://pizzachili.dcc.uchile.cl/texts/dna/dna.50MB.gz).