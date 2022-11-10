# br-index tests

## About

br-index unit test suite.

## How to run

In order to build and run all tests, move to *build* directory and execute
```
cmake ..
make unit-tests
```

## List of tests

- SparseSdVectorTest
- HuffmanStringTest
- RleStringTest
- PermutedLcpTest
- BrIndexInmemoryTest

## Dataset used for integration tests

In addition to the dataset described below, you can add texts to *dataset* directory.
The br-index will be tested also on those texts.

Testing on large files (>1GB) is not recommended; It'll take too much time and run out of memory.

### einstein.txt

100000 bytes prefix of einstein.en.txt taken from (http://pizzachili.dcc.uchile.cl/repcorpus/real/einstein.en.txt.gz).

### influenza.txt

100000 bytes prefix of influenza DNA sequence taken from (http://pizzachili.dcc.uchile.cl/repcorpus/real/influenza.gz).

### dna.txt

100000 bytes prefix of dna.50MB taken from (http://pizzachili.dcc.uchile.cl/texts/dna/dna.50MB.gz).