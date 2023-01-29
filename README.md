# full-br-index
br-index (https://github.com/U-Ar/br-index) also supporting pattern contractions.

## About

This repository provides the __fully functional__ bi-directional r-index (_full br-index_).


The r-index is the compressed text index which supports efficient count(P) and locate(P).
Its uses O(r) words of space, whose r is the number of equal-letter runs in BWT of the text.

The br-index is the bi-directional text index based on the r-index, supporting _left-extension_, _right-extension_, and _locate_ at any step of the search. It consists of the r-index on the text, the r-index on the reversed text, and the permuted LCP array (PLCP).

The full br-index in this repository is the enhanced version of the br-index with additional data structures.
In addition to extensions and locate, it supports _left-contraction_ and _right-contraction_, which are the inverse operations of _left-extension_ and _right-extension_. There are time-space trade-offs based on the value of the parameter *bl*. The recommended value is around 15, and it is set to 8 by default. The supports of the five key operations enable complex searches. As an example we implement the _Maximal Exact Matches_ query in _src/bri-query.cpp_.

The index is constructed by the _Prefix-Free Parsing_ method. Suitable for highly repetitive huge text collections.
The rather simple in-memory construction is also supported (-i option for _bri-build_), but it consumes memory.

## System Requirements

- This project is based on [sdsl-lite](https://github.com/simongog/sdsl-lite) library.
Install sdsl-lite beforehand and modify variables SDSL_INCLUDE and SDSL_LIB in _CMakeLists.txt_.

- This project has been tested under gcc 7.5.0, gcc 8.4.0 and gcc 12.1.0.

## How to Use

Firstly, clone the repository. Since a submodule is used ([iutest](https://github.com/srz-zumix/iutest)), recursive cloning is necessary.
```bash
git clone --recursive https://github.com/U-Ar/full-br-index.git
```
In order to build, execute following commands: (This project is using CMake)
```bash
mkdir build
cd build
cmake ..
make
```
You can also run unit tests and integration tests. 
```
make unit-tests
make integration-tests
```
For the detail of the tests, refer to _test/README.md_.

By default 5 entry-point executables will be created in the _build_ directory.
<dl>
	<dt>bri-build (Python script)</dt>
	<dd>Builds the br-index on the input text file using Prefix-Free Parsing. Using -t option is not recommended now: it causes errors during the computation of multi-threaded PFP.</dd>
    <dt>bri-query</dt>
	<dd>Computes searching queries on the index. (count, locate, MEMs, full-task)</dd>
	<dt>bri-space</dt>
	<dd>Shows the statistics of the text and the breakdown of the index space usage.</dd>
	<dt>unit-test</dt>
	<dd>runs unit tests.</dd>
	<dt>integration-test (Python script)</dt>
	<dd>runs integration tests.</dd>
</dl>

The codes used for experiments are archived in the _experiment\_archive_ directory. 

## Citations

For the br-index (without contractions), cite 
- Arakawa, Y., Navarro, G., & Sadakane, K. (2022). Bi-Directional r-Indexes. In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022). Schloss Dagstuhl-Leibniz-Zentrum für Informatik.

The fully-functional br-index: not published.
