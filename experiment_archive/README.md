# Archive of experimental codes for paper

## About

Codes for experiments (little messed) are put here.
You can build this directory independently.

## Brief description

- _src/br_index_fixed.hpp_ is used to measure index sizes of the br-index and the modified r-index with varying *bl*.
- _src/br_index_st.hpp_ is the full br-index implementation used to compute MEMs and suffix tree operations.
- _src/window_slide.cpp_ is the code to measure the time to slide fixed length window on a given text with various methods and *bl*.
- _src/k_vs_space.cpp_ is the code to measure the index size with varying parameter *bl*.
- _src/bri-mems.cpp_ is the code to evaluate Maximal Exact Matches.
- _src/x_treeop.cpp_ is the code to evaluate basic suffix tree operations.

_src/k_vs_time.cpp_, _src/fixedtest.cpp_, and _src/br_index_limited.hpp_ were not used.
